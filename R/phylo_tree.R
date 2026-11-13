# Tree handling for the phylogenetic module. Method-agnostic: used by every HyPhy method.
#
# Branches are identified by their bipartition of the taxa present in an orthogroup, never
# by node names, so a hypothesis keeps one meaning after pruning. Rules are documented in
# docs/phylo-methods.md and docs/dev/phylo-design.md.

.phylo_name_pattern <- "^[A-Za-z][A-Za-z0-9_]*$"

#' Stop unless ape is installed
#' @keywords internal
.need_ape <- function() {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("The phylogenetic commands need the R package 'ape'. Install it with install.packages('ape'), ",
         "or run inside the dndsR container.", call. = FALSE)
  }
}

#' Read and check a species tree
#'
#' @param path Newick file with one tree.
#' @return An `ape` phylo object, with branch lengths and node labels removed.
#' @keywords internal
.phylo_read_species_tree <- function(path) {
  .need_ape()
  if (!file.exists(path)) stop("Species tree not found: ", path, call. = FALSE)
  tree <- tryCatch(ape::read.tree(path), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(tree) || inherits(tree, "multiPhylo")) {
    stop("Could not read one Newick tree from ", path, call. = FALSE)
  }
  .phylo_check_tip_names(tree$tip.label, path)
  if (ape::Ntip(tree) < 3L) stop("The species tree needs at least 3 taxa: ", path, call. = FALSE)
  tree$edge.length <- NULL
  labels <- tree$node.label
  usable <- !is.null(labels) && all(grepl(.phylo_name_pattern, labels)) && !anyDuplicated(labels)
  tree$node.label <- NULL
  attr(tree, "node_names") <- if (usable) labels else NULL
  tree
}

#' Readable names for the internal nodes of a species tree
#'
#' The tree's own node labels if it has valid, unique ones; otherwise `mrca(x,y)`, with
#' the first tip under each of the node's first two children.
#' @keywords internal
.phylo_internal_names <- function(tree) {
  given <- attr(tree, "node_names")
  if (!is.null(given)) return(given)
  n <- ape::Ntip(tree)
  first_tip <- function(v) {
    while (v > n) v <- tree$edge[tree$edge[, 1] == v, 2][1]
    tree$tip.label[v]
  }
  vapply(seq_len(tree$Nnode) + n, function(v) {
    kids <- tree$edge[tree$edge[, 1] == v, 2]
    sprintf("mrca(%s,%s)", first_tip(kids[1]), first_tip(kids[2]))
  }, character(1))
}

#' @keywords internal
.phylo_check_tip_names <- function(tips, source) {
  bad <- tips[!grepl(.phylo_name_pattern, tips)]
  if (length(bad)) {
    stop("Taxon names must start with a letter and contain only letters, digits, and underscores ",
         "(HyPhy renames other characters). In ", source, ": ", paste(utils::head(bad, 5), collapse = ", "),
         call. = FALSE)
  }
  dup <- unique(tips[duplicated(tips)])
  if (length(dup)) stop("Duplicate taxon names in ", source, ": ", paste(dup, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Descendant tips of every node, as tip labels
#' @keywords internal
.phylo_descendants <- function(tree) {
  n <- ape::Ntip(tree)
  parts <- ape::prop.part(tree)
  c(as.list(tree$tip.label), lapply(parts, function(p) tree$tip.label[p]))
}

#' Canonical key of a bipartition of the present taxa
#'
#' The side written is the one without the first present taxon in C-locale order, so the
#' same split gets the same key from any tree or rooting.
#' @param side Taxa on one side.
#' @param present All taxa present.
#' @return A string, or `NA` if the split is trivial (empty or all taxa).
#' @keywords internal
.phylo_split_key <- function(side, present) {
  side <- intersect(side, present)
  if (!length(side) || length(side) == length(present)) return(NA_character_)
  first <- sort(present, method = "radix")[1L]
  if (first %in% side) side <- setdiff(present, side)
  paste(sort(side, method = "radix"), collapse = ",")
}

#' Readable side of a split: the smaller side, or the canonical side when both are equal
#' @keywords internal
.phylo_split_label <- function(side, present) {
  side <- intersect(side, present)
  other <- setdiff(present, side)
  if (length(side) == length(other)) return(.phylo_split_key(side, present))
  paste(sort(if (length(side) < length(other)) side else other, method = "radix"), collapse = ",")
}

#' Read a foreground file and resolve it on the species tree
#'
#' The file is tab-separated with columns `scope` and `taxa` (comma-separated). Scopes:
#' `terminal` (the terminal branch of each listed taxon), `stem` (the branch leading to the
#' most recent common ancestor of the listed taxa, which must form a clade), and `clade`
#' (that stem and every branch below it). Rows are combined.
#'
#' @param path Foreground file.
#' @param tree Species tree from [.phylo_read_species_tree()].
#' @return A list describing the hypothesis on the reference tree.
#' @keywords internal
.phylo_read_foreground <- function(path, tree) {
  if (!file.exists(path)) stop("Foreground file not found: ", path, call. = FALSE)
  spec <- utils::read.delim(path, colClasses = "character", strip.white = TRUE, comment.char = "#")
  names(spec) <- tolower(names(spec))
  if (!all(c("scope", "taxa") %in% names(spec)) || !nrow(spec)) {
    stop("The foreground file needs a header and columns 'scope' and 'taxa': ", path, call. = FALSE)
  }
  spec$scope <- tolower(spec$scope)
  bad_scope <- setdiff(spec$scope, c("terminal", "stem", "clade"))
  if (length(bad_scope)) {
    stop("Unknown foreground scope: ", paste(bad_scope, collapse = ", "), ". Use terminal, stem, or clade.", call. = FALSE)
  }
  desc <- .phylo_descendants(tree)
  n_nodes <- length(desc)
  root <- ape::Ntip(tree) + 1L
  child_of_edge <- tree$edge[, 2]
  node_rows <- vector("list", n_nodes) # hypothesis rows asking for the edge above each node
  node_req <- logical(n_nodes)         # must also exist in a gene tree: terminal and stem branches
  clades <- list()
  if (any(spec$scope %in% c("stem", "clade")) && !ape::is.rooted(tree)) {
    stop("Stem and clade hypotheses need a rooted species tree, because a clade is only defined with a root.", call. = FALSE)
  }
  ask <- function(nodes, i) for (v in nodes) node_rows[[v]] <<- c(node_rows[[v]], i)
  for (i in seq_len(nrow(spec))) {
    taxa <- unique(trimws(strsplit(spec$taxa[i], ",", fixed = TRUE)[[1]]))
    taxa <- taxa[nzchar(taxa)]
    missing <- setdiff(taxa, tree$tip.label)
    if (!length(taxa) || length(missing)) {
      stop("Foreground row ", i, " names taxa not in the species tree: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    if (spec$scope[i] == "terminal") {
      ask(match(taxa, tree$tip.label), i)
      node_req[match(taxa, tree$tip.label)] <- TRUE
      next
    }
    if (length(taxa) < 2L) {
      stop("Foreground row ", i, ": a ", spec$scope[i], " needs at least two taxa. Use scope 'terminal' for one taxon.", call. = FALSE)
    }
    mrca <- ape::getMRCA(tree, taxa)
    extra <- setdiff(desc[[mrca]], taxa)
    if (length(extra)) {
      stop("Foreground row ", i, ": ", paste(taxa, collapse = ", "), " do not form a clade in the species tree. ",
           "Their smallest clade also contains ", paste(extra, collapse = ", "), ".", call. = FALSE)
    }
    if (mrca == root) stop("Foreground row ", i, ": the root has no stem branch.", call. = FALSE)
    ask(mrca, i)
    node_req[mrca] <- TRUE
    if (spec$scope[i] == "clade") {
      ask(setdiff(which(vapply(desc, function(d) all(d %in% taxa), logical(1))), root), i)
      clades[[length(clades) + 1L]] <- taxa
    }
  }
  edges <- data.frame(
    child = child_of_edge,
    name = c(tree$tip.label, .phylo_internal_names(tree))[child_of_edge],
    foreground = lengths(node_rows)[child_of_edge] > 0L,
    required = node_req[child_of_edge],
    stringsAsFactors = FALSE
  )
  edges$rows <- node_rows[child_of_edge]
  edges$descendants <- desc[child_of_edge]
  list(spec = spec, tree = tree, edges = edges, clades = clades,
       taxa = tree$tip.label, foreground_taxa = sort(unique(unlist(edges$descendants[edges$foreground]))))
}

#' Basis of each exclusion or failure reason
#'
#' `mathematical`: the requested branches cannot be separated, or do not exist, in the
#' topology for this orthogroup. `no_data`: no foreground taxon is present. `policy`:
#' testable in principle, but excluded by a conservative rule (single-copy only, minimum
#' taxa, no narrowing of the hypothesis). `data_quality`: the sequences fail checks.
#' `input`: the inputs do not match. `software`: a step failed.
#' @keywords internal
.phylo_reason_basis <- function(reason) {
  map <- c(
    foreground_unidentifiable = "mathematical", foreground_not_in_gene_tree = "mathematical",
    foreground_absent = "no_data",
    foreground_partly_unidentifiable = "policy", foreground_incomplete = "policy", foreground_branch_zero_length = "policy",
    multi_copy = "policy", too_few_taxa = "policy",
    cds_invalid = "data_quality", internal_stop = "data_quality", split_codon_gap = "data_quality",
    empty_sequence = "data_quality", terminal_stop = "data_quality", length_not_multiple_of_3 = "data_quality",
    invalid_character = "data_quality",
    gene_not_found = "input", unknown_taxon = "input", empty_alignment = "input", not_aligned = "input",
    gene_tree_missing = "input", gene_tree_unreadable = "input", gene_tree_taxa_mismatch = "input",
    hyphy_error = "software", hyphy_output_invalid = "software", label_mismatch = "software",
    alignment_failed = "software", backtranslation_mismatch = "software", prepared_tree_mismatch = "software"
  )
  out <- unname(map[reason])
  out[is.na(reason)] <- NA_character_
  out
}

#' Splits of an unrooted topology, with adjacency to labelled branches
#' @keywords internal
.phylo_topology <- function(tree, present) {
  if (ape::is.rooted(tree) && ape::Ntip(tree) > 2L) tree <- ape::unroot(tree)
  desc <- .phylo_descendants(tree)
  child <- tree$edge[, 2]
  list(tree = tree, desc = desc, child = child, parent = tree$edge[, 1],
       keys = vapply(desc[child], .phylo_split_key, character(1), present = present))
}

#' Map a foreground hypothesis onto the tree for one orthogroup
#'
#' Branches of the species tree that join after pruning become one branch. A requested
#' branch joined to a background branch cannot be separated from it. If no requested branch
#' can be separated, the hypothesis is `foreground_unidentifiable`; if only some can, testing
#' those would narrow the hypothesis, so it is `foreground_partly_unidentifiable`. If a
#' hypothesis row has no taxa present, or a terminal row lacks a listed taxon, it is
#' `foreground_incomplete`. Missing members of a clade-scope row are allowed and counted.
#' Nothing is relabelled.
#'
#' With a gene tree, discordance with the pruned species tree is measured in both modes.
#' In `"topology"` mode the gene tree is given to HyPhy: every terminal and stem foreground
#' branch, and each clade's stem, must exist in it (`foreground_not_in_gene_tree`
#' otherwise), and branches inside a clade-scope foreground are foreground. In `"report"`
#' mode the species tree is used, and the gene tree is only compared.
#'
#' @param hyp Result of [.phylo_read_foreground()].
#' @param present Taxa present in the orthogroup.
#' @param gene_tree Optional `ape` phylo with tips named by taxon.
#' @param gene_tree_mode `"topology"` or `"report"`.
#' @return A list with `ok`, `reason`, `detail`, `newick` (labelled, unrooted, no branch
#'   lengths), `branches` (one row per branch), and `info` (foreground and discordance
#'   measures).
#' @keywords internal
.phylo_map_foreground <- function(hyp, present, gene_tree = NULL, gene_tree_mode = c("topology", "report")) {
  gene_tree_mode <- match.arg(gene_tree_mode)
  present <- sort(intersect(hyp$taxa, present), method = "radix")
  e <- hyp$edges
  e$key <- vapply(e$descendants, .phylo_split_key, character(1), present = present)
  fg_ref <- e[e$foreground, , drop = FALSE]
  info <- list(fg_requested_branches = nrow(fg_ref), fg_absent_branches = sum(is.na(fg_ref$key)))
  fail <- function(reason, detail) {
    list(ok = FALSE, reason = reason, detail = detail, newick = NA_character_, branches = NULL, info = info)
  }
  e <- e[!is.na(e$key), , drop = FALSE]
  groups <- split(e, e$key)
  status <- vapply(groups, function(g) if (all(g$foreground)) "fg" else if (any(g$foreground)) "mixed" else "bg", character(1))
  merged <- vapply(groups, function(g) paste(g$name, collapse = "+"), character(1))

  if (!any(status %in% c("fg", "mixed"))) {
    return(fail("foreground_absent", "no foreground taxon is present"))
  }
  if (any(status == "mixed")) {
    k <- names(status)[status == "mixed"]
    joined <- paste0("these foreground and background branches of the species tree become one branch: ", paste(merged[k], collapse = "; "))
    if (!any(status == "fg")) {
      return(fail("foreground_unidentifiable", paste0("no requested branch can be separated from background history; ", joined)))
    }
    return(fail("foreground_partly_unidentifiable", paste0(
      "some requested branches remain separable, but testing only them would narrow the hypothesis; ", joined)))
  }
  for (i in seq_len(nrow(hyp$spec))) {
    mine <- fg_ref[vapply(fg_ref$rows, function(r) i %in% r, logical(1)), , drop = FALSE]
    gone <- is.na(mine$key)
    if (all(gone)) {
      return(fail("foreground_incomplete", sprintf("no taxon of foreground row %d (%s %s) is present", i, hyp$spec$scope[i], hyp$spec$taxa[i])))
    }
    if (hyp$spec$scope[i] == "terminal" && any(gone)) {
      return(fail("foreground_incomplete", sprintf("terminal foreground row %d lacks: %s", i, paste(mine$name[gone], collapse = ", "))))
    }
  }
  fg_keys <- names(status)[status == "fg"]
  required_keys <- names(groups)[vapply(groups, function(g) any(g$required), logical(1))]

  species <- .phylo_topology(ape::drop.tip(hyp$tree, setdiff(hyp$taxa, present), collapse.singles = TRUE), present)
  if (!setequal(species$keys, names(status))) {
    stop("Internal error: pruned species tree branches do not match the reference mapping.", call. = FALSE)
  }
  lab <- species
  source <- "species_tree"
  if (!is.null(gene_tree)) {
    gene_tree$edge.length <- NULL
    gene_tree$node.label <- NULL
    gene <- .phylo_topology(gene_tree, present)
    missing_species <- setdiff(species$keys, gene$keys)
    fg_edge <- species$keys %in% fg_keys
    touched <- unique(c(species$parent[fg_edge], species$child[fg_edge]))
    near <- species$keys[species$parent %in% touched | species$child %in% touched]
    info <- c(info, list(
      gene_tree_mode = gene_tree_mode,
      gt_missing_species_splits = length(missing_species),
      gt_extra_splits = length(setdiff(gene$keys, species$keys)),
      gt_discordant_near_foreground = any(missing_species %in% near),
      gt_missing_splits = paste(vapply(missing_species, .phylo_split_label, character(1), present = present), collapse = ";")
    ))
    if (gene_tree_mode == "topology") {
      absent <- setdiff(required_keys, gene$keys)
      if (length(absent)) {
        return(fail("foreground_not_in_gene_tree", paste0(
          "the gene tree has no branch splitting off these requested taxa: ",
          paste(vapply(absent, .phylo_split_label, character(1), present = present), collapse = "; "))))
      }
      lab <- gene
      source <- "gene_tree"
    }
  }
  fg <- lab$keys %in% fg_keys
  if (source == "gene_tree") {
    for (clade in hyp$clades) {
      k <- intersect(clade, present)
      inside <- vapply(lab$desc[lab$child], function(d) {
        d <- intersect(d, present)
        all(d %in% k) || all(setdiff(present, d) %in% k)
      }, logical(1))
      fg <- fg | inside
    }
  }

  n_tip <- ape::Ntip(lab$tree)
  node_names <- c(lab$tree$tip.label, paste0("Node", seq_len(lab$tree$Nnode)))
  node_fg <- logical(length(node_names))
  node_fg[lab$child] <- fg
  newick <- .phylo_write_newick(lab$tree, node_names, node_fg, "FG")
  branches <- data.frame(
    branch = node_names[lab$child],
    terminal = lab$child <= n_tip,
    taxa_side = vapply(lab$desc[lab$child], .phylo_split_label, character(1), present = present),
    foreground = fg,
    species_tree_branches = ifelse(lab$keys %in% names(merged), merged[lab$keys], NA_character_),
    in_species_tree = lab$keys %in% species$keys,
    topology_source = source,
    stringsAsFactors = FALSE
  )
  list(ok = TRUE, reason = NA_character_, detail = NA_character_, newick = newick,
       branches = branches[order(!branches$foreground, branches$branch), , drop = FALSE], info = info)
}

#' Write a topology in Newick with HyPhy branch labels
#'
#' @param tree `ape` phylo.
#' @param node_names Names for all nodes (tips first).
#' @param node_fg Whether the branch above each node is labelled.
#' @param label HyPhy label written in braces.
#' @return A Newick string without branch lengths.
#' @keywords internal
.phylo_write_newick <- function(tree, node_names, node_fg, label = "FG") {
  root <- ape::Ntip(tree) + 1L
  kids <- split(tree$edge[, 2], factor(tree$edge[, 1], levels = seq_along(node_names)))
  write_node <- function(v) {
    tag <- if (node_fg[v]) paste0("{", label, "}") else ""
    if (v <= ape::Ntip(tree)) return(paste0(node_names[v], tag))
    inner <- paste(vapply(kids[[v]], write_node, character(1)), collapse = ",")
    if (v == root) paste0("(", inner, ");") else paste0("(", inner, ")", node_names[v], tag)
  }
  write_node(root)
}

#' Read an optional gene tree for one orthogroup
#'
#' Tips may be named by taxon, by sequence name (`taxon|gene`), by gene ID, or as
#' OrthoFinder names them, `<proteome>_<gene>`, where the proteome name is the taxon or an
#' alias given in the genome table.
#'
#' @param path Newick file, or `NULL`.
#' @param seqs Data frame with `seq_name`, `taxon`, and `gene_id`.
#' @param aliases Named character vector: alias to taxon.
#' @return An `ape` phylo with tips renamed to taxa, or `NULL`. Stops with a reason as the
#'   condition message when the tree does not match the sequences.
#' @keywords internal
.phylo_read_gene_tree <- function(path, seqs, aliases = NULL) {
  if (is.null(path)) return(NULL)
  tree <- tryCatch(ape::read.tree(path), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(tree) || inherits(tree, "multiPhylo")) stop("gene_tree_unreadable")
  lookup <- c(stats::setNames(seqs$taxon, seqs$seq_name), stats::setNames(seqs$taxon, seqs$taxon))
  gene <- !is.na(seqs$gene_id)
  if (any(gene)) {
    lookup <- c(lookup, stats::setNames(seqs$taxon[gene], seqs$gene_id[gene]),
                stats::setNames(seqs$taxon[gene], paste0(seqs$taxon[gene], "_", seqs$gene_id[gene])))
    for (a in names(aliases)) {
      hit <- gene & seqs$taxon == aliases[[a]]
      lookup <- c(lookup, stats::setNames(seqs$taxon[hit], paste0(a, "_", seqs$gene_id[hit])))
    }
  }
  lookup <- lookup[!duplicated(names(lookup))]
  mapped <- unname(lookup[tree$tip.label])
  if (anyNA(mapped) || anyDuplicated(mapped) || !setequal(mapped, seqs$taxon)) stop("gene_tree_taxa_mismatch")
  tree$tip.label <- mapped
  tree
}
