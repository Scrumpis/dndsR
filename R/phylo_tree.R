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
  node_fg <- logical(n_nodes) # foreground status of the edge above each node
  node_req <- logical(n_nodes) # must also exist in a gene tree: terminal and stem branches
  clades <- list()
  needs_root <- any(spec$scope %in% c("stem", "clade"))
  if (needs_root && !ape::is.rooted(tree)) {
    stop("Stem and clade hypotheses need a rooted species tree, because a clade is only defined with a root.", call. = FALSE)
  }
  for (i in seq_len(nrow(spec))) {
    taxa <- unique(trimws(strsplit(spec$taxa[i], ",", fixed = TRUE)[[1]]))
    taxa <- taxa[nzchar(taxa)]
    missing <- setdiff(taxa, tree$tip.label)
    if (!length(taxa) || length(missing)) {
      stop("Foreground row ", i, " names taxa not in the species tree: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    if (spec$scope[i] == "terminal") {
      node_fg[match(taxa, tree$tip.label)] <- TRUE
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
    node_fg[mrca] <- TRUE
    node_req[mrca] <- TRUE
    if (spec$scope[i] == "clade") {
      below <- which(vapply(desc, function(d) all(d %in% taxa), logical(1)))
      node_fg[below] <- TRUE
      clades[[length(clades) + 1L]] <- taxa
    }
  }
  edges <- data.frame(
    child = child_of_edge,
    name = c(tree$tip.label, .phylo_internal_names(tree))[child_of_edge],
    foreground = node_fg[child_of_edge],
    required = node_req[child_of_edge],
    stringsAsFactors = FALSE
  )
  edges$descendants <- desc[child_of_edge]
  list(spec = spec, tree = tree, edges = edges, clades = clades,
       taxa = tree$tip.label, foreground_taxa = sort(unique(unlist(edges$descendants[edges$foreground]))))
}

#' Map a foreground hypothesis onto the tree for one orthogroup
#'
#' Branches of the species tree that join after pruning become one branch. If such a join
#' combines foreground and background branches, the hypothesis is not identifiable for this
#' orthogroup and nothing is tested. With a gene tree, every foreground branch must exist in
#' it; branches inside a clade-scope foreground are foreground.
#'
#' @param hyp Result of [.phylo_read_foreground()].
#' @param present Taxa present in the orthogroup.
#' @param gene_tree Optional `ape` phylo with tips named by taxon.
#' @return A list with `ok`, `reason`, `detail`, `newick` (labelled, unrooted, no branch
#'   lengths), and `branches` (one row per branch).
#' @keywords internal
.phylo_map_foreground <- function(hyp, present, gene_tree = NULL) {
  fail <- function(reason, detail) list(ok = FALSE, reason = reason, detail = detail, newick = NA_character_, branches = NULL)
  present <- sort(intersect(hyp$taxa, present), method = "radix")
  e <- hyp$edges
  e$key <- vapply(e$descendants, .phylo_split_key, character(1), present = present)
  e <- e[!is.na(e$key), , drop = FALSE]
  groups <- split(e, e$key)
  status <- vapply(groups, function(g) if (all(g$foreground)) "fg" else if (any(g$foreground)) "mixed" else "bg", character(1))
  merged <- vapply(groups, function(g) paste(g$name, collapse = "+"), character(1))
  if (any(status == "mixed")) {
    k <- names(status)[status == "mixed"]
    return(fail("foreground_not_identifiable", paste0(
      "with the taxa present, these foreground and background branches of the species tree become one branch: ",
      paste(merged[k], collapse = "; "))))
  }
  if (!any(status == "fg")) {
    return(fail("foreground_absent", "no foreground branch remains with the taxa present"))
  }
  fg_keys <- names(status)[status == "fg"]
  required_keys <- names(groups)[vapply(groups, function(g) any(g$required), logical(1))]

  if (is.null(gene_tree)) {
    lab <- ape::drop.tip(hyp$tree, setdiff(hyp$taxa, present), collapse.singles = TRUE)
    source <- "species_tree"
  } else {
    lab <- gene_tree
    lab$edge.length <- NULL
    lab$node.label <- NULL
    source <- "gene_tree"
  }
  if (ape::is.rooted(lab) && ape::Ntip(lab) > 2L) lab <- ape::unroot(lab)
  desc <- .phylo_descendants(lab)
  child <- lab$edge[, 2]
  keys <- vapply(desc[child], .phylo_split_key, character(1), present = present)
  fg <- keys %in% fg_keys
  if (!is.null(gene_tree)) {
    absent <- setdiff(required_keys, keys)
    if (length(absent)) {
      return(fail("foreground_discordant", paste0(
        "the gene tree lacks these foreground branches of the species tree (taxa on one side): ",
        paste(absent, collapse = "; "))))
    }
    for (clade in hyp$clades) {
      k <- intersect(clade, present)
      inside <- vapply(desc[child], function(d) {
        d <- intersect(d, present)
        all(d %in% k) || all(setdiff(present, d) %in% k)
      }, logical(1))
      fg <- fg | inside
    }
  } else if (!setequal(keys, names(status))) {
    stop("Internal error: pruned species tree branches do not match the reference mapping.", call. = FALSE)
  }

  n_tip <- ape::Ntip(lab)
  node_names <- c(lab$tip.label, paste0("Node", seq_len(lab$Nnode)))
  node_fg <- logical(length(node_names))
  node_fg[child] <- fg
  newick <- .phylo_write_newick(lab, node_names, node_fg, "FG")
  branches <- data.frame(
    branch = node_names[child],
    terminal = child <= n_tip,
    taxa_side = vapply(desc[child], .phylo_split_label, character(1), present = present),
    foreground = fg,
    species_tree_branches = ifelse(keys %in% names(merged), merged[keys], NA_character_),
    topology_source = source,
    stringsAsFactors = FALSE
  )
  list(ok = TRUE, reason = NA_character_, detail = NA_character_, newick = newick,
       branches = branches[order(!branches$foreground, branches$branch), , drop = FALSE])
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
#' @param path Newick file, or `NULL`.
#' @param seq_species Named vector mapping sequence names to taxa.
#' @return An `ape` phylo with tips renamed to taxa, or `NULL`. Stops with a reason string
#'   attached as the condition message when the tree does not match the alignment.
#' @keywords internal
.phylo_read_gene_tree <- function(path, seq_species) {
  if (is.null(path)) return(NULL)
  tree <- tryCatch(ape::read.tree(path), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(tree) || inherits(tree, "multiPhylo")) stop("gene_tree_unreadable")
  tips <- tree$tip.label
  mapped <- ifelse(tips %in% names(seq_species), seq_species[tips], tips)
  if (anyDuplicated(mapped) || !setequal(mapped, unname(seq_species))) stop("gene_tree_taxa_mismatch")
  tree$tip.label <- unname(mapped)
  tree
}
