# Simulator for the lineage-level calibration (docs/dev/phylo-lineage-calibration.md).
# Builds on phylo-sim.R (codon model, F3x4 frequencies and kappa from params.json), adding
# omega classes for every branch. Independent of dndsR's package code.
#
#   Rscript --vanilla phylo-lineage-sim.R <params.json> <cell> <tree> <n genes> <seed> <out dir>
# Writes <out>/alignments/<gene>.fasta, species_tree.nwk (internal nodes labelled n1, n2,
# ..., which dndsR uses as branch names), foreground.tsv (the clade of the first two taxa),
# genes.tsv, and branch_truth.tsv (true omega classes of every branch).
here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(f) else "tools/dev/validation"
})
source(file.path(here, "phylo-sim.R")) # its own entry point runs only with 5 arguments; this script takes 6

balanced <- function(n, term, inner) {
  build <- function(tips) {
    if (length(tips) == 1L) return(sprintf("%s:%s", tips, term))
    h <- length(tips) / 2
    sprintf("(%s,%s):%s", build(tips[1:h]), build(tips[(h + 1):length(tips)]), inner)
  }
  tips <- paste0("t", seq_len(n))
  h <- n / 2
  sprintf("(%s,%s);", build(tips[1:h]), build(tips[(h + 1):n]))
}

lineage_trees <- function(par) list(
  T5 = par[["tree"]],
  T8b = balanced(8, 0.03, 0.02),
  T8p = "(((((((t1:0.005,t2:0.01):0.02,t3:0.03):0.02,t4:0.06):0.02,t5:0.1):0.02,t6:0.15):0.02,t7:0.25):0.02,t8:0.03);",
  T16 = balanced(16, 0.03, 0.015),
  T16m = balanced(16, 0.03, 0.015)
)

#' Omega classes and weights of every branch, by the node below it
cell_omegas <- function(cell, tree, fg_nodes) {
  n_nodes <- ape::Ntip(tree) + tree$Nnode
  one <- function(w) rep(list(w), n_nodes)
  prop3 <- c(0.6, 0.3, 0.1)
  third_tip <- match(tree$tip.label[3], tree$tip.label)
  out <- switch(cell,
    A0 = list(om = one(0.2), prop = 1),
    A1 = { o <- one(0.2); o[fg_nodes] <- list(0.6); list(om = o, prop = 1) },
    A2 = {
      old <- if (exists(".Random.seed", globalenv())) .Random.seed else NULL
      set.seed(20261120) # one fixed set of branch omegas per tree, not per gene
      w <- exp(log(0.2) + 0.7 * stats::rnorm(n_nodes))
      if (!is.null(old)) assign(".Random.seed", old, globalenv())
      list(om = as.list(w), prop = 1)
    },
    A3 = list(om = one(c(0.05, 0.3, 1)), prop = prop3),
    A4 = { o <- one(0.2); o[fg_nodes] <- list(1.5); list(om = o, prop = 1) },
    B0 = list(om = one(c(0.05, 0.3, 1)), prop = prop3),
    B1 = { o <- one(c(0.05, 0.3, 1)); o[third_tip] <- list(c(0.05, 0.3, 10)); list(om = o, prop = prop3) },
    C0 = list(om = one(c(0.05, 0.3, 1.5)), prop = prop3),
    C1 = { o <- one(c(0.05, 0.3, 1.5)); o[fg_nodes] <- list(c(0.05, 0.3, 1.5)^0.5); list(om = o, prop = prop3) },
    C2 = { o <- one(c(0.05, 0.3, 1.5)); o[fg_nodes] <- list(c(0.05, 0.3, 1.5)^2); list(om = o, prop = prop3) },
    # Follow-up study (docs/dev/phylo-lineage-followup.md).
    D0 = { o <- one(0.2); o[fg_nodes] <- list(1); list(om = o, prop = 1) },
    C0i = list(om = one(c(0.05, 0.3, 1.5)), prop = prop3, per_branch = TRUE),
    stop("Unknown cell: ", cell)
  )
  if (is.null(out[["per_branch"]])) out[["per_branch"]] <- FALSE
  out
}

#' One gene: every branch with its own omega classes; no synonymous rate variation
#'
#' A site's omega class is drawn once and kept on every branch, unless `per_branch`: then it
#' is drawn afresh on each branch, as HyPhy's random-effects branch-site models assume.
sim_gene_nodes <- function(tree, n_sites, pi, kappa, om, prop, per_branch = FALSE) {
  n_tip <- ape::Ntip(tree); root <- n_tip + 1L
  cls <- sample.int(length(prop), n_sites, TRUE, prob = prop)
  states <- vector("list", max(tree$edge))
  states[[root]] <- sample.int(n_sense, n_sites, TRUE, prob = pi)
  tr <- ape::reorder.phylo(tree, "cladewise")
  for (e in seq_len(nrow(tr$edge))) {
    par <- tr$edge[e, 1]; ch <- tr$edge[e, 2]
    out <- integer(n_sites)
    if (per_branch) cls <- sample.int(length(prop), n_sites, TRUE, prob = prop)
    for (k in unique(cls)) {
      idx <- which(cls == k)
      p <- p_matrix(kappa, om[[ch]][k], 1, 0, pi, tr$edge.length[e])
      from <- states[[par]][idx]
      for (a in unique(from)) { w <- idx[from == a]; out[w] <- sample.int(n_sense, length(w), TRUE, prob = p[a, ]) }
    }
    states[[ch]] <- out
  }
  aln <- t(vapply(seq_len(n_tip), function(i) sense[states[[i]]], character(n_sites)))
  rownames(aln) <- tree$tip.label
  aln
}

run_lineage_cell <- function(params_file, cell, tree_name, n, seed, out) {
  par <- jsonlite::read_json(params_file, simplifyVector = TRUE)
  RNGversion("4.0.0"); set.seed(seed)
  pi <- f3x4(matrix(unlist(par[["nuc_f3x4"]]), 3, byrow = TRUE))
  tree <- ape::read.tree(text = lineage_trees(par)[[tree_name]])
  tree$edge.length <- tree$edge.length * par[["c"]]
  tree$node.label <- paste0("n", seq_len(tree$Nnode))
  fg_tips <- tree$tip.label[1:2]
  fg_nodes <- c(match(fg_tips, tree$tip.label), ape::getMRCA(tree, fg_tips))
  spec <- cell_omegas(cell, tree, fg_nodes)
  dir.create(file.path(out, "alignments"), recursive = TRUE, showWarnings = FALSE)
  topo <- tree; topo$edge.length <- NULL
  ape::write.tree(topo, file.path(out, "species_tree.nwk"))
  writeLines(c("scope\ttaxa", paste0("clade\t", paste(fg_tips, collapse = ","))), file.path(out, "foreground.tsv"))
  names_all <- c(tree$tip.label, tree$node.label)
  child <- tree$edge[, 2]
  truth <- data.frame(branch_id = names_all[child], length_codon = tree$edge.length,
                      omega_classes = vapply(spec[["om"]][child], function(w) paste(signif(w, 4), collapse = ","), character(1)),
                      omega_mean = vapply(spec[["om"]][child], function(w) sum(w * spec[["prop"]][seq_along(w)]), numeric(1)),
                      foreground = child %in% fg_nodes, selected = cell == "B1" & child == 3L)
  utils::write.table(truth, file.path(out, "branch_truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  genes <- list()
  for (i in seq_len(n)) {
    id <- sprintf("%s_%s_%04d", cell, tree_name, i)
    aln <- sim_gene_nodes(tree, 300L, pi, par[["kappa"]], spec[["om"]], spec[["prop"]], spec[["per_branch"]])
    dropped <- character()
    if (tree_name == "T16m") dropped <- sample(tree$tip.label, 4L)
    aln <- aln[setdiff(tree$tip.label, dropped), , drop = FALSE]
    seqs <- apply(aln, 1, paste, collapse = "")
    writeLines(as.vector(rbind(paste0(">", names(seqs)), seqs)), file.path(out, "alignments", paste0(id, ".fasta")))
    genes[[i]] <- data.frame(gene = id, cell = cell, tree = tree_name, dropped = paste(dropped, collapse = ","))
  }
  utils::write.table(do.call(rbind, genes), file.path(out, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(out)
}

if (!interactive() && length(commandArgs(TRUE)) == 6L) {
  a <- commandArgs(TRUE)
  run_lineage_cell(a[1], a[2], a[3], as.integer(a[4]), as.integer(a[5]), a[6])
}
