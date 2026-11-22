# Simulator for studies L, B, and R-D (docs/dev/lineage-summary-study.md,
# lineage-robustness-study.md, relax-bootstrap-study.md). Reuses phylo-lineage-sim.R (trees, codon
# model, the first studies' cells) without changing it, so the earlier studies stay reproducible.
# Independent of dndsR's package code.
#
#   Rscript --vanilla phylo-lineage-sim2.R <params.json> <cell> <tree> <n genes> <seed> <out dir> <codons>
#
# Always 7 arguments (phylo-lineage-sim.R's own entry point runs only with 6).
# Writes, as phylo-lineage-sim.R does: alignments/<gene>.fasta, species_tree.nwk, foreground.tsv,
# genes.tsv, branch_truth.tsv; and gene_branch_truth.tsv (each gene's true mean omega on each
# branch, for estimands over a particular gene set). Indel cells write cds/, cds_table.tsv,
# orthogroup_table.tsv, and alignments/<gene>.true.fasta instead of alignments/<gene>.fasta.
#
# Cells
#   Study L (omega classes scale with gene and branch effects; truth = opportunity-weighted mean):
#     LS0  one omega 0.2                       LS1  classes 0.05/0.3/1 at 0.6/0.3/0.1
#     LS2  gamma shape 0.5, mean 0.2, 8 classes (conditional means)
#     LG1  LS2 x gene effect (lognormal, sd 0.7, mean 1) x branch effect (lognormal sd 0.5, one fixed draw per tree)
#     LA1  LG1, and genes with a gene effect in the top third of its distribution lose t1 with probability 0.5
#   Study B (null; site classes 0.05/0.3/1 at 0.6/0.3/0.1 on every branch, as B0):
#     BE0  no error                            BE3  codon indels (rate 0.02 per codon per unit length, mean 2 codons)
#     BE3h as BE3, with insertions homologous across the inserting branch's descendants (add_indels(homologous = TRUE))
#     BE3hS as BE3h, with every gene selected on the foreground (classes 0.05/0.3/10 there; alignment pilot)
#     BMIXh BMIX's selection (10% of genes on the foreground) with BE3h's indels, no corruption (multiple-testing study)
#     BE4  t1 replaced by a paralog that split from the species lineage at the root
#     BE5  gene tree discordant: t2 and t3 exchanged (one NNI at the stem of the foreground clade)
#     BE7  recombination: the second half of the codons has BE5's gene tree
#     BMIX genome-like mixture: 10% of genes selected on the foreground (classes 0.05/0.3/10 there)
#     (BE1, BE2, and BMIX's errors are added afterwards by phylo-windows.py; BE6 is BE0 on T16m)
#   Study R-D: any cell of phylo-lineage-sim.R (for example C0), at the given number of codons.
here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(f) else "tools/dev/validation"
})
source(file.path(here, "phylo-lineage-sim.R"))

gamma_classes <- function(shape, mean, k) {
  rate <- shape / mean
  q <- stats::qgamma(seq_len(k - 1L) / k, shape, rate)
  mean * k * diff(c(0, stats::pgamma(q, shape + 1, rate), 1))
}

#' Fixed per-tree branch effects (one draw, the same for every gene and seed)
branch_effects <- function(n_nodes) {
  old <- if (exists(".Random.seed", globalenv())) .Random.seed else NULL
  set.seed(20261122)
  e <- exp(0.5 * stats::rnorm(n_nodes))
  if (!is.null(old)) assign(".Random.seed", old, globalenv())
  e
}

b0 <- list(om = c(0.05, 0.3, 1), prop = c(0.6, 0.3, 0.1), gene_sd = 0, branch_sd = FALSE)

#' Base classes, gene-effect model, and error of a study L or B cell
cell2 <- function(cell) {
  g8 <- gamma_classes(0.5, 0.2, 8L)
  switch(cell,
    LS0 = list(om = 0.2, prop = 1, gene_sd = 0, branch_sd = FALSE),
    LS1 = list(om = b0$om, prop = b0$prop, gene_sd = 0, branch_sd = FALSE),
    LS2 = list(om = g8, prop = rep(1 / 8, 8), gene_sd = 0, branch_sd = FALSE),
    LG1 = list(om = g8, prop = rep(1 / 8, 8), gene_sd = 0.7, branch_sd = TRUE),
    LA1 = list(om = g8, prop = rep(1 / 8, 8), gene_sd = 0.7, branch_sd = TRUE, loss = TRUE),
    BE0 = c(b0, error = "none"), BE3 = c(b0, error = "indel"), BE3h = c(b0, error = "indel_h"), BE3hS = c(b0, error = "indel_h", select_all = TRUE), BMIXh = c(b0, error = "indel_h", mix = TRUE), BE4 = c(b0, error = "paralog"),
    BE5 = c(b0, error = "discordant"), BE7 = c(b0, error = "recombination"), BMIX = c(b0, error = "mixture"),
    NULL)
}

write_fasta <- function(aln, path) {
  seqs <- apply(aln, 1, paste, collapse = "")
  writeLines(as.vector(rbind(paste0(">", names(seqs)), seqs)), path)
}

swap_rows <- function(aln, a, b, cols = seq_len(ncol(aln))) {
  tmp <- aln[a, cols]; aln[a, cols] <- aln[b, cols]; aln[b, cols] <- tmp; aln
}

#' A gene on `tree` plus a paralog of tip 1 that split before the root: a new root joins the
#' species tree (stem half of t1's depth) and the paralog (t1's depth plus that stem). The stem was
#' 10% of the depth in the pilot, where t1's branch reached 3 times the median in only 75% of genes
#' (gate: 90%); changed before the main runs (docs/dev/lineage-robustness-study.md).
paralog_tree <- function(tree) {
  depth <- ape::node.depth.edgelength(tree)[1]
  stem <- 0.5 * depth
  sp <- ape::write.tree(tree); sp <- sub(";$", "", sp)
  ape::read.tree(text = sprintf("(%s:%.10g,t1p:%.10g);", sp, stem, depth + stem))
}

run_cell2 <- function(params_file, cell, tree_name, n, seed, out, codons) {
  par <- jsonlite::read_json(params_file, simplifyVector = TRUE)
  RNGversion("4.0.0"); set.seed(seed)
  pi <- f3x4(matrix(unlist(par$nuc_f3x4), 3, byrow = TRUE))
  tree <- ape::read.tree(text = lineage_trees(par)[[tree_name]])
  tree$edge.length <- tree$edge.length * par$c
  tree$node.label <- paste0("n", seq_len(tree$Nnode))
  n_nodes <- ape::Ntip(tree) + tree$Nnode
  fg_tips <- tree$tip.label[1:2]
  fg_nodes <- c(match(fg_tips, tree$tip.label), ape::getMRCA(tree, fg_tips))
  spec <- cell2(cell)
  legacy <- is.null(spec)
  if (legacy) spec0 <- cell_omegas(cell, tree, fg_nodes)
  dir.create(file.path(out, "alignments"), recursive = TRUE, showWarnings = FALSE)
  topo <- tree; topo$edge.length <- NULL
  ape::write.tree(topo, file.path(out, "species_tree.nwk"))
  writeLines(c("scope\ttaxa", paste0("clade\t", paste(fg_tips, collapse = ","))), file.path(out, "foreground.tsv"))
  names_all <- c(tree$tip.label, tree$node.label)
  child <- tree$edge[, 2]
  eff <- if (!legacy && isTRUE(spec$branch_sd)) branch_effects(n_nodes) else rep(1, n_nodes)
  base_mean <- if (legacy) vapply(spec0$om, function(w) sum(w * spec0$prop[seq_along(w)]), numeric(1)) else rep(sum(spec$om * spec$prop), n_nodes)
  # Population truth: the gene effect has mean 1, so the opportunity-weighted mean omega of a
  # branch over the gene population is its base mean times its branch effect.
  truth <- data.frame(branch_id = names_all[child], length_codon = tree$edge.length,
                      omega_mean = base_mean[child] * eff[child], branch_effect = eff[child],
                      foreground = child %in% fg_nodes)
  utils::write.table(truth, file.path(out, "branch_truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  loss_cut <- if (!legacy && isTRUE(spec$loss)) stats::qlnorm(2 / 3, -spec$gene_sd^2 / 2, spec$gene_sd) else Inf
  ptree <- if (!legacy && identical(spec$error, "paralog")) paralog_tree(tree) else NULL
  genes <- list(); gbt <- list(); cds <- list()
  for (i in seq_len(n)) {
    id <- sprintf("%s_%s_%04d", cell, tree_name, i)
    m <- 1; selected <- FALSE
    if (legacy) {
      om <- spec0$om; prop <- spec0$prop; per_branch <- spec0$per_branch
    } else {
      if (spec$gene_sd > 0) m <- exp(stats::rnorm(1, -spec$gene_sd^2 / 2, spec$gene_sd))
      om <- lapply(seq_len(n_nodes), function(k) spec$om * m * eff[k]); prop <- spec$prop; per_branch <- FALSE
      if ((identical(spec$error, "mixture") || isTRUE(spec$mix)) && stats::runif(1) < 0.1) {
        selected <- TRUE; om[fg_nodes] <- list(c(0.05, 0.3, 10))
      }
      if (isTRUE(spec$select_all)) { selected <- TRUE; om[fg_nodes] <- list(c(0.05, 0.3, 10)) }
    }
    if (!is.null(ptree)) {
      aln <- sim_gene_nodes(ptree, codons, pi, par$kappa, rep(list(spec$om), ape::Ntip(ptree) + ptree$Nnode), prop)
      aln["t1", ] <- aln["t1p", ]; aln <- aln[tree$tip.label, , drop = FALSE]
    } else {
      aln <- sim_gene_nodes(tree, codons, pi, par$kappa, om, prop, per_branch)
    }
    if (!legacy && identical(spec$error, "discordant")) aln <- swap_rows(aln, "t2", "t3")
    if (!legacy && identical(spec$error, "recombination")) aln <- swap_rows(aln, "t2", "t3", (codons %/% 2 + 1):codons)
    dropped <- character()
    if (tree_name == "T16m") dropped <- sample(tree$tip.label, 4L)
    if (m > loss_cut && stats::runif(1) < 0.5) dropped <- union(dropped, "t1")
    if (!legacy && spec$error %in% c("indel", "indel_h")) {
      tr <- ape::reorder.phylo(tree, "cladewise")
      aln <- add_indels(aln, tree, tr$edge, tr$edge.length, 0.02, 2, pi, homologous = identical(spec$error, "indel_h"))
    }
    aln <- aln[setdiff(tree$tip.label, dropped), , drop = FALSE]
    if (!legacy && spec$error %in% c("indel", "indel_h")) {
      write_fasta(aln, file.path(out, "alignments", paste0(id, ".true.fasta")))
      seqs <- apply(aln, 1, paste, collapse = "")
      for (tx in names(seqs)) cds[[length(cds) + 1L]] <- data.frame(orthogroup = id, taxon = tx, gene_id = paste0(tx, "_", id), seq = gsub("-", "", seqs[[tx]]))
    } else {
      write_fasta(aln, file.path(out, "alignments", paste0(id, ".fasta")))
    }
    gm <- vapply(om[child], function(w) sum(w * prop[seq_along(w)]), numeric(1))
    gbt[[i]] <- data.frame(gene = id, branch_id = names_all[child], omega_mean = gm)
    genes[[i]] <- data.frame(gene = id, cell = cell, tree = tree_name, codons = codons, gene_effect = m,
                             selected = selected, dropped = paste(dropped, collapse = ","))
  }
  utils::write.table(do.call(rbind, genes), file.path(out, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(do.call(rbind, gbt), file.path(out, "gene_branch_truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (length(cds)) {
    cds <- do.call(rbind, cds)
    dir.create(file.path(out, "cds"), showWarnings = FALSE)
    tab <- "taxon\tcds"
    for (tx in unique(cds$taxon)) {
      writeLines(as.vector(rbind(paste0(">", cds$gene_id[cds$taxon == tx]), cds$seq[cds$taxon == tx])), file.path(out, "cds", paste0(tx, ".fasta")))
      tab <- c(tab, paste0(tx, "\tcds/", tx, ".fasta"))
    }
    writeLines(tab, file.path(out, "cds_table.tsv"))
    utils::write.table(cds[, c("orthogroup", "taxon", "gene_id")], file.path(out, "orthogroup_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(out)
}

if (!interactive() && length(commandArgs(TRUE)) == 7L) {
  a <- commandArgs(TRUE)
  run_cell2(a[1], a[2], a[3], as.integer(a[4]), as.integer(a[5]), a[6], as.integer(a[7]))
}
