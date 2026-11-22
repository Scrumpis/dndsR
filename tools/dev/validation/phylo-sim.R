# Codon simulator for calibrating dndsR's phylogenetic workflow. Independent of HyPhy and of
# dndsR's package code; uses ape to read trees and codon-model.R's codon table.
# Design: docs/dev/phylo-calibration-design.md.
#
#   Rscript --vanilla phylo-sim.R <params.json> <cell> <n genes> <seed> <out dir>
# Writes <out>/alignments/<gene>.fasta (or, for indel cells, CDS and tables for phylo_prepare),
# <out>/species_tree.nwk, <out>/foreground.tsv, and <out>/truth.tsv.
here <- local({
  sourced <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL) # when source()d
  if (!is.null(sourced)) return(dirname(sourced))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(f) else "tools/dev/validation"
})
source(file.path(here, "codon-model.R")) # bases, codons, code, sense
n_sense <- length(sense)
pos <- do.call(rbind, strsplit(sense, ""))
syn_pair <- outer(code[sense], code[sense], "==")
ndiff <- outer(seq_len(n_sense), seq_len(n_sense), Vectorize(function(i, j) sum(pos[i, ] != pos[j, ])))
is_ts <- function(a, b) paste(sort(c(a, b)), collapse = "") %in% c("AG", "CT")
# Product over differing positions of the HKY factor (kappa for a transition, 1 otherwise).
kfac_cache <- new.env()
kappa_factor <- function(kappa) {
  key <- format(kappa)
  if (!is.null(kfac_cache[[key]])) return(kfac_cache[[key]])
  f <- matrix(0, n_sense, n_sense)
  for (i in seq_len(n_sense)) for (j in seq_len(n_sense)) {
    d <- which(pos[i, ] != pos[j, ])
    if (length(d) %in% 1:2) f[i, j] <- prod(vapply(d, function(k) if (is_ts(pos[i, k], pos[j, k])) kappa else 1, numeric(1)))
  }
  kfac_cache[[key]] <- f
  f
}

#' Codon frequencies from per-position nucleotide frequencies (F3x4), over sense codons
f3x4 <- function(nuc) { # 3 x 4 matrix, columns T C A G
  colnames(nuc) <- bases
  p <- vapply(sense, function(cd) prod(nuc[cbind(1:3, match(strsplit(cd, "")[[1]], bases))]), numeric(1))
  p / sum(p)
}

#' Rate matrix: HKY mutation, frequencies pi, omega on nonsynonymous changes, synonymous rate r,
#' and double changes at rate delta. Scaled so that omega = 1, r = 1, delta = 0 gives one
#' expected substitution per codon per unit of time.
rate_q <- function(kappa, omega, r, delta, pi) {
  f <- kappa_factor(kappa)
  q <- matrix(0, n_sense, n_sense)
  single <- ndiff == 1L; double <- ndiff == 2L
  q[single] <- f[single]
  q[double] <- delta * f[double]
  q <- q * ifelse(syn_pair, r, omega)
  q <- sweep(q, 2, pi, "*")
  base <- f; base[!single] <- 0
  base <- sweep(base, 2, pi, "*")
  scale <- sum(pi * rowSums(base))
  q <- q / scale
  diag(q) <- -rowSums(q)
  q
}

p_cache <- new.env()
p_matrix <- function(kappa, omega, r, delta, pi, t) {
  key <- paste(kappa, omega, r, delta, t, sep = "|")
  if (!is.null(p_cache[[key]])) return(p_cache[[key]])
  ekey <- paste("e", kappa, omega, r, delta, sep = "|")
  if (is.null(p_cache[[ekey]])) {
    q <- rate_q(kappa, omega, r, delta, pi)
    # Symmetrize for a stable decomposition of a reversible matrix.
    s <- sqrt(pi); sym <- diag(s) %*% q %*% diag(1 / s); sym <- (sym + t(sym)) / 2
    e <- eigen(sym, symmetric = TRUE)
    p_cache[[ekey]] <- list(v = e$vectors, d = e$values, s = s)
  }
  e <- p_cache[[ekey]]
  p <- diag(1 / e$s) %*% e$v %*% diag(exp(e$d * t)) %*% t(e$v) %*% diag(e$s)
  p[p < 0] <- 0
  p <- p / rowSums(p)
  p_cache[[key]] <- p
  p
}

#' Simulate one gene on a tree.
#' @param tree ape phylo with edge lengths (expected substitutions per codon).
#' @param fg_nodes Nodes whose parent edge is foreground.
#' @param omega_bg,omega_fg,prop Per-class omega on background and foreground branches, and weights.
#' @param srv Multipliers and weights of synonymous rate classes.
#' @return A character matrix of codons, taxa by sites (true alignment, "---" for gaps).
sim_gene <- function(tree, n_sites, pi, kappa, omega_bg, omega_fg, prop, srv, delta = 0,
                     fg_nodes = integer(), indel_rate = 0, indel_mean = 2, cluster = FALSE) {
  n_tip <- ape::Ntip(tree); root <- n_tip + 1L
  block <- NULL
  if (cluster) { # stage 2B: the last class as one contiguous block of round(prop * sites) codons
    k <- length(prop); m <- round(prop[k] * n_sites)
    cls <- sample.int(k - 1L, n_sites, TRUE, prob = prop[-k])
    at <- sample.int(n_sites - m + 1L, 1L); cls[at:(at + m - 1L)] <- k; block <- c(at, at + m - 1L)
  } else cls <- sample.int(length(prop), n_sites, TRUE, prob = prop)
  rcls <- sample.int(length(srv$rate), n_sites, TRUE, prob = srv$weight)
  states <- vector("list", max(tree$edge))
  states[[root]] <- sample.int(n_sense, n_sites, TRUE, prob = pi)
  tr <- ape::reorder.phylo(tree, "cladewise") # parents before children
  ord <- tr$edge
  lens <- tr$edge.length
  for (e in seq_len(nrow(ord))) {
    par <- ord[e, 1]; ch <- ord[e, 2]
    om <- if (ch %in% fg_nodes) omega_fg else omega_bg
    out <- integer(n_sites)
    for (k in unique(cls)) for (m in unique(rcls)) {
      idx <- which(cls == k & rcls == m)
      if (!length(idx)) next
      p <- p_matrix(kappa, om[k], srv$rate[m], delta, pi, lens[e])
      from <- states[[par]][idx]
      for (a in unique(from)) {
        w <- idx[from == a]
        out[w] <- sample.int(n_sense, length(w), TRUE, prob = p[a, ])
      }
    }
    states[[ch]] <- out
  }
  aln <- t(vapply(seq_len(n_tip), function(i) sense[states[[i]]], character(n_sites)))
  rownames(aln) <- tree$tip.label
  if (indel_rate > 0) aln <- add_indels(aln, tree, ord, lens, indel_rate, indel_mean, pi)
  attr(aln, "block") <- block
  aln
}

#' Codon indels with the true alignment kept. Insertions add columns present only below the
#' branch where they happen; deletions gap existing columns below their branch.
# homologous = FALSE (the stage 2 M4 generator, kept for its frozen studies) gives each descendant
# of an inserting branch its own random codons, so an insertion on an internal branch is not
# homologous across those descendants. homologous = TRUE draws one inserted sequence and copies it
# to every descendant (identical by descent; it does not evolve further).
add_indels <- function(aln, tree, ord, lens, rate, mean_len, pi, homologous = FALSE) {
  n_tip <- ape::Ntip(tree)
  desc <- function(node) {
    if (node <= n_tip) return(node)
    unlist(lapply(tree$edge[tree$edge[, 1] == node, 2], desc))
  }
  for (e in seq_len(nrow(ord))) {
    below <- desc(ord[e, 2])
    n_ev <- stats::rpois(1, rate * ncol(aln) * lens[e])
    for (k in seq_len(n_ev)) {
      len <- stats::rgeom(1, 1 / mean_len) + 1L
      at <- sample.int(ncol(aln), 1L)
      if (stats::runif(1) < 0.5) {
        cols <- at:min(ncol(aln), at + len - 1L)
        aln[below, cols] <- "---"
      } else {
        new <- matrix("---", nrow(aln), len)
        new[below, ] <- if (homologous) matrix(sample(sense, len, TRUE, prob = pi), length(below), len, byrow = TRUE) else
          sample(sense, length(below) * len, TRUE, prob = pi)
        aln <- cbind(aln[, seq_len(at), drop = FALSE], new, aln[, -seq_len(at), drop = FALSE])
      }
    }
  }
  keep <- colSums(aln != "---") > 0
  aln[, keep, drop = FALSE]
}

#' A misaligned window: 12 codons of one taxon replaced by random sense codons.
misalign <- function(aln, taxon, width = 12L, pi) {
  at <- sample.int(ncol(aln) - width + 1L, 1L)
  aln[taxon, at:(at + width - 1L)] <- sample(sense, width, TRUE, prob = pi)
  aln
}

srv_default <- function() { r <- c(0.3, 1, 2.2); w <- c(0.3, 0.5, 0.2); list(rate = r / sum(r * w), weight = w) }
srv_strong <- function() { # gamma shape 0.4, 8 equal-weight classes by their conditional means
  a <- 0.4; q <- stats::qgamma(seq(0, 1, length.out = 9), a, a)
  m <- vapply(1:8, function(i) 8 * (stats::pgamma(q[i + 1], a + 1, a) - stats::pgamma(q[i], a + 1, a)), numeric(1))
  list(rate = m / mean(m), weight = rep(1 / 8, 8))
}

#' Cell definitions. All take the shared parameters `par` (from params.json).
cell_spec <- function(cell, par) {
  n0 <- list(omega_bg = c(0.05, 0.3, 1.0), omega_fg = c(0.05, 0.3, 1.0), prop = c(0.6, 0.3, 0.1))
  base <- c(n0, list(length = par$L50, scale = 1, srv = srv_default(), delta = 0, indel = 0,
                     misalign = FALSE, gene_tree = NULL, alternative = FALSE, cluster = FALSE))
  alt <- function(w, p) { s <- base; s$omega_fg <- c(0.05, 0.3, w); s$prop <- c(0.6, 0.3, 0.1) * (1 - p) / 0.9; s$prop[3] <- p; s$alternative <- TRUE; s }
  switch(cell,
    N0_L25 = modifyList(base, list(length = par$L25)),
    N0_L50 = base,
    N0_L75 = modifyList(base, list(length = par$L75)),
    N0_L50_div05 = modifyList(base, list(scale = 0.5)),
    N0_L50_div2 = modifyList(base, list(scale = 2)),
    N1_L50 = modifyList(base, list(omega_bg = c(0.02, 0.15, 0.5), omega_fg = c(0.02, 0.15, 0.5))),
    A3_p05 = alt(3, 0.05), A3_p20 = alt(3, 0.2), A10_p05 = alt(10, 0.05), A10_p20 = alt(10, 0.2), A10_p10 = alt(10, 0.1),
    # Added for stage 2 (phylo-calibration-stage2.md); existing cells are unchanged.
    A5_p10 = alt(5, 0.1), M1_mnp_d01 = modifyList(base, list(delta = 0.1)), M1_mnp_d03 = modifyList(base, list(delta = 0.3)),
    M1_mnp = modifyList(base, list(delta = 0.2)),
    # Added for stage 2B (phylo-calibration-stage2b.md); existing cells are unchanged.
    A2_p20 = alt(2, 0.2), A3_p10 = alt(3, 0.1), A5_p20 = alt(5, 0.2),
    N0_L25_div05 = modifyList(base, list(length = par$L25, scale = 0.5)), N0_L25_div2 = modifyList(base, list(length = par$L25, scale = 2)),
    N0_L75_div05 = modifyList(base, list(length = par$L75, scale = 0.5)), N0_L75_div2 = modifyList(base, list(length = par$L75, scale = 2)),
    AC10_p05 = modifyList(alt(10, 0.05), list(cluster = TRUE)), AC5_p10 = modifyList(alt(5, 0.1), list(cluster = TRUE)),
    A10_p10_indel = modifyList(alt(10, 0.1), list(indel = 0.02)),
    M2_srv = modifyList(base, list(srv = srv_strong())),
    M3_misalign = modifyList(base, list(misalign = TRUE)),
    M4_indel = modifyList(base, list(indel = 0.02)),
    M5_discordant = modifyList(base, list(gene_tree = par$discordant_tree)),
    stop("Unknown cell: ", cell)
  )
}

run_cell <- function(params_file, cell, n, seed, out) {
  par <- jsonlite::read_json(params_file, simplifyVector = TRUE)
  spec <- cell_spec(cell, par)
  RNGversion("4.0.0"); set.seed(seed)
  pi <- f3x4(matrix(unlist(par$nuc_f3x4), 3, byrow = TRUE))
  species <- ape::read.tree(text = par$tree)
  species$edge.length <- species$edge.length * par$c * spec$scale
  sim_tree <- if (is.null(spec$gene_tree)) species else {
    g <- ape::read.tree(text = spec$gene_tree); g$edge.length <- g$edge.length * par$c * spec$scale; g
  }
  fg_tips <- par$foreground
  fg_nodes <- c(match(fg_tips, sim_tree$tip.label), ape::getMRCA(sim_tree, fg_tips))
  if (!ape::is.monophyletic(sim_tree, fg_tips)) fg_nodes <- match(fg_tips, sim_tree$tip.label)
  dir.create(file.path(out, "alignments"), recursive = TRUE, showWarnings = FALSE)
  writeLines(par$species_topology, file.path(out, "species_tree.nwk"))
  writeLines(c("scope\ttaxa", paste0("clade\t", paste(fg_tips, collapse = ","))), file.path(out, "foreground.tsv"))
  truth <- list()
  cds <- list()
  for (i in seq_len(n)) {
    id <- sprintf("%s_%04d", cell, i)
    aln <- sim_gene(sim_tree, spec$length, pi, par$kappa, spec$omega_bg, spec$omega_fg, spec$prop, spec$srv,
                    spec$delta, fg_nodes, spec$indel, cluster = isTRUE(spec$cluster))
    block <- attr(aln, "block")
    if (spec$misalign) aln <- misalign(aln, fg_tips[1], 12L, pi)
    aln <- aln[species$tip.label, , drop = FALSE]
    seqs <- apply(aln, 1, paste, collapse = "")
    if (spec$indel > 0) {
      for (tx in names(seqs)) cds[[length(cds) + 1L]] <- data.frame(orthogroup = id, taxon = tx, gene_id = paste0(tx, "_", id),
                                                                      seq = gsub("-", "", seqs[[tx]]))
      writeLines(as.vector(rbind(paste0(">", names(seqs)), seqs)), file.path(out, "alignments", paste0(id, ".true.fasta")))
    } else {
      writeLines(as.vector(rbind(paste0(">", names(seqs)), seqs)), file.path(out, "alignments", paste0(id, ".fasta")))
    }
    truth[[i]] <- data.frame(gene = id, cell = cell, alternative = spec$alternative, codons = ncol(aln))
    if (isTRUE(spec$cluster)) { truth[[i]]$selected_start <- block[1]; truth[[i]]$selected_end <- block[2] }
  }
  write.table(do.call(rbind, truth), file.path(out, "truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (length(cds)) {
    cds <- do.call(rbind, cds)
    dir.create(file.path(out, "cds"), showWarnings = FALSE)
    tab <- c("taxon\tcds")
    for (tx in unique(cds$taxon)) {
      f <- file.path(out, "cds", paste0(tx, ".fasta"))
      writeLines(as.vector(rbind(paste0(">", cds$gene_id[cds$taxon == tx]), cds$seq[cds$taxon == tx])), f)
      tab <- c(tab, paste0(tx, "\tcds/", tx, ".fasta"))
    }
    writeLines(tab, file.path(out, "cds_table.tsv"))
    write.table(cds[, c("orthogroup", "taxon", "gene_id")], file.path(out, "orthogroup_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(out)
}

if (!interactive() && length(commandArgs(TRUE)) == 5L) {
  a <- commandArgs(TRUE)
  run_cell(a[1], a[2], as.integer(a[3]), as.integer(a[4]), a[5])
}
