# Multiple-testing mixture study summary (docs/dev/lineage-multiple-testing.md; measures and
# criteria fixed there).
#   Rscript --vanilla lineage-mt-summary.R <study-Q root> <out prefix>
# For each unit and run (absrel, branch_omega, and absrel_true where present): discoveries at 0.05,
# false-discovery proportion with a gene-bootstrap 95% interval (2,000 resamples, seed 1), power,
# and, for p_bh_lineage, per-family results and the families with no true alternative.
a <- commandArgs(TRUE); root <- a[1]; prefix <- a[2]
B <- 2000L
judge <- function(lo, hi) if (is.na(hi)) "no discoveries" else if (hi <= 0.075) "controls" else if (lo > 0.05) "fails" else "inconclusive"
fdp_boot <- function(gene, disc, false) {
  ids <- unique(gene); f <- factor(gene, levels = ids)
  d <- tapply(disc, f, sum); fa <- tapply(disc & false, f, sum)
  set.seed(1L)
  bs <- vapply(seq_len(B), function(i) { s <- sample(seq_along(ids), replace = TRUE); D <- sum(d[s]); if (D == 0) NA_real_ else sum(fa[s]) / D }, 0)
  c(lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)), hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)))
}
row <- function(unit, run, output, gene, p, truth) {
  ok <- !is.na(p); gene <- gene[ok]; p <- p[ok]; truth <- truth[ok]
  disc <- p <= 0.05; D <- sum(disc); F <- sum(disc & !truth)
  ci <- if (D > 0) fdp_boot(gene, disc, !truth) else c(lo = NA, hi = NA)
  data.frame(unit = unit, run = run, output = output, tests = length(p), true_tests = sum(truth), discoveries = D, false = F,
             fdp = if (D) F / D else NA_real_, fdp_lo = ci[["lo"]], fdp_hi = ci[["hi"]], judgement = judge(ci[["lo"]], ci[["hi"]]),
             power = if (sum(truth)) sum(disc & truth) / sum(truth) else NA_real_, stringsAsFactors = FALSE)
}
rd <- function(u, run, what) { f <- file.path(root, "units", u, run, sprintf("%s_%s.tsv", run, what)); if (file.exists(f)) utils::read.delim(f, stringsAsFactors = FALSE) else NULL }
units <- basename(dirname(Sys.glob(file.path(root, "units", "*", "DONE"))))
out <- list(); fam <- list()
for (u in units) {
  genes <- utils::read.delim(file.path(root, "units", u, "sim", "genes.tsv"), stringsAsFactors = FALSE)
  sel <- genes$gene[genes$selected]
  for (run in intersect(c("absrel", "absrel_true"), list.files(file.path(root, "units", u)))) {
    r <- rd(u, run, "results"); b <- rd(u, run, "branches")
    t <- r[r$status == "tested", ]; b <- b[b$orthogroup %in% t$orthogroup & !is.na(b$p_raw), ]
    # Truth per branch comes from the simulator (sim/branch_truth.tsv), not from the branches table's
    # `foreground` column, which is TRUE for every tested branch when aBSREL tests all branches.
    bt <- utils::read.delim(file.path(root, "units", u, "sim", "branch_truth.tsv"), stringsAsFactors = FALSE)
    truth_b <- b$orthogroup %in% sel & b$species_tree_branches %in% bt$branch_id[bt$foreground]
    out[[length(out) + 1L]] <- row(u, run, "p_raw per branch (uncorrected, for context)", b$orthogroup, b$p_raw, truth_b)
    out[[length(out) + 1L]] <- row(u, run, "p_holm per branch", b$orthogroup, b$p_holm, truth_b)
    out[[length(out) + 1L]] <- row(u, run, "p_bh_lineage, all families pooled", b$orthogroup, b$p_bh_lineage, truth_b)
    pooled <- stats::p.adjust(b$p_raw, "BH")
    out[[length(out) + 1L]] <- row(u, run, "M3 exploratory: one BH family over all branch tests", b$orthogroup, pooled, truth_b)
    out[[length(out) + 1L]] <- row(u, run, "gene p_value (unadjusted)", t$orthogroup, t$p_value, t$orthogroup %in% sel)
    out[[length(out) + 1L]] <- row(u, run, "gene p_adj", t$orthogroup, t$p_adj, t$orthogroup %in% sel)
    for (g in unique(b$species_tree_branches)) {
      i <- b$species_tree_branches == g
      x <- row(u, run, paste("p_bh_lineage family", g), b$orthogroup[i], b$p_bh_lineage[i], truth_b[i])
      x$family <- g; x$has_alternative <- any(truth_b[i]); fam[[length(fam) + 1L]] <- x
    }
  }
  bo <- rd(u, "branch_omega", "results")
  if (!is.null(bo)) {
    bo <- bo[bo$status == "tested", ]
    for (x in c("p_adj_two_vs_m0", "p_adj_fg_omega_gt1"))
      out[[length(out) + 1L]] <- row(u, "branch_omega", x, bo$orthogroup, bo[[x]], bo$orthogroup %in% sel)
  }
}
res <- do.call(rbind, out); fam <- do.call(rbind, fam)
w <- function(x, what) utils::write.table(x, paste0(prefix, "_", what, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w(res, "fdr"); w(fam, "families")
print(res[, c("unit", "run", "output", "discoveries", "false", "fdp", "fdp_lo", "fdp_hi", "judgement", "power")], digits = 3, row.names = FALSE)
nullfam <- fam[!fam$has_alternative, ]
for (k in unique(paste(nullfam$unit, nullfam$run))) {
  x <- nullfam[paste(nullfam$unit, nullfam$run) == k, ]
  ci <- stats::binom.test(sum(x$discoveries > 0), nrow(x))$conf.int
  cat(sprintf("%s: families with no alternative that rejected anything: %d of %d (95%% CI %.0f-%.0f%%)\n", k, sum(x$discoveries > 0), nrow(x), 100 * ci[1], 100 * ci[2]))
}
