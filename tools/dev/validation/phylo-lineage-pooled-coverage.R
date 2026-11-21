# Exploratory (docs/dev/validation-audit.md, E1). Usage: Rscript --vanilla phylo-lineage-pooled-coverage.R <first study root> <out.tsv>
# lineage calibration study's outputs (known per-branch truth). Same estimator and interval as
# R/phylo_summarize.R (high_ds excluded, 1000 resamples of orthogroups). Disjoint gene subsets.
root <- commandArgs(TRUE)[1]; out <- commandArgs(TRUE)[2]
pooled <- function(NdN, N, SdS, S) { ok <- is.finite(NdN + N + SdS + S); if (!any(ok) || sum(SdS[ok]) <= 0) return(NA); (sum(NdN[ok]) / sum(N[ok])) / (sum(SdS[ok]) / sum(S[ok])) }
set.seed(20261121); rows <- list()
for (u in c("A0_T5","A0_T8b","A0_T8p","A0_T16","A0_T16m","A1_T5","A1_T8b","A1_T16","A2_T8b","A2_T8p","A2_T16","A3_T5","A3_T8b")) {
  d <- file.path(root, "units", u); tr <- read.delim(file.path(d, "sim", "branch_truth.tsv"))
  r <- read.delim(file.path(d, "branch_omega", "branch_omega_results.tsv")); b <- read.delim(file.path(d, "branch_omega", "branch_omega_branches.tsv"))
  b <- b[b$orthogroup %in% r$orthogroup[r$status == "tested"] & b$in_species_tree & !grepl("+", b$species_tree_branches, fixed = TRUE) & !grepl("high_ds", b$flags), ]
  genes <- sample(unique(b$orthogroup))
  for (n in c(25, 50, 100)) {
    groups <- split(genes, ceiling(seq_along(genes) / n)); groups <- groups[lengths(groups) == n]
    for (g in seq_along(groups)) for (id in tr$branch_id) {
      x <- b[b$orthogroup %in% groups[[g]] & b$species_tree_branches == id, ]
      if (nrow(x) < 5) next
      est <- pooled(x$NdN_free, x$N_free, x$SdS_free, x$S_free)
      bo <- replicate(1000, { i <- sample.int(nrow(x), replace = TRUE); pooled(x$NdN_free[i], x$N_free[i], x$SdS_free[i], x$S_free[i]) })
      ci <- quantile(bo, c(.025, .975), na.rm = TRUE); t <- tr$omega_mean[tr$branch_id == id]
      rows[[length(rows) + 1]] <- data.frame(unit = u, n = n, group = g, branch_id = id, terminal = !grepl("^n", id), length = tr$length_codon[tr$branch_id == id],
        truth = t, est = est, lo = ci[1], hi = ci[2], covered = ci[1] <= t & t <= ci[2], below = ci[2] < t, above = ci[1] > t, n_genes = nrow(x))
    }
  }
}
x <- do.call(rbind, rows); write.table(x, out, sep = "\t", quote = FALSE, row.names = FALSE)
cp <- function(k, n) sprintf("%.1f%% (%.1f-%.1f)", 100 * k / n, 100 * binom.test(k, n)$conf.int[1], 100 * binom.test(k, n)$conf.int[2])
cat("By cell family and n (intervals: Clopper-Pearson, ignoring dependence among branches of one gene set):\n")
x$fam <- sub("_.*", "", x$unit)
for (f in unique(x$fam)) for (n in c(25, 50, 100)) { y <- x[x$fam == f & x$n == n, ]; if (nrow(y)) cat(f, "n =", n, ":", nrow(y), "intervals, covered", cp(sum(y$covered), nrow(y)), "; truth above interval", sum(y$below), ", below", sum(y$above), "; median log2 bias", round(median(log2(y$est / y$truth)), 3), "\n") }
cat("By branch length tercile (all, n = 100):\n"); y <- x[x$n == 100, ]; y$lt <- cut(y$length, quantile(y$length, 0:3 / 3), include.lowest = TRUE)
for (l in levels(y$lt)) { z <- y[y$lt == l, ]; cat(l, ":", cp(sum(z$covered), nrow(z)), " median log2 bias", round(median(log2(z$est / z$truth)), 3), "\n") }
