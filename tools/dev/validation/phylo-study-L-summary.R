# Study L summary (docs/dev/lineage-summary-study.md): estimators, intervals, bias, ascertainment,
# and the trait-model consequence, from phylo-lineage-study2.sb units.
#   Rscript --vanilla phylo-study-L-summary.R <results root> <out prefix> [units, comma-separated] [sizes]
# Same branch filter as R/phylo_summarize.R: tested orthogroups, species-tree branches that are not
# merged ("+"), high_ds excluded. Disjoint random gene sets of each size (seed 20261123).
# Writes <prefix>_intervals.tsv (one row per gene set, branch, and estimator), <prefix>_coverage.tsv,
# <prefix>_bias.tsv, <prefix>_ascertainment.tsv, <prefix>_trait.tsv, and prints the criteria.
a <- commandArgs(TRUE)
root <- a[1]; prefix <- a[2]
units <- if (length(a) > 2 && nzchar(a[3])) strsplit(a[3], ",")[[1]] else sub("/DONE$", "", sub(paste0(root, "/units/"), "", Sys.glob(file.path(root, "units", "*", "DONE"))))
sizes <- if (length(a) > 3) as.integer(strsplit(a[4], ",")[[1]]) else c(10L, 25L, 50L, 100L, 200L)
B <- 1000L
set.seed(20261123)
ratio <- function(NdN, N, SdS, S) { s <- sum(SdS) / sum(S); if (!is.finite(s) || s <= 0) NA_real_ else (sum(NdN) / sum(N)) / s }

load_unit <- function(u) {
  d <- file.path(root, "units", u)
  r <- read.delim(file.path(d, "branch_omega", "branch_omega_results.tsv"), stringsAsFactors = FALSE)
  b <- read.delim(file.path(d, "branch_omega", "branch_omega_branches.tsv"), stringsAsFactors = FALSE)
  tested <- r$orthogroup[r$status == "tested"]
  b_all <- b[b$orthogroup %in% tested, ]
  b <- b_all[b_all$in_species_tree & !grepl("+", b_all$species_tree_branches, fixed = TRUE) & !grepl("high_ds", b_all$flags), ]
  list(unit = u, b = b, b_all = b_all, tested = tested, n_genes = nrow(r),
       truth = read.delim(file.path(d, "sim", "branch_truth.tsv"), stringsAsFactors = FALSE),
       gtruth = read.delim(file.path(d, "sim", "gene_branch_truth.tsv"), stringsAsFactors = FALSE),
       genes = read.delim(file.path(d, "sim", "genes.tsv"), stringsAsFactors = FALSE))
}

#' Estimators on one branch's rows x (one per gene); g0 holds each of those genes' M0 totals over all branches
estimates <- function(x, g0) {
  p <- ratio(x$NdN_free, x$N_free, x$SdS_free, x$S_free)
  c(pool = p, median = stats::median(x$omega_free), rel = p / ratio(g0$NdN, g0$N, g0$SdS, g0$S), ds = sum(x$SdS_free) / sum(x$S_free))
}

intervals <- list(); trait_sets <- list()
for (u in units) {
  U <- load_unit(u)
  tr <- U$truth
  # Truth of `rel`: the branch's truth over the tree-wide truth, length-weighted, as M0's pooled counts weight branches
  tw <- sum(tr$length_codon * tr$omega_mean) / sum(tr$length_codon)
  g0 <- stats::aggregate(cbind(NdN = NdN_m0, N = N_m0, SdS = SdS_m0, S = S_m0) ~ orthogroup, U$b_all, sum)
  genes <- sample(U$tested)
  for (n in sizes) {
    groups <- split(genes, ceiling(seq_along(genes) / n)); groups <- groups[lengths(groups) == n]
    for (g in seq_along(groups)) {
      gs <- groups[[g]]
      bb <- U$b[U$b$orthogroup %in% gs, ]
      for (id in tr$branch_id) {
        x <- bb[bb$species_tree_branches == id, ]
        if (nrow(x) < 5) next
        gx <- g0[match(x$orthogroup, g0$orthogroup), ]
        est <- estimates(x, gx)
        boot <- replicate(B, { i <- sample.int(nrow(x), replace = TRUE); estimates(x[i, ], gx[i, ])[c("pool", "median", "rel")] })
        t_pop <- tr$omega_mean[tr$branch_id == id]
        t_smp <- mean(U$gtruth$omega_mean[U$gtruth$branch_id == id & U$gtruth$gene %in% x$orthogroup])
        for (e in c("pool", "median", "rel")) {
          bo <- boot[e, ]; bo <- bo[is.finite(bo)]
          pct <- stats::quantile(bo, c(0.025, 0.975), names = FALSE)
          se <- stats::sd(log(bo[bo > 0])); tq <- stats::qt(0.975, nrow(x) - 1)
          logt <- if (is.finite(est[[e]]) && est[[e]] > 0) exp(log(est[[e]]) + c(-1, 1) * tq * se) else c(NA, NA)
          tp <- if (e == "rel") t_pop / tw else t_pop
          ts <- if (e == "rel") NA else t_smp
          intervals[[length(intervals) + 1L]] <- data.frame(unit = u, n = n, set = g, branch_id = id, length = tr$length_codon[tr$branch_id == id],
            estimator = e, n_genes = nrow(x), est = est[[e]], truth_pop = tp, truth_sample = ts,
            pct_lo = pct[1], pct_hi = pct[2], logt_lo = logt[1], logt_hi = logt[2], ds = est[["ds"]])
        }
      }
    }
  }
}
iv <- do.call(rbind, intervals)
write.table(iv, paste0(prefix, "_intervals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Coverage, with a cluster bootstrap over gene sets (branches of one set share genes)
cov_rows <- list()
for (key in unique(paste(iv$unit, iv$n, iv$estimator))) {
  x <- iv[paste(iv$unit, iv$n, iv$estimator) == key, ]
  for (truth in c("truth_pop", "truth_sample")) for (ci in c("pct", "logt")) {
    t <- x[[truth]]; lo <- x[[paste0(ci, "_lo")]]; hi <- x[[paste0(ci, "_hi")]]
    ok <- is.finite(t) & is.finite(lo) & is.finite(hi)
    if (!any(ok)) next
    hit <- (lo <= t & t <= hi)[ok]; sets <- x$set[ok]
    cl <- split(hit, sets)
    bs <- replicate(2000, { k <- sample(names(cl), length(cl), replace = TRUE); mean(unlist(cl[k])) })
    cov_rows[[length(cov_rows) + 1L]] <- data.frame(unit = x$unit[1], n = x$n[1], estimator = x$estimator[1], truth = truth, interval = ci,
      intervals = sum(ok), sets = length(cl), coverage = mean(hit), lo = stats::quantile(bs, 0.025, names = FALSE), hi = stats::quantile(bs, 0.975, names = FALSE),
      truth_above = mean((t > hi)[ok]), truth_below = mean((t < lo)[ok]),
      adequate = mean(hit) >= 0.93 & stats::quantile(bs, 0.025, names = FALSE) >= 0.91)
  }
}
cv <- do.call(rbind, cov_rows)
write.table(cv, paste0(prefix, "_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Bias, from the largest gene sets of each unit, by branch-length class within the unit
bias_rows <- list()
for (u in unique(iv$unit)) {
  x <- iv[iv$unit == u & iv$n == max(iv$n[iv$unit == u]), ]
  cls <- cut(x$length, unique(stats::quantile(unique(x$length), c(0, 1 / 3, 2 / 3, 1))), include.lowest = TRUE, labels = FALSE)
  if (length(unique(x$length)) < 3) cls <- match(x$length, sort(unique(x$length)))
  x$class <- cls
  for (e in unique(x$estimator)) for (k in sort(unique(x$class))) {
    y <- x[x$estimator == e & x$class == k, ]
    rb <- y$est / y$truth_pop - 1
    bias_rows[[length(bias_rows) + 1L]] <- data.frame(unit = u, estimator = e, length_class = k, min_length = min(y$length), max_length = max(y$length),
      n_sets = length(unique(y$set)), n = y$n[1], median_rel_bias = stats::median(rb), mean_rel_bias = mean(rb),
      median_rel_bias_sample = if (e == "rel") NA else stats::median(y$est / y$truth_sample - 1))
  }
}
bs <- do.call(rbind, bias_rows)
write.table(bs, paste0(prefix, "_bias.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# L4: ascertainment, t1 against t3 in L-A1-type units (genes lose t1); unmatched and matched
asc <- list()
for (u in units) {
  U <- load_unit(u)
  if (!any(grepl("t1", U$genes$dropped))) next
  t13 <- log(U$truth$omega_mean[U$truth$branch_id == "t1"] / U$truth$omega_mean[U$truth$branch_id == "t3"])
  x1 <- U$b[U$b$species_tree_branches == "t1", ]; x3 <- U$b[U$b$species_tree_branches == "t3", ]
  shared <- intersect(x1$orthogroup, x3$orthogroup)
  boot_lr <- function(a1, a3, matched) replicate(B, {
    if (matched) { k <- sample(shared, replace = TRUE); i1 <- match(k, a1$orthogroup); i3 <- match(k, a3$orthogroup) }
    else { i1 <- sample.int(nrow(a1), replace = TRUE); i3 <- sample.int(nrow(a3), replace = TRUE) }
    log(ratio(a1$NdN_free[i1], a1$N_free[i1], a1$SdS_free[i1], a1$S_free[i1]) / ratio(a3$NdN_free[i3], a3$N_free[i3], a3$SdS_free[i3], a3$S_free[i3]))
  })
  for (matched in c(FALSE, TRUE)) {
    a1 <- if (matched) x1[x1$orthogroup %in% shared, ] else x1; a3 <- if (matched) x3[x3$orthogroup %in% shared, ] else x3
    lr <- log(ratio(a1$NdN_free, a1$N_free, a1$SdS_free, a1$S_free) / ratio(a3$NdN_free, a3$N_free, a3$SdS_free, a3$S_free))
    bo <- boot_lr(a1, a3, matched)
    asc[[length(asc) + 1L]] <- data.frame(unit = u, comparison = if (matched) "matched" else "unmatched", genes_t1 = nrow(a1), genes_t3 = nrow(a3),
      log_ratio = lr, truth = t13, bias = lr - t13, lo = stats::quantile(bo, 0.025, names = FALSE) - t13, hi = stats::quantile(bo, 0.975, names = FALSE) - t13,
      adequate = abs(lr - t13) <= 0.05)
  }
}
if (length(asc)) write.table(do.call(rbind, asc), paste0(prefix, "_ascertainment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# L5: trait regressions across branches; trait correlated 0.5 with the part of log length orthogonal to log truth
tr_rows <- list()
for (u in unique(iv$unit)) {
  x <- iv[iv$unit == u & iv$n == 100, ]
  if (!nrow(x)) next
  for (e in c("pool", "median", "rel")) {
    y <- x[x$estimator == e & is.finite(x$est) & x$est > 0, ]
    fp <- c(est = 0, est_adj = 0, truth = 0); tot <- 0
    for (g in unique(y$set)) {
      z <- y[y$set == g, ]
      if (nrow(z) < 6) next
      lt <- log(z$truth_pop); ll <- log(z$length)
      r <- stats::residuals(stats::lm(ll ~ lt)); if (stats::sd(r) == 0) next
      r <- (r - mean(r)) / stats::sd(r)
      for (k in seq_len(10000 %/% length(unique(y$set)) + 1L)) {
        trait <- 0.5 * r + sqrt(0.75) * stats::rnorm(nrow(z))
        p1 <- summary(stats::lm(log(z$est) ~ trait))$coefficients["trait", 4]
        p2 <- summary(stats::lm(log(z$est) ~ trait + log(z$ds)))$coefficients["trait", 4]
        p3 <- summary(stats::lm(lt ~ trait))$coefficients["trait", 4]
        fp <- fp + c(p1 < 0.05, p2 < 0.05, p3 < 0.05); tot <- tot + 1
      }
    }
    if (tot) tr_rows[[length(tr_rows) + 1L]] <- data.frame(unit = u, estimator = e, draws = tot, fp_est = fp[["est"]] / tot,
      fp_est_adj_ds = fp[["est_adj"]] / tot, fp_truth = fp[["truth"]] / tot)
  }
}
if (length(tr_rows)) write.table(do.call(rbind, tr_rows), paste0(prefix, "_trait.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Coverage, percentile interval, population truth:\n")
print(cv[cv$interval == "pct" & cv$truth == "truth_pop", c("unit", "n", "estimator", "intervals", "coverage", "lo", "hi", "truth_above", "adequate")], row.names = FALSE)
cat("\nBias:\n"); print(bs, row.names = FALSE)
if (length(asc)) { cat("\nAscertainment:\n"); print(do.call(rbind, asc), row.names = FALSE) }
if (length(tr_rows)) { cat("\nTrait false positives:\n"); print(do.call(rbind, tr_rows), row.names = FALSE) }
