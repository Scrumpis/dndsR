# Summarizes the calibration pilot as fixed in docs/dev/phylo-calibration-design.md.
#   Rscript --vanilla phylo-calibration-summary.R <results dir> <out dir>
# Writes rates.tsv (per run and cell), null_distribution.tsv, paired.tsv, genomewide.tsv,
# qc_flags.tsv, and pvalue_histograms.svg.
a <- commandArgs(TRUE); res_dir <- a[1]; out <- a[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
runs <- list.dirs(res_dir, recursive = FALSE, full.names = FALSE)
read_run <- function(r) {
  f <- file.path(res_dir, r, paste0(r, "_results.tsv"))
  if (!file.exists(f)) return(NULL)
  x <- utils::read.delim(f, na.strings = "NA")
  t <- utils::read.delim(file.path(res_dir, r, "truth.tsv"))
  x$alternative <- t$alternative[match(x$orthogroup, t$gene)]
  x$run <- r
  x
}
dat <- Filter(Negate(is.null), lapply(runs, read_run))
names(dat) <- vapply(dat, function(x) x$run[1], character(1))
cp <- function(k, n) if (n > 0) stats::binom.test(k, n)$conf.int else c(NA, NA)
fmt <- function(k, n) { ci <- cp(k, n); sprintf("%.3f (%.3f-%.3f)", k / n, ci[1], ci[2]) }

# Rates per run (each run is one cell under one setting). Never pooled across cells.
rates <- do.call(rbind, lapply(dat, function(x) {
  t <- x$status == "tested"
  n <- sum(t); s05 <- sum(x$p_value[t] < 0.05); s01 <- sum(x$p_value[t] < 0.01)
  ci05 <- cp(s05, n); ci01 <- cp(s01, n)
  data.frame(run = x$run[1], kind = if (all(x$alternative, na.rm = TRUE)) "power" else if (any(x$alternative, na.rm = TRUE)) "mixed" else "null",
             genes = nrow(x), tested = n, excluded = sum(x$status == "excluded"), failed = sum(x$status == "failed"),
             rate_05 = s05 / n, lo_05 = ci05[1], hi_05 = ci05[2], rate_01 = s01 / n, lo_01 = ci01[1], hi_01 = ci01[2],
             lrt_zero = mean(x$lrt[t] < 1e-8), refits = sum(grepl("hyphy_refit", x$qc_flags[t])),
             flag = if (all(!x$alternative, na.rm = TRUE)) {
               if (ci05[1] > 0.05) "inflated (lower bound > 0.05)" else if (ci05[2] <= 0.075) "within 0.075" else "uncertain"
             } else "")
}))
utils::write.table(rates, file.path(out, "rates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Null distribution: share of LRT = 0, and positive LRTs against chi-square with 2 df.
nulls <- rates$run[rates$kind == "null"]
nd <- do.call(rbind, lapply(nulls, function(r) {
  x <- dat[[r]]; l <- x$lrt[x$status == "tested"]
  pos <- l[l > 1e-8]
  ks <- if (length(pos) > 5) suppressWarnings(stats::ks.test(pos, "pchisq", 2)$p.value) else NA
  data.frame(run = r, n = length(l), share_lrt_zero = mean(l <= 1e-8), mean_positive_lrt = mean(pos),
             ks_positive_vs_chisq2_p = ks,
             share_positive_above_chisq2_q95 = mean(pos > stats::qchisq(0.95, 2)))
}))
utils::write.table(nd, file.path(out, "null_distribution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# P-value histograms with the mass at exactly 0.5 shown separately.
grDevices::svg(file.path(out, "pvalue_histograms.svg"), width = 10, height = 2.2 * ceiling(length(nulls) / 4))
graphics::par(mfrow = c(ceiling(length(nulls) / 4), 4), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.5, 0))
for (r in nulls) {
  p <- dat[[r]]$p_value[dat[[r]]$status == "tested"]
  graphics::hist(p[p < 0.4999], breaks = seq(0, 0.5, 0.025), main = sprintf("%s (p = 0.5: %d of %d)", r, sum(p >= 0.4999), length(p)),
                 xlab = "p-value below 0.5", cex.main = 0.75, col = "grey80", border = "white")
  graphics::abline(h = length(p) * 0.5 / 20, lty = 2) # expected per bin if positive LRTs follow chi-square 2
}
grDevices::dev.off()

# Paired comparisons: same genes, another strategy, setting, or processor.
pairs <- list(c("S0_N0_L50", "S1_N0_L50"), c("S0_N0_L50", "S2_N0_L50"), c("S0_A10_p05", "S1_A10_p05"), c("S0_A10_p05", "S2_A10_p05"),
              c("S0_N0_L50", "MH_N0_L50"), c("S0_M1_mnp", "MH_M1_mnp"), c("S0_N0_L50", "E_N0_L50"), c("S0_M3_misalign", "E_M3_misalign"),
              c("S0_N0_L50", "AVX2_N0_L50"))
timing <- if (file.exists(file.path(res_dir, "timing.tsv"))) utils::read.delim(file.path(res_dir, "timing.tsv"), header = FALSE, col.names = c("run", "seconds", "threads", "host")) else NULL
paired <- do.call(rbind, lapply(pairs, function(pr) {
  if (!all(pr %in% names(dat))) return(NULL)
  a <- dat[[pr[1]]]; b <- dat[[pr[2]]]
  m <- merge(a[a$status == "tested", ], b[b$status == "tested", ], by = "orthogroup", suffixes = c(".a", ".b"))
  call_a <- m$p_value.a < 0.05; call_b <- m$p_value.b < 0.05
  best <- if ("lrt_best_fits.b" %in% names(m)) m$lrt_best_fits.b else NA
  call_best <- if (!all(is.na(best))) 0.5 * stats::pchisq(best, 2, lower.tail = FALSE) < 0.05 else NA
  secs <- function(r) if (is.null(timing)) NA else timing$seconds[timing$run == r][1]
  data.frame(reference = pr[1], comparison = pr[2], genes = nrow(m),
             rate_reference = mean(call_a), rate_comparison = mean(call_b),
             calls_changed = sum(call_a != call_b), lrt_identical = sum(abs(m$lrt.a - m$lrt.b) < 1e-12),
             max_abs_lrt_difference = max(abs(m$lrt.a - m$lrt.b)),
             rate_best_fits = if (all(is.na(call_best))) NA else mean(call_best),
             calls_changed_best_fits = if (all(is.na(call_best))) NA else sum(call_best != call_a),
             primary_alt_not_best = if ("primary_unconstrained_is_best.b" %in% names(m)) sum(!m$primary_unconstrained_is_best.b) else NA,
             primary_null_not_best = if ("primary_null_is_best.b" %in% names(m)) sum(!m$primary_null_is_best.b) else NA,
             seconds_reference = secs(pr[1]), seconds_comparison = secs(pr[2]))
}))
if (!is.null(paired)) utils::write.table(paired, file.path(out, "paired.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Genome-wide: BH within the run.
if ("S0_genomewide" %in% names(dat)) {
  g <- dat[["S0_genomewide"]]; g <- g[g$status == "tested", ]
  gw <- do.call(rbind, lapply(c(0.05, 0.1), function(q) {
    d <- g$p_adj < q; fd <- sum(d & !g$alternative)
    data.frame(q = q, tested = nrow(g), true_alternatives = sum(g$alternative), discoveries = sum(d), false_discoveries = fd,
               fdp = if (sum(d)) fd / sum(d) else 0, power = sum(d & g$alternative) / sum(g$alternative))
  }))
  utils::write.table(gw, file.path(out, "genomewide.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# QC flags among false positives and true negatives in misspecified null cells.
flag_rows <- list()
for (r in intersect(c("S0_M3_misalign", "S0_M4_indel", "S0_N0_L50", "S0_M1_mnp", "S0_M5_discordant"), names(dat))) {
  x <- dat[[r]]; x <- x[x$status == "tested", ]
  pos <- x$p_value < 0.05
  flags <- unique(sub(":.*", "", unlist(strsplit(ifelse(is.na(x$qc_flags), "", x$qc_flags), ";"))))
  for (f in setdiff(flags, "")) {
    has <- grepl(paste0("(^|;)", f), ifelse(is.na(x$qc_flags), "", x$qc_flags))
    flag_rows[[length(flag_rows) + 1L]] <- data.frame(run = r, flag = f, false_positives = sum(pos), flagged_false_positives = sum(has & pos),
                                                      true_negatives = sum(!pos), flagged_true_negatives = sum(has & !pos),
                                                      rate_if_flagged = if (sum(has)) mean(pos[has]) else NA,
                                                      rate_if_not_flagged = if (sum(!has)) mean(pos[!has]) else NA)
  }
}
if (length(flag_rows)) utils::write.table(do.call(rbind, flag_rows), file.path(out, "qc_flags.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(rates[, c("run", "tested", "failed", "rate_05", "lo_05", "hi_05", "rate_01", "lrt_zero", "flag")], digits = 3, row.names = FALSE)
