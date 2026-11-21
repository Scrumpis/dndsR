# Summarize the lineage follow-up study (docs/dev/phylo-lineage-followup.md) against its criteria.
#   Rscript --vanilla phylo-lineage-followup-summary.R <follow-up root> <first study root> [<direct RELAX folder>]
# Reads <root>/units/<cell>_<tree>[-<settings>]/; writes <root>/summary/{size,diagnostics,relax_pairs,direct}.tsv.
# Base R only.
a <- commandArgs(TRUE); root <- a[1]; first <- a[2]; direct <- if (length(a) > 2) a[3] else NA
dir.create(file.path(root, "summary"), showWarnings = FALSE)
cp <- function(x, n) if (n == 0) c(NA, NA) else stats::binom.test(x, n)$conf.int[1:2]
first_verdict <- function(x, n) {
  ci <- cp(x, n); r <- x / n
  if (n == 0) "no data" else if (ci[2] <= 0.05) "calibrated" else if (ci[1] > 0.05) "anticonservative" else if (r <= 0.05) "consistent with nominal" else "above nominal, not significantly"
}
bradley <- function(x, n) {
  ci <- cp(x, n)
  if (n == 0) "no data" else if (ci[1] >= 0.025 && ci[2] <= 0.075) "within tolerance" else if (ci[1] > 0.075 || ci[2] < 0.025) "outside tolerance" else "inconclusive"
}
rate_row <- function(unit, test, p) {
  p <- p[!is.na(p)]; x <- sum(p <= 0.05); n <- length(p); ci <- cp(x, n)
  data.frame(unit = unit, test = test, n = n, positives = x, rate = if (n) x / n else NA, lo = ci[1], hi = ci[2],
             rate_01 = if (n) mean(p <= 0.01) else NA, verdict = first_verdict(x, n), bradley = bradley(x, n))
}
read_res <- function(d, m) {
  f <- file.path(d, m, sprintf("%s_results.tsv", m))
  if (!file.exists(f)) return(NULL)
  x <- utils::read.delim(f, stringsAsFactors = FALSE)
  x$qc_flags <- ifelse(is.na(x$qc_flags), "", as.character(x$qc_flags))
  x[x$status == "tested", ]
}
units <- list.dirs(file.path(root, "units"), recursive = FALSE)
units <- units[file.exists(file.path(units, "DONE"))]
size <- list(); diag <- list(); relax <- list()
for (d in units) {
  u <- basename(d); cell <- sub("_.*", "", u)
  if (cell == "D0") {
    r <- read_res(d, "branch_omega")
    ok <- !is.na(r$p_fg_omega_gt1)
    size[[length(size) + 1]] <- rate_row(u, "p_fg_omega_gt1 (one-sided, omega_fg = 1)", r$p_fg_omega_gt1)
    size[[length(size) + 1]] <- rate_row(u, "lrt_fg_omega_eq1 two-sided, chi-square 1", stats::pchisq(r$lrt_fg_omega_eq1, 1, lower.tail = FALSE))
    size[[length(size) + 1]] <- rate_row(u, "p_adj_fg_omega_gt1", r$p_adj_fg_omega_gt1)
    sig <- ok & r$p_fg_omega_gt1 <= 0.05
    diag[[length(diag) + 1]] <- data.frame(unit = u, n_tested = nrow(r), n_two_ratio = sum(ok),
      frac_omega_fg_gt1 = mean(r$omega_fg[ok] > 1), median_omega_fg = stats::median(r$omega_fg[ok]),
      mean_lrt_fg_eq1 = mean(r$lrt_fg_omega_eq1[ok]), var_lrt_fg_eq1 = stats::var(r$lrt_fg_omega_eq1[ok]),
      p_two_vs_m0_rate = mean(r$p_two_vs_m0[ok] <= 0.05),
      sig_with_codeml_warning = sum(sig & grepl("codeml_warning", r$qc_flags)),
      sig_with_lrt_negative = sum(sig & grepl("lrt_negative", r$qc_flags)),
      any_lrt_negative_fg = sum(grepl("lrt_negative:fg_omega_eq1", r$qc_flags)),
      excluded_two_ratio = sum(!ok))
  }
  if (cell %in% c("C0", "C0i")) {
    r <- read_res(d, "relax")
    size[[length(size) + 1]] <- rate_row(u, "relax p_value", r$p_value)
    relax[[u]] <- r
    diag[[length(diag) + 1]] <- data.frame(unit = u, n_tested = nrow(r), mean_lrt = mean(pmax(r$lrt, 0)), median_lrt = stats::median(pmax(r$lrt, 0)),
      q95_lrt = unname(stats::quantile(pmax(r$lrt, 0), 0.95)), negative_lrt = sum(r$lrt < 0), k_at_bound = sum(r$k_at_bound),
      sig_relaxed = sum(r$p_value <= 0.05 & r$k < 1), sig_intensified = sum(r$p_value <= 0.05 & r$k > 1))
  }
  if (cell == "A0") {
    r <- read_res(d, "branch_omega")
    size[[length(size) + 1]] <- rate_row(u, "p_free_vs_m0 (replicate)", r$p_free_vs_m0)
    f <- utils::read.delim(file.path(first, "units", "A0_T16m", "branch_omega", "branch_omega_results.tsv"))
    size[[length(size) + 1]] <- rate_row(u, "p_free_vs_m0 (first study + replicate, 600)", c(r$p_free_vs_m0, f$p_free_vs_m0))
    b <- utils::read.delim(file.path(d, "branch_omega", "branch_omega_branches.tsv"))
    nm <- tapply(grepl("+", b$species_tree_branches, fixed = TRUE), b$orthogroup, sum)
    r$n_merged <- as.integer(nm[r$orthogroup])
    for (k in sort(unique(r$n_merged))) size[[length(size) + 1]] <- rate_row(u, sprintf("p_free_vs_m0, %d merged branches", k), r$p_free_vs_m0[r$n_merged == k])
    diag[[length(diag) + 1]] <- data.frame(unit = u, n_tested = nrow(r), df = paste(unique(r$df_free_vs_m0), collapse = ","),
      mean_lrt_over_df = mean(r$lrt_free_vs_m0 / r$df_free_vs_m0), var_ratio = stats::var(r$lrt_free_vs_m0) / (2 * mean(r$df_free_vs_m0)),
      codeml_warning_free = sum(grepl("codeml_warning:free", r$qc_flags)), starts_disagree_free = sum(grepl("starts_disagree:free", r$qc_flags)))
  }
}
out <- function(x, f) { # rows may have different columns (one set per cell); missing ones are NA
  cols <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(r) { for (c in setdiff(cols, names(r))) r[[c]] <- NA; r[cols] })
  utils::write.table(do.call(rbind, x), file.path(root, "summary", f), sep = "\t", quote = FALSE, row.names = FALSE)
}
out(size, "size.tsv"); out(diag, "diagnostics.tsv")
# R2 through dndsR: the same genes (same seed) with 1 and 10 starting points.
one <- relax[["C0_T16"]]; ten <- relax[["C0_T16-starting_points_10"]]
if (!is.null(one) && !is.null(ten)) {
  m <- merge(one[, c("orthogroup", "lrt", "p_value", "lnl_null", "lnl_alternative")], ten[, c("orthogroup", "lrt", "p_value", "lnl_null", "lnl_alternative")],
             by = "orthogroup", suffixes = c("_1", "_10"))
  m$d_lnl_null <- m$lnl_null_10 - m$lnl_null_1; m$d_lnl_alt <- m$lnl_alternative_10 - m$lnl_alternative_1
  utils::write.table(m, file.path(root, "summary", "relax_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  s1 <- m$p_value_1 <= 0.05
  cat(sprintf("R2 (dndsR, new genes): significant with 1 start %d, with 10 %d; of the 1-start positives, still significant %d; null lnL improved by >0.01 in %d genes, alternative in %d; median LRT change among 1-start positives %.3f\n",
              sum(s1), sum(m$p_value_10 <= 0.05), sum(s1 & m$p_value_10 <= 0.05), sum(m$d_lnl_null > 0.01), sum(m$d_lnl_alt > 0.01),
              stats::median(m$lrt_10[s1] - m$lrt_1[s1])))
}
# R1 and R2 on the first study's genes, run directly through HyPhy.
if (!is.na(direct)) {
  ref <- utils::read.delim(file.path(first, "units", "C0_T16", "relax", "relax_results.tsv"))
  rd <- function(mode, s) {
    f <- file.path(direct, sprintf("direct_%s_seed%d.tsv", mode, s))
    if (!file.exists(f)) return(NULL)
    x <- utils::read.delim(f); x[match(ref$orthogroup, x$orthogroup), ]
  }
  rows <- list()
  add <- function(label, lrt, p, extra = list()) {
    r <- cbind(rate_row(label, "relax p_value", p), mean_lrt = mean(pmax(lrt, 0), na.rm = TRUE),
      n_missing = sum(is.na(p)), max_abs_lrt_diff_vs_dndsR = max(abs(lrt - ref$lrt), na.rm = TRUE),
      n_lrt_diff_gt_0.01 = sum(abs(lrt - ref$lrt) > 0.01, na.rm = TRUE))
    rows[[length(rows) + 1]] <<- if (length(extra)) cbind(r, as.data.frame(extra)) else r
  }
  add("dndsR, first study (seed 1)", ref$lrt, ref$p_value)
  x <- rd("dndsr", 1); if (!is.null(x)) add("direct, dndsR's inputs, seed 1 (R1)", x$lrt, x$p_value)
  hs <- Filter(Negate(is.null), lapply(1:3, function(s) rd("hand", s)))
  for (i in seq_along(hs)) add(sprintf("direct, hand-written inputs, seed %d", i), hs[[i]]$lrt, hs[[i]]$p_value)
  if (length(hs) == 3) {
    nul <- do.call(pmax, c(lapply(hs, `[[`, "lnl_null"), na.rm = TRUE))
    alt <- do.call(pmax, c(lapply(hs, `[[`, "lnl_alternative"), na.rm = TRUE))
    lrt <- pmax(2 * (alt - nul), 0)
    spread_null <- do.call(pmax, c(lapply(hs, `[[`, "lnl_null"), na.rm = TRUE)) - do.call(pmin, c(lapply(hs, `[[`, "lnl_null"), na.rm = TRUE))
    spread_alt <- do.call(pmax, c(lapply(hs, `[[`, "lnl_alternative"), na.rm = TRUE)) - do.call(pmin, c(lapply(hs, `[[`, "lnl_alternative"), na.rm = TRUE))
    add("direct, hand-written inputs, best of seeds 1-3 per model (R2)", lrt, stats::pchisq(lrt, 1, lower.tail = FALSE),
        list(genes_null_spread_gt_0.01 = sum(spread_null > 0.01, na.rm = TRUE), genes_alt_spread_gt_0.01 = sum(spread_alt > 0.01, na.rm = TRUE),
             median_null_spread = stats::median(spread_null, na.rm = TRUE), median_alt_spread = stats::median(spread_alt, na.rm = TRUE)))
  }
  if (length(rows)) {
    cols <- unique(unlist(lapply(rows, names)))
    rows <- lapply(rows, function(r) { for (c in setdiff(cols, names(r))) r[[c]] <- NA; r[cols] })
    out(rows, "direct.tsv")
  }
}
cat("units:", length(units), "\n")
