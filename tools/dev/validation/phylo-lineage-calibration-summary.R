# Summarize the lineage calibration against the criteria fixed in
# docs/dev/phylo-lineage-calibration.md. Reads <root>/units/<cell>_<tree>/; writes
# <root>/summary/{size,power,estimation,pooled,timing}.tsv. Base R only; the pooled omega is
# recomputed here, not taken from phylo_summarize.
#   Rscript --vanilla phylo-lineage-calibration-summary.R <results root>
root <- commandArgs(TRUE)[1]
units <- list.dirs(file.path(root, "units"), recursive = FALSE)
units <- units[file.exists(file.path(units, "DONE"))]
dir.create(file.path(root, "summary"), showWarnings = FALSE)
cp <- function(x, n) if (n == 0) c(NA, NA) else stats::binom.test(x, n)$conf.int[1:2]
verdict <- function(x, n) {
  ci <- cp(x, n); r <- x / n
  if (n == 0) "no data" else if (ci[2] <= 0.05) "calibrated" else if (ci[1] > 0.05) "anticonservative" else if (r <= 0.05) "consistent with nominal" else "above nominal, not significantly"
}
rate_row <- function(unit, test, p, alpha = 0.05) {
  p <- p[!is.na(p)]; x <- sum(p <= alpha); n <- length(p); ci <- cp(x, n)
  data.frame(unit = unit, test = test, n = n, positives = x, rate = if (n) x / n else NA, lo = ci[1], hi = ci[2], verdict = verdict(x, n))
}
read <- function(d, m, what) {
  f <- file.path(d, m, sprintf("%s_%s.tsv", m, what))
  if (file.exists(f)) utils::read.delim(f, stringsAsFactors = FALSE) else NULL
}
size <- list(); power <- list(); est <- list(); pooled <- list(); timing <- list()
for (d in units) {
  unit <- basename(d); cell <- sub("_.*", "", unit)
  timing[[unit]] <- utils::read.delim(file.path(d, "timing.tsv"))
  truth <- utils::read.delim(file.path(d, "sim", "branch_truth.tsv"), stringsAsFactors = FALSE)
  if (substr(cell, 1, 1) == "A") {
    r <- read(d, "branch_omega", "results"); b <- read(d, "branch_omega", "branches")
    t <- r[r$status == "tested", ]
    if (cell %in% c("A0", "A3")) {
      for (p in c("p_two_vs_m0", "p_free_vs_m0", "p_fg_omega_gt1", "p_adj_two_vs_m0")) size[[length(size) + 1]] <- rate_row(unit, p, t[[p]])
    }
    if (cell == "A1") power[[length(power) + 1]] <- rate_row(unit, "p_two_vs_m0 (fg 0.6 vs bg 0.2)", t$p_two_vs_m0)
    if (cell == "A4") power[[length(power) + 1]] <- rate_row(unit, "p_fg_omega_gt1 (fg 1.5)", t$p_fg_omega_gt1)
    # Per-branch free-ratio omega against the truth, by dndsR's reliability flag.
    b <- b[b$orthogroup %in% t$orthogroup & b$in_species_tree & !grepl("+", b$species_tree_branches, fixed = TRUE), ]
    b$truth <- truth$omega_mean[match(b$species_tree_branches, truth$branch_id)]
    ok <- is.finite(b$omega_free) & b$omega_free > 0 & is.finite(b$truth)
    err <- abs(log2(b$omega_free[ok] / b$truth[ok]))
    rel <- b$reliable_free[ok] %in% TRUE
    est[[length(est) + 1]] <- data.frame(unit = unit, branches = sum(ok), frac_reliable = mean(rel),
      mae_log2_all = stats::median(err), mae_log2_reliable = stats::median(err[rel]), mae_log2_flagged = if (any(!rel)) stats::median(err[!rel]) else NA,
      frac_within_2fold_reliable = mean(err[rel] <= 1))
    if (cell %in% c("A1", "A4")) {
      tr_fg <- if (cell == "A1") 0.6 else 1.5
      est[[length(est)]]$two_ratio_fg_mae_log2 <- stats::median(abs(log2(t$omega_fg / tr_fg)), na.rm = TRUE)
      est[[length(est)]]$two_ratio_bg_mae_log2 <- stats::median(abs(log2(t$omega_bg / 0.2)), na.rm = TRUE)
    }
    if (cell == "A2") {
      bb <- b[ok & !grepl("high_ds", b$flags), ]
      for (id in unique(bb$species_tree_branches)) {
        x <- bb[bb$species_tree_branches == id, ]
        po <- (sum(x$NdN_free) / sum(x$N_free)) / (sum(x$SdS_free) / sum(x$S_free))
        xr <- x[x$reliable_free %in% TRUE, ]
        pooled[[length(pooled) + 1]] <- data.frame(unit = unit, branch_id = id, truth = x$truth[1], n = nrow(x), pooled = po,
          median_reliable = if (nrow(xr)) stats::median(xr$omega_free) else NA, mean_all = mean(x$omega_free[x$omega_free < 998]))
      }
    }
  }
  if (substr(cell, 1, 1) == "B") {
    r <- read(d, "absrel", "results"); b <- read(d, "absrel", "branches")
    t <- r[r$status == "tested", ]; b <- b[b$orthogroup %in% t$orthogroup & !is.na(b$p_raw), ]
    sel <- truth$branch_id[truth$selected]
    is_sel <- b$species_tree_branches %in% sel
    if (cell == "B0") {
      size[[length(size) + 1]] <- rate_row(unit, "absrel p_raw per branch", b$p_raw)
      size[[length(size) + 1]] <- rate_row(unit, "absrel p_value per gene (Holm min)", t$p_value)
      size[[length(size) + 1]] <- rate_row(unit, "absrel p_bh_lineage per branch", b$p_bh_lineage)
    } else {
      power[[length(power) + 1]] <- rate_row(unit, "absrel p_holm, selected branch", b$p_holm[is_sel])
      power[[length(power) + 1]] <- rate_row(unit, "absrel p_bh_lineage, selected branch", b$p_bh_lineage[is_sel])
      size[[length(size) + 1]] <- rate_row(unit, "absrel p_raw, unselected branches", b$p_raw[!is_sel])
    }
  }
  if (substr(cell, 1, 1) == "C") {
    r <- read(d, "relax", "results"); t <- r[r$status == "tested", ]
    if (cell == "C0") size[[length(size) + 1]] <- rate_row(unit, "relax p_value", t$p_value) else {
      want <- if (cell == "C1") "relaxation" else "intensification"
      power[[length(power) + 1]] <- rate_row(unit, sprintf("relax p_value (k = %s)", if (cell == "C1") "0.5" else "2"), t$p_value)
      sig <- t[t$p_value <= 0.05, ]
      power[[length(power)]]$correct_direction <- if (nrow(sig)) mean(sig$direction == want) else NA
    }
  }
}
bind <- function(x) {
  cols <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(r) { for (c in setdiff(cols, names(r))) r[[c]] <- NA; r[cols] }))
}
w <- function(x, f) if (length(x)) utils::write.table(bind(x), file.path(root, "summary", f), sep = "\t", quote = FALSE, row.names = FALSE)
w(size, "size.tsv"); w(power, "power.tsv"); w(est, "estimation.tsv"); w(timing, "timing.tsv")
if (length(pooled)) w(pooled, "pooled.tsv")
cat("units:", length(units), "\n")
