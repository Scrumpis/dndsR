# Summaries for stage 2 (docs/dev/phylo-calibration-stage2.md).
#   Rscript --vanilla phylo-stage2-summary.R <results root> [out dir]
# Writes units.tsv, paired.tsv, starts.tsv, starts_summary.tsv, filtering.tsv, cost.tsv,
# and incomplete.tsv (units not finished, with the reason). Every unit in the table appears,
# finished or not; failed and excluded fits are counted by reason, never dropped.
a <- commandArgs(TRUE); root <- a[1]; out <- if (length(a) > 1) a[2] else file.path(root, "summary")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
src <- jsonlite::read_json(file.path(root, "freeze.json"))$source
v <- file.path(src, "tools/dev/validation")
sc <- utils::read.delim(file.path(v, "phylo-stage2-scenarios.tsv"), stringsAsFactors = FALSE)
un <- utils::read.delim(file.path(v, "phylo-stage2-units.tsv"), stringsAsFactors = FALSE)
cp <- function(x, n) if (n == 0) c(NA, NA) else c(if (x == 0) 0 else stats::qbeta(0.025, x, n - x + 1),
                                                  if (x == n) 1 else stats::qbeta(0.975, x + 1, n - x))
counts <- function(v) if (!length(v)) "" else paste(sprintf("%s:%d", names(table(v)), as.integer(table(v))), collapse = "; ")

read_unit <- function(u) {
  d <- file.path(root, "units", u)
  f <- file.path(d, paste0(u, "_results.tsv"))
  meta <- if (file.exists(file.path(d, "meta.tsv"))) utils::read.delim(file.path(d, "meta.tsv"), stringsAsFactors = FALSE) else NULL
  res <- if (file.exists(f)) utils::read.delim(f, stringsAsFactors = FALSE) else NULL
  list(res = res, meta = meta, done = file.exists(file.path(d, "DONE")))
}
dat <- lapply(stats::setNames(un$unit, un$unit), read_unit)

# Per unit.
rows <- lapply(seq_len(nrow(un)), function(i) {
  u <- un[i, ]; s <- sc[sc$scenario == u$scenario, ]; d <- dat[[u$unit]]
  base <- data.frame(unit = u$unit, block = u$block, scenario = u$scenario, kind = s$kind, config = u$config,
                     filter = u$filter, seeds = u$seeds, simulated = s$n, finished = d$done, stringsAsFactors = FALSE)
  if (is.null(d$res)) return(cbind(base, tested = NA, excluded = NA, failed = NA, excluded_reasons = "", failed_reasons = "",
                                   positives_05 = NA, rate_05 = NA, lo_05 = NA, hi_05 = NA, rate_01 = NA, lo_01 = NA, hi_01 = NA,
                                   bh_q05 = NA, lrt_zero = NA, negative_lrt = NA, refits = NA, seconds = NA, threads = NA,
                                   core_seconds_per_fit = NA, cpu = NA))
  r <- d$res; t <- r[r$status == "tested", ]
  x05 <- sum(t$p_value < 0.05); x01 <- sum(t$p_value < 0.01); n <- nrow(t)
  m <- d$meta
  fits <- nrow(r) * length(strsplit(u$seeds, ",")[[1]])
  cbind(base, tested = n, excluded = sum(r$status == "excluded"), failed = sum(r$status == "failed"),
        excluded_reasons = counts(r$reason[r$status == "excluded"]), failed_reasons = counts(r$reason[r$status == "failed"]),
        positives_05 = x05, rate_05 = x05 / n, lo_05 = cp(x05, n)[1], hi_05 = cp(x05, n)[2],
        rate_01 = x01 / n, lo_01 = cp(x01, n)[1], hi_01 = cp(x01, n)[2],
        bh_q05 = sum(!is.na(t$p_adj) & t$p_adj < 0.05), lrt_zero = mean(t$lrt <= 1e-8), negative_lrt = sum(t$lrt < 0),
        refits = sum(grepl("hyphy_refit_after_negative_lrt", t$qc_flags)),
        seconds = if (is.null(m)) NA else m$seconds, threads = if (is.null(m)) NA else m$threads,
        core_seconds_per_fit = if (is.null(m)) NA else m$seconds * m$threads / fits, cpu = if (is.null(m)) NA else m$cpu)
})
units <- do.call(rbind, rows)
utils::write.table(units, file.path(out, "units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
inc <- units[!units$finished, c("unit", "block", "scenario", "config")]
inc$reason <- vapply(inc$unit, function(u) { m <- dat[[u]]$meta
  if (is.null(m)) "not run" else sprintf("exit code %s", m$exit_code) }, "")
utils::write.table(inc, file.path(out, "incomplete.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("units:", nrow(units), " finished:", sum(units$finished), " incomplete:", nrow(inc), "\n")
cpus <- unique(stats::na.omit(units$cpu)); cat("processors:", paste(cpus, collapse = " | "), "\n")

# Paired comparisons on the same genes.
pair <- function(ua, ub, label) {
  A <- dat[[ua]]$res; B <- dat[[ub]]$res
  if (is.null(A) || is.null(B)) return(NULL)
  A <- A[A$status == "tested", c("orthogroup", "p_value")]; B <- B[B$status == "tested", c("orthogroup", "p_value")]
  m <- merge(A, B, by = "orthogroup", suffixes = c(".a", ".b"))
  ca <- m$p_value.a < 0.05; cb <- m$p_value.b < 0.05
  gain <- sum(!ca & cb); loss <- sum(ca & !cb)
  data.frame(comparison = label, a = ua, b = ub, genes = nrow(m), rate_a = mean(ca), rate_b = mean(cb),
             a_only = loss, b_only = gain, calls_changed = gain + loss, share_changed = (gain + loss) / nrow(m),
             mcnemar_p = if (gain + loss) stats::binom.test(gain, gain + loss)$p.value else 1, stringsAsFactors = FALSE)
}
id <- function(block, s, c, f = "none", seeds = "s1") paste(block, s, c, f, seeds, sep = "__")
pr <- list()
for (s in sc$scenario) {
  for (c in c("C1", "C2", "C3")) pr[[length(pr) + 1]] <- pair(id("config", s, "C0"), id("config", s, c), paste0("config ", c, " vs C0"))
  pr[[length(pr) + 1]] <- pair(id("config", s, "C1"), id("config", s, "C3"), "config C3 vs C1")
  pr[[length(pr) + 1]] <- pair(id("config", s, "C2"), id("config", s, "C3"), "config C3 vs C2")
}
for (s in c("N0_L50", "M1_mnp", "A10_p10", "A10_p20")) pr[[length(pr) + 1]] <- pair(id("config", s, "C1"), id("triple", s, "C4"), "C4 vs C1")
for (u in un$unit[un$block == "filter"]) {
  x <- un[un$unit == u, ]; pr[[length(pr) + 1]] <- pair(id("config", x$scenario, x$config), u, paste0("filter ", x$filter, " vs none, ", x$config))
}
# Corruption against the clean simulation of the same genes (same seed: the first genes match).
clean <- c(R1_internal_At = "N0_L50", R2_terminal_At = "N0_L50", R3_frameshift_At = "N0_L50", R4_internal_Br = "N0_L50", A10_p10_R1 = "A10_p10")
for (s in names(clean)) for (c in c("C0", "C1", "C2", "C3"))
  pr[[length(pr) + 1]] <- pair(id("config", clean[[s]], c), id("config", s, c), paste0("corrupted vs clean, ", c))
paired <- do.call(rbind, pr)
utils::write.table(paired, file.path(out, "paired.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Starts: best fits over seeds {1}, {1,2}, {1,2,3}, read from the kept HyPhy JSON.
lnl <- function(j) {
  if (!file.exists(j)) return(c(NA, NA))
  x <- tryCatch(jsonlite::read_json(j), error = function(e) NULL)
  if (is.null(x)) return(c(NA, NA))
  unc <- x$fits[["Unconstrained model"]][["Log Likelihood"]]; con <- x$fits[["Constrained model"]][["Log Likelihood"]]
  c(if (is.null(unc)) NA else unc, if (is.null(con)) NA else con)
}
st <- list()
for (u in un$unit[un$block == "starts"]) {
  tarf <- file.path(root, "units", u, "hyphy_raw.tar.gz"); if (!file.exists(tarf)) next
  tmp <- tempfile(); dir.create(tmp); utils::untar(tarf, exdir = tmp)
  res <- dat[[u]]$res; res <- res[res$status == "tested", ]
  for (og in res$orthogroup) {
    d <- file.path(tmp, "hyphy", og)
    L <- rbind(lnl(file.path(d, "busted.json")), lnl(file.path(d, "start_seed2", "busted.json")), lnl(file.path(d, "start_seed3", "busted.json")))
    nullp <- ifelse(is.na(L[, 2]), L[, 1], L[, 2])
    best <- function(k) { lrt <- max(0, 2 * (max(L[1:k, 1]) - max(nullp[1:k]))); c(lrt, 0.5 * stats::pchisq(lrt, 2, lower.tail = FALSE)) }
    b <- sapply(1:3, best)
    st[[length(st) + 1]] <- data.frame(unit = u, orthogroup = og, primary_lrt = res$lrt[res$orthogroup == og],
      lrt_1 = b[1, 1], lrt_2 = b[1, 2], lrt_3 = b[1, 3], p_1 = b[2, 1], p_2 = b[2, 2], p_3 = b[2, 3],
      alt_gain_2 = max(L[1:2, 1]) - L[1, 1], alt_gain_3 = max(L[, 1]) - L[1, 1],
      null_gain_2 = max(nullp[1:2]) - nullp[1], null_gain_3 = max(nullp) - nullp[1], stringsAsFactors = FALSE)
  }
  unlink(tmp, recursive = TRUE)
}
if (length(st)) {
  st <- do.call(rbind, st)
  utils::write.table(st, file.path(out, "starts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  ss <- do.call(rbind, lapply(split(st, st$unit), function(x) data.frame(unit = x$unit[1], genes = nrow(x),
    rate_1 = mean(x$p_1 < 0.05), rate_2 = mean(x$p_2 < 0.05), rate_3 = mean(x$p_3 < 0.05),
    changed_1_vs_3 = sum((x$p_1 < 0.05) != (x$p_3 < 0.05)), changed_2_vs_3 = sum((x$p_2 < 0.05) != (x$p_3 < 0.05)),
    improved_01_by_2 = mean(pmax(x$alt_gain_2, x$null_gain_2) > 0.1), improved_05_by_2 = mean(pmax(x$alt_gain_2, x$null_gain_2) > 0.5),
    improved_01_by_3 = mean(pmax(x$alt_gain_3, x$null_gain_3) > 0.1), improved_05_by_3 = mean(pmax(x$alt_gain_3, x$null_gain_3) > 0.5),
    # Seed 1 alone must reproduce the reported LRT (a check of this reconstruction).
    max_abs_lrt1_vs_primary = max(abs(pmax(x$primary_lrt, 0) - x$lrt_1)))))
  utils::write.table(ss, file.path(out, "starts_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Filtering: what was removed, and how often clean genes are flagged.
fl <- do.call(rbind, lapply(sc$scenario, function(s) {
  d <- file.path(root, "sim", s); if (!file.exists(file.path(d, "DONE"))) return(NULL)
  p <- utils::read.delim(file.path(d, "windows_profile.tsv"), stringsAsFactors = FALSE)
  g <- tapply(p$windows > 0 & p$taxon %in% c("At", "Al"), p$alignment, any)
  gb <- tapply(p$windows > 0, p$alignment, any)
  out <- data.frame(scenario = s, genes = length(g), flagged_foreground = mean(g), flagged_any = mean(gb))
  for (f in c("mask", "mask+gappy")) {
    lg <- utils::read.delim(file.path(d, paste0("filtered_", f), "filter_log.tsv"))
    out[[paste0(f, "_codons_masked_mean")]] <- mean(lg$codons_masked)
    out[[paste0(f, "_columns_dropped_mean")]] <- mean(lg$columns_dropped)
  }
  out
}))
utils::write.table(fl, file.path(out, "filtering.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Cost by configuration and block.
ok <- units[units$finished & is.finite(units$core_seconds_per_fit), ]
cost <- stats::aggregate(core_seconds_per_fit ~ config + kind + seeds, ok, stats::median)
utils::write.table(cost, file.path(out, "cost.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("summary written to", out, "\n")
