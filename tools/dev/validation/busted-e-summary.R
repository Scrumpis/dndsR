# BUSTED-E pilot summary (docs/dev/busted-e-pilot.md; measures and criteria fixed there).
#   Rscript --vanilla busted-e-summary.R <busted-e root> <study-B3 root> <out prefix>
# BUSTED and masked BUSTED come from study B (study-B3, unchanged); BUSTED-E from this pilot, on the
# same genes. Writes <prefix>_{size,mixture,cost,criteria}.tsv.
a <- commandArgs(TRUE); root <- a[1]; b3 <- a[2]; prefix <- a[3]
cp <- function(x, n) if (n == 0) c(NA, NA) else stats::binom.test(x, n)$conf.int[1:2]
category <- function(x, n) { ci <- cp(x, n); if (n == 0) return("no data")
  if (ci[2] <= 0.05) "calibrated" else if (ci[1] >= 0.025 && ci[2] <= 0.075) "within tolerance" else if (ci[1] > 0.05) "anticonservative" else "inconclusive" }
rdres <- function(f) { r <- utils::read.delim(f, stringsAsFactors = FALSE); r[r$status == "tested", ] }
arm_file <- function(unit, src, arm) if (arm == "busted_e") file.path(root, "units", unit, "busted_e", "busted_e_results.tsv") else
  file.path(b3, "units", src, arm, sprintf("%s_results.tsv", arm))
units <- utils::read.delim(file.path(root, "src", "tools", "dev", "validation", "busted-e-pilot-units.tsv"), stringsAsFactors = FALSE)
size <- list(); mix <- list(); cost <- list(); crit <- list()
for (i in seq_len(nrow(units))) {
  u <- units$name[i]; src <- units$source_unit[i]
  if (!file.exists(file.path(root, "units", u, "DONE"))) next
  e <- rdres(arm_file(u, src, "busted_e")); keep <- e$orthogroup
  arms <- list(busted = rdres(arm_file(u, src, "busted")), busted_mask = rdres(arm_file(u, src, "busted_mask")), busted_e = e)
  arms <- lapply(arms, function(r) r[r$orthogroup %in% keep, ])
  if (u != "E-MIX") {
    for (k in names(arms)) { p <- arms[[k]]$p_value; x <- sum(p <= 0.05, na.rm = TRUE); n <- sum(!is.na(p)); ci <- cp(x, n)
      size[[length(size) + 1L]] <- data.frame(unit = u, arm = k, genes = n, positives = x, rate = x / n, lo = ci[1], hi = ci[2], category = category(x, n)) }
  } else {
    genes <- utils::read.delim(file.path(root, "units", u, "genes.tsv"), stringsAsFactors = FALSE); sel <- genes$gene[genes$selected]
    # p_adj of masked and unmasked BUSTED were corrected over the same 1,000 genes; BUSTED-E's over its own run (the same genes).
    calls <- lapply(arms, function(r) stats::setNames(r$p_adj <= 0.05, r$orthogroup))
    ogs <- keep; ntrue <- sum(ogs %in% sel)
    for (k in names(arms)) { d <- calls[[k]][ogs]; d[is.na(d)] <- FALSE; D <- sum(d); F <- sum(d & !ogs %in% sel); ci <- if (D) cp(F, D) else c(NA, NA)
      mix[[k]] <- data.frame(arm = k, discoveries = D, false = F, fdp = if (D) F / D else NA, fdp_lo = ci[1], fdp_hi = ci[2], power = sum(d & ogs %in% sel) / ntrue, true_genes = ntrue) }
    # power difference BUSTED-E minus masked BUSTED, gene bootstrap (2,000, seed 20261122, as phylo-study-B-exclude.R)
    set.seed(20261122)
    pw <- function(k, s) { d <- calls[[k]][s]; d[is.na(d)] <- FALSE; sum(d & s %in% sel) / sum(s %in% sel) }
    bs <- replicate(2000, { s <- sample(ogs, replace = TRUE); pw("busted_e", s) - pw("busted_mask", s) })
    mix[["busted_e"]]$power_diff_vs_masked <- mix[["busted_e"]]$power - mix[["busted_mask"]]$power
    mix[["busted_e"]]$diff_lo <- stats::quantile(bs, 0.025, names = FALSE); mix[["busted_e"]]$diff_hi <- stats::quantile(bs, 0.975, names = FALSE)
  }
  pv <- utils::read.delim(file.path(root, "units", u, "provenance.tsv"), stringsAsFactors = FALSE)
  cost[[u]] <- data.frame(unit = u, genes = length(keep), seconds = pv$seconds, cpus = pv$cpus, core_s_per_gene = pv$seconds * pv$cpus / length(keep))
}
size <- do.call(rbind, size); cost <- do.call(rbind, cost)
mix <- if (length(mix)) do.call(rbind, lapply(mix, function(x) { for (c in c("power_diff_vs_masked", "diff_lo", "diff_hi")) if (is.null(x[[c]])) x[[c]] <- NA; x })) else NULL
add <- function(c, ok, detail) crit[[length(crit) + 1L]] <<- data.frame(criterion = c, pass = ok, detail = detail)
if (!is.null(size)) {
  e1 <- size[size$unit == "E-E1" & size$arm == "busted_e", ]; e0 <- size[size$unit == "E-E0" & size$arm == "busted_e", ]
  if (nrow(e1)) add("1 gate: E-E1 upper bound <= 7.5%", e1$hi <= 0.075, sprintf("%.1f%% (%.1f-%.1f)", 100 * e1$rate, 100 * e1$lo, 100 * e1$hi))
  if (nrow(e0)) add("1 gate: E-E0 not anticonservative", e0$lo <= 0.05, sprintf("%.1f%% (%.1f-%.1f)", 100 * e0$rate, 100 * e0$lo, 100 * e0$hi))
}
if (!is.null(mix)) {
  x <- mix[mix$arm == "busted_e", ]
  add("2 FDR: FDP upper bound <= 7.5%", !is.na(x$fdp_hi) && x$fdp_hi <= 0.075, sprintf("FDP %.3f (%.3f-%.3f), %d discoveries", x$fdp, x$fdp_lo, x$fdp_hi, x$discoveries))
  add("3 power above masked BUSTED, interval excluding 0", x$diff_lo > 0, sprintf("%.3f against %.3f; difference %.3f (%.3f to %.3f)", x$power, mix$power[mix$arm == "busted_mask"], x$power_diff_vs_masked, x$diff_lo, x$diff_hi))
}
crit <- if (length(crit)) do.call(rbind, crit) else NULL
w <- function(x, what) if (!is.null(x)) utils::write.table(x, paste0(prefix, "_", what, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w(size, "size"); w(mix, "mixture"); w(cost, "cost"); w(crit, "criteria")
for (x in list(size, mix, cost, crit)) if (!is.null(x)) print(x, digits = 3, row.names = FALSE)
