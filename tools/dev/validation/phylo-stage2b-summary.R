# Summaries for stage 2B (docs/dev/phylo-calibration-stage2b.md), with its decision criteria.
#   Rscript --vanilla phylo-stage2b-summary.R <results root> <out dir> <unit table> [unit table ...]
# Unit tables are names in the frozen tools/dev/validation. Writes, in <out dir>:
#   genes.tsv.gz    one row per unit and gene: calls with 1 and 3 starts, likelihoods, AICc,
#                   alignment checksum, and the unit the fit came from
#   units.tsv       per unit: denominators, failures by reason, rates with exact intervals,
#                   D1 class, likelihood consistency, cost
#   audit.tsv       D4: one start against three in the audited null-like cells, and the decision
#   noise.tsv       replicate discordance (C1, seeds 1-3 against 4-6), the D3 noise floor
#   paired.tsv      paired comparisons: rates, discordance, McNemar, Newcombe interval for the
#                   difference, D3 and the non-inferiority class (margin 5 points)
#   rules.tsv       D8: C1, C3, C1 and C3, model averaging with 1 and 3 starts; mean BUSTED-E weight
#   mixture.tsv     D9: the genome-like mixture after BH
#   filtering.tsv   positions removed and genes changed by each filter
#   starts.tsv      per three-start unit: calls changed and fits improved
#   cost.tsv, incomplete.tsv
# Calls use p < 0.05. "p3" combines the best alternative and best null fits over the unit's
# starts (dndsR's p_best_fits, never the smallest p); "p1" is the first start as reported.
a <- commandArgs(TRUE); root <- a[1]; out <- a[2]; tables <- a[-(1:2)]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
src <- jsonlite::read_json(file.path(root, "freeze.json"))$source
v <- file.path(src, "tools/dev/validation")
sc <- utils::read.delim(file.path(v, "phylo-stage2b-scenarios.tsv"), stringsAsFactors = FALSE)
un <- do.call(rbind, lapply(tables, function(t) utils::read.delim(file.path(v, t), stringsAsFactors = FALSE)))
ALPHA <- 0.05; MARGIN <- 0.05; FG <- c("At", "Al")
nulllike <- c("null", "misspec", "corrupt")
wr <- function(x, f) utils::write.table(x, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE)
cp <- function(x, n) if (is.na(n) || n == 0) c(NA, NA) else c(if (x == 0) 0 else stats::qbeta(0.025, x, n - x + 1),
                                                                  if (x == n) 1 else stats::qbeta(0.975, x + 1, n - x))
d1 <- function(x, n) { b <- cp(x, n); if (is.na(b[1])) NA else if (b[2] <= 0.05) "calibrated" else if (b[1] > 0.05) "unsafe" else "inconclusive" }
counts <- function(v) if (!length(v)) "" else paste(sprintf("%s:%d", names(table(v)), as.integer(table(v))), collapse = "; ")

# --- Per-gene results, composed for units that fitted only the alignments their filter changed.
aicc_of <- function(tarf, ogs) {
  res <- data.frame(orthogroup = ogs, aicc_1 = NA_real_, aicc_best = NA_real_, stringsAsFactors = FALSE)
  if (!file.exists(tarf) || !length(ogs)) return(res)
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  utils::untar(tarf, exdir = tmp)
  for (i in seq_along(ogs)) {
    d <- file.path(tmp, "hyphy", ogs[i])
    js <- c(file.path(d, "busted.json"), Sys.glob(file.path(d, "start_seed*", "busted.json")))
    f <- t(vapply(js, function(j) {
      x <- tryCatch(jsonlite::read_json(j), error = function(e) NULL)
      u <- if (is.null(x)) NULL else x$fits[["Unconstrained model"]]
      if (is.null(u) || is.null(u[["AIC-c"]])) c(NA_real_, NA_real_) else c(u[["AIC-c"]], u[["Log Likelihood"]])
    }, numeric(2)))
    if (all(is.na(f[, 1]))) next
    res$aicc_1[i] <- f[1, 1]
    res$aicc_best[i] <- f[which.max(f[, 2]), 1]
  }
  res
}
keep <- c("orthogroup", "status", "reason", "alignment_md5", "lrt", "p_value", "lnl_unconstrained", "lnl_constrained",
          "n_starts", "lrt_best_fits", "p_best_fits", "lnl_unconstrained_best", "lnl_null_best", "qc_flags")
own <- list(); meta <- list()
for (u in un$unit) {
  d <- file.path(root, "units", u)
  f <- file.path(d, paste0(u, "_results.tsv"))
  meta[[u]] <- if (file.exists(file.path(d, "meta.tsv"))) utils::read.delim(file.path(d, "meta.tsv"), stringsAsFactors = FALSE) else NULL
  if (!file.exists(f)) { own[[u]] <- NULL; next }
  r <- utils::read.delim(f, stringsAsFactors = FALSE, colClasses = c(orthogroup = "character"))
  for (k in setdiff(keep, names(r))) r[[k]] <- NA
  r <- r[, keep]
  r$p1 <- r$p_value
  r$p3 <- ifelse(!is.na(r$p_best_fits), r$p_best_fits, r$p_value)
  r$reported_minus_lnl <- abs(pmax(r$lrt, 0) - pmax(0, 2 * (r$lnl_unconstrained - r$lnl_constrained)))
  r$refit <- grepl("hyphy_refit_after_negative_lrt", r$qc_flags)
  r <- merge(r, aicc_of(file.path(d, "hyphy_raw.tar.gz"), r$orthogroup[r$status == "tested"]), by = "orthogroup", all.x = TRUE)
  r$source_unit <- u
  own[[u]] <- r
}
done <- vapply(un$unit, function(u) file.exists(file.path(root, "units", u, "DONE")), TRUE)
genes <- list()
compose <- function(u) {
  if (!is.null(genes[[u]])) return(genes[[u]])
  x <- un[un$unit == u, ]; fit <- if (is.null(x$fit) || is.na(x$fit)) "all" else x$fit
  if (fit == "all") { g <- own[[u]] } else {
    b <- un$unit[un$scenario == x$scenario & un$config == x$config & un$seeds == x$seeds & un$filter == sub("^changed:", "", fit)]
    base <- if (length(b) == 1) compose(b) else NULL
    if (is.null(base) || !done[[u]]) g <- NULL else {
      mine <- own[[u]]
      g <- if (is.null(mine)) base else rbind(base[!base$orthogroup %in% mine$orthogroup, ], mine)
    }
  }
  genes[[u]] <<- g
  g
}
for (u in un$unit) compose(u)

truth <- do.call(rbind, lapply(sc$scenario, function(s) {
  d <- file.path(root, "sim", s); f <- file.path(d, "truth.tsv"); if (!file.exists(f)) return(NULL)
  t <- utils::read.delim(f, stringsAsFactors = FALSE)
  cl <- file.path(d, "corrupt_log.tsv")
  cg <- if (file.exists(cl)) utils::read.delim(cl, stringsAsFactors = FALSE) else data.frame(alignment = character(), taxon = character())
  cg$gene <- sub("\\.fasta$", "", cg$alignment)
  data.frame(scenario = s, orthogroup = t$gene, alternative = as.logical(t$alternative),
             corrupted = t$gene %in% cg$gene, corrupted_fg = t$gene %in% cg$gene[cg$taxon %in% FG], stringsAsFactors = FALSE)
}))

allg <- do.call(rbind, lapply(names(genes), function(u) if (is.null(genes[[u]])) NULL else
  cbind(unit = u, scenario = un$scenario[un$unit == u], genes[[u]], stringsAsFactors = FALSE)))
allg <- merge(allg, truth, by = c("scenario", "orthogroup"), all.x = TRUE)
gz <- gzfile(file.path(out, "genes.tsv.gz"), "w"); utils::write.table(allg, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

# --- Per unit.
tested <- function(u) { g <- genes[[u]]; if (is.null(g)) NULL else g[g$status == "tested", ] }
rows <- lapply(seq_len(nrow(un)), function(i) {
  x <- un[i, ]; s <- sc[sc$scenario == x$scenario, ]; g <- genes[[x$unit]]; m <- meta[[x$unit]]; o <- own[[x$unit]]
  base <- data.frame(unit = x$unit, block = x$block, scenario = x$scenario, kind = s$kind, config = x$config, filter = x$filter,
                     seeds = x$seeds, fit = x$fit, simulated = s$n, finished = done[[x$unit]],
                     fitted_here = if (is.null(o)) 0L else nrow(o), stringsAsFactors = FALSE)
  if (is.null(g)) return(base)
  t <- g[g$status == "tested", ]; n <- nrow(t)
  if (s$kind == "mixture") t <- t[0, ]  # rates for the mixture are in mixture.tsv
  k1 <- sum(t$p1 < ALPHA); k3 <- sum(t$p3 < ALPHA); n2 <- nrow(t)
  cbind(base, tested = n, excluded = sum(g$status == "excluded"), failed = sum(g$status == "failed"),
        excluded_reasons = counts(g$reason[g$status == "excluded"]), failed_reasons = counts(g$reason[g$status == "failed"]),
        positives_1 = k1, rate_1 = k1 / n2, lo_1 = cp(k1, n2)[1], hi_1 = cp(k1, n2)[2],
        positives_3 = k3, rate_3 = k3 / n2, lo_3 = cp(k3, n2)[1], hi_3 = cp(k3, n2)[2],
        rate_01_3 = mean(t$p3 < 0.01),
        class_1 = if (s$kind %in% nulllike) d1(k1, n2) else NA, class_3 = if (s$kind %in% nulllike) d1(k3, n2) else NA,
        reported_vs_lnl_gt_0.01 = sum(t$reported_minus_lnl > 0.01, na.rm = TRUE),
        reported_vs_lnl_gt_0.5 = sum(t$reported_minus_lnl > 0.5, na.rm = TRUE),
        max_reported_vs_lnl = suppressWarnings(max(t$reported_minus_lnl, na.rm = TRUE)),
        refits = sum(t$refit), negative_lrt = sum(t$lrt < 0, na.rm = TRUE),
        core_hours = if (is.null(m)) NA else m$seconds * m$threads / 3600, cpu = if (is.null(m)) NA else m$cpu,
        stringsAsFactors = FALSE)
})
units <- do.call(rbind, lapply(rows, function(r) { for (k in setdiff(names(rows[[which.max(sapply(rows, ncol))]]), names(r))) r[[k]] <- NA; r }))
wr(units, "units.tsv")
inc <- units[!units$finished, c("unit", "block", "scenario", "config")]
inc$reason <- vapply(inc$unit, function(u) { m <- meta[[u]]; if (is.null(m)) "not run" else sprintf("exit code %s", m$exit_code) }, "")
wr(inc, "incomplete.tsv")
cat("units:", nrow(units), " finished:", sum(units$finished), " incomplete:", nrow(inc), "\n")
cat("processors:", paste(unique(stats::na.omit(units$cpu)), collapse = " | "), "\n")
cat("core-hours spent:", round(sum(units$core_hours, na.rm = TRUE), 1), "\n")
cost <- stats::aggregate(core_hours ~ block, units[is.finite(units$core_hours), ], sum)
wr(rbind(cost, data.frame(block = "total", core_hours = sum(cost$core_hours))), "cost.tsv")

find <- function(s, c, f = "none", seeds = "1,2,3") { u <- un$unit[un$scenario == s & un$config == c & un$filter == f & un$seeds == seeds]; if (length(u)) u[1] else NA }

# --- D4: one start against three in the audited null-like cells.
aud <- un[un$block %in% c("B1", "B6") & un$scenario %in% c("b_N0_L25", "b_N0_L50_div2", "b_N1_L50", "b_M1_d02", "b_R1_At", "b_M4") &
          un$filter %in% c("none", "mask") & !(un$scenario == "b_M4" & un$filter != "none"), ]
audit <- do.call(rbind, lapply(aud$unit, function(u) {
  t <- tested(u); if (is.null(t)) return(data.frame(unit = u, genes = NA, positives_1 = NA, positives_3 = NA, changed = NA, class_1 = NA, class_3 = NA, pass = FALSE))
  ch <- sum((t$p1 < ALPHA) != (t$p3 < ALPHA)); c1 <- d1(sum(t$p1 < ALPHA), nrow(t)); c3 <- d1(sum(t$p3 < ALPHA), nrow(t))
  data.frame(unit = u, genes = nrow(t), positives_1 = sum(t$p1 < ALPHA), positives_3 = sum(t$p3 < ALPHA), changed = ch,
             class_1 = c1, class_3 = c3, pass = ch <= 2 && identical(c1, c3))
}))
decision <- if (nrow(audit) && all(audit$pass) && all(done[aud$unit])) "1start" else "3start"
wr(audit, "audit.tsv"); writeLines(decision, file.path(out, "audit_decision.txt"))
cat("D4 audit:", sum(audit$pass), "of", nrow(audit), "cells pass; B2 uses", decision, "\n")

# --- Noise floor: C1 with seeds 1-3 against 4-6 on the same genes.
callpair <- function(ua, ub, pa = "p3", pb = "p3") {
  A <- tested(ua); B <- tested(ub); if (is.null(A) || is.null(B)) return(NULL)
  m <- merge(A[, c("orthogroup", pa)], B[, c("orthogroup", pb)], by = "orthogroup", suffixes = c(".a", ".b"))
  list(a = m[[2]] < ALPHA, b = m[[3]] < ALPHA, n = nrow(m))
}
noise <- do.call(rbind, lapply(c("b_A10_p10", "b_A10_p20"), function(s) {
  p <- callpair(find(s, "C1"), find(s, "C1", seeds = "4,5,6")); if (is.null(p)) return(NULL)
  data.frame(scenario = s, genes = p$n, discordant = sum(p$a != p$b), rate_123 = mean(p$a), rate_456 = mean(p$b))
}))
if (!is.null(noise)) {
  dr <- sum(noise$discordant); nr <- sum(noise$genes); up <- cp(dr, nr)[2]
  noise <- rbind(noise, data.frame(scenario = "pooled", genes = nr, discordant = dr, rate_123 = NA, rate_456 = NA))
  noise$discordance_rate <- noise$discordant / noise$genes; noise$upper_95 <- c(rep(NA, nrow(noise) - 1), up)
  wr(noise, "noise.tsv")
} else up <- NA
threshold <- function(n) if (is.na(up)) NA else ceiling(n * up)

# --- Paired comparisons. Difference is b - a; for non-inferiority, b is the treatment.
wilson <- function(x, n, z = stats::qnorm(0.975)) { p <- x / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c0 - h, c0 + h) }
newcombe <- function(a, b) { # Newcombe (1998), method 10, paired proportions
  n <- length(a); n11 <- sum(a & b); n10 <- sum(a & !b); n01 <- sum(!a & b); n00 <- sum(!a & !b)
  p1 <- mean(a); p2 <- mean(b); w1 <- wilson(sum(a), n); w2 <- wilson(sum(b), n)
  den <- (n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00); phi <- if (den == 0) 0 else (n11 * n00 - n10 * n01) / sqrt(den)
  D <- p2 - p1
  c(D - sqrt(max(0, (p2 - w2[1])^2 - 2 * phi * (p2 - w2[1]) * (w1[2] - p1) + (w1[2] - p1)^2)),
    D + sqrt(max(0, (w2[2] - p2)^2 - 2 * phi * (w2[2] - p2) * (p1 - w1[1]) + (p1 - w1[1])^2)))
}
pr <- list()
pair <- function(label, question, ua, ub, pa = "p3", pb = "p3", ni = FALSE) {
  if (is.na(ua) || is.na(ub)) return(invisible())
  p <- callpair(ua, ub, pa, pb); if (is.null(p) || p$n == 0) return(invisible())
  lose <- sum(p$a & !p$b); gain <- sum(!p$a & p$b); ci <- newcombe(p$a, p$b); thr <- threshold(p$n)
  mc <- if (gain + lose) stats::binom.test(gain, gain + lose)$p.value else 1
  pr[[length(pr) + 1]] <<- data.frame(question = question, comparison = label, a = ua, b = ub, genes = p$n,
    rate_a = mean(p$a), rate_b = mean(p$b), a_only = lose, b_only = gain, difference = mean(p$b) - mean(p$a),
    diff_lo = ci[1], diff_hi = ci[2], mcnemar_p = mc, noise_threshold = thr,
    d3_real_difference = !is.na(thr) && mc < 0.05 && (gain + lose) > thr,
    noninferiority = if (!ni) NA else if (ci[1] > -MARGIN) "non-inferior" else if (ci[2] < -MARGIN) "inferior" else "inconclusive (precision)",
    stringsAsFactors = FALSE)
}
power <- c("A2_p20", "A3_p05", "A3_p10", "A3_p20", "A5_p10", "A5_p20", "A10_p05", "A10_p10", "A10_p20")
for (s in c(paste0("b_", power), "b_AC10_p05", "b_AC5_p10", "b_N0_L25", "b_N0_L50_div2", "b_N1_L50", "b_M1_d02", "b_M4"))
  pair("C3 vs C1", "configuration", find(s, "C1"), find(s, "C3"))
for (s in paste0("b_", c("A3_p20", "A10_p05", "A10_p10", "A10_p20"))) pair("C1 vs C0", "configuration", find(s, "C0"), find(s, "C1"))
for (s in c("b_A10_p10", "b_A10_p20", "b_M1_d03")) pair("C4 vs C1", "D5 triple hits", find(s, "C1"), find(s, "C4"))
for (s in c(paste0("b_", power), "b_AC10_p05", "b_AC5_p10", "b_A10_p10_indel")) {
  for (c in c("C1", "C3")) pair(paste("mask vs none,", c), "D6 masking, clean or clustered selection", find(s, c), find(s, c, "mask"), ni = TRUE)
  pair("mask_fg vs none, C1", "D6 masking, clean or clustered selection", find(s, "C1"), find(s, "C1", "mask_fg"), ni = TRUE)
  pair("mask_fg vs mask, C1", "D6 masking scope", find(s, "C1", "mask"), find(s, "C1", "mask_fg"))
}
for (c0 in c("A3_p20", "A5_p10", "A10_p10", "A10_p20")) {
  s <- paste0("b_", c0); r <- paste0(s, "_R1")
  for (c in c("C1", "C3")) {
    pair(paste("corrupted vs clean,", c), "D6 corruption effect", find(s, c), find(r, c))
    pair(paste("corrupted+mask vs clean,", c), "D6 masking, selection with corruption", find(s, c), find(r, c, "mask"), ni = TRUE)
  }
  pair("corrupted+mask_fg vs clean, C1", "D6 masking, selection with corruption", find(s, "C1"), find(r, "C1", "mask_fg"), ni = TRUE)
}
pair("gapflank vs none, C1", "D7 gap-flank mask, selection with indels", find("b_A10_p10_indel", "C1"), find("b_A10_p10_indel", "C1", "gapflank"), ni = TRUE)
for (f in c("true", "mask", "gapflank")) pair(paste(f, "vs realigned, C1"), "D7 indels, null", find("b_M4", "C1"), find("b_M4", "C1", f))
paired <- do.call(rbind, pr); wr(paired, "paired.tsv")

# --- D8: conservative rules, with model averaging (Selberg et al. 2024): p = w_B p_B + w_E p_E,
# Akaike weights from the alternative models' AICc. With 3 starts, each model's AICc is that of
# its best alternative fit, and its p is p3.
ma <- function(A, B, which) {
  m <- merge(A, B, by = "orthogroup", suffixes = c(".b", ".e"))
  ab <- if (which == 1) m$aicc_1.b else m$aicc_best.b; ae <- if (which == 1) m$aicc_1.e else m$aicc_best.e
  pb <- if (which == 1) m$p1.b else m$p3.b; pe <- if (which == 1) m$p1.e else m$p3.e
  mn <- pmin(ab, ae); wb <- exp(0.5 * (mn - ab)); we <- exp(0.5 * (mn - ae)); z <- wb + we
  data.frame(orthogroup = m$orthogroup, p_b = pb, p_e = pe, w_e = we / z, p_ma = (wb * pb + we * pe) / z)
}
rules <- list()
for (s in unique(un$scenario)) for (f in c("none", "mask")) {
  u1 <- find(s, "C1", f); u3 <- find(s, "C3", f); if (is.na(u1) || is.na(u3)) next
  A <- tested(u1); B <- tested(u3); if (is.null(A) || is.null(B)) next
  kind <- sc$kind[sc$scenario == s]
  for (w in c(1, 3)) {
    M <- ma(A, B, w); M <- merge(M, truth[truth$scenario == s, ], by = "orthogroup", all.x = TRUE)
    if (kind == "mixture") M <- M[!M$alternative, ]  # the mixture's null genes; its BH results are in mixture.tsv
    calls <- list(C1 = M$p_b < ALPHA, C3 = M$p_e < ALPHA, C1_and_C3 = M$p_b < ALPHA & M$p_e < ALPHA, MA = M$p_ma < ALPHA)
    for (r in names(calls)) { k <- sum(calls[[r]], na.rm = TRUE); n <- sum(!is.na(calls[[r]]))
      rules[[length(rules) + 1]] <- data.frame(scenario = s, kind = kind, filter = f, starts = w, rule = r, genes = n, positives = k,
        rate = k / n, lo = cp(k, n)[1], hi = cp(k, n)[2], class = if (kind %in% c(nulllike, "mixture")) d1(k, n) else NA,
        mean_w_busted_e = mean(M$w_e, na.rm = TRUE), share_w_busted_e_gt_0.5 = mean(M$w_e > 0.5, na.rm = TRUE)) }
  }
}
wr(do.call(rbind, rules), "rules.tsv")

# --- D9: the genome-like mixture after BH (empirical false discovery proportion of this one mixture).
mixrows <- list()
mix <- truth[truth$scenario == "b_MIX", ]
procs <- list()
for (f in c("none", "mask", "mask_fg")) for (c in c("C1", "C3")) {
  u <- find("b_MIX", c, f); t <- if (is.na(u)) NULL else tested(u)
  if (!is.null(t)) procs[[paste(c, f)]] <- stats::setNames(t$p3, t$orthogroup)
}
for (f in c("none", "mask")) {
  A <- tested(find("b_MIX", "C1", f)); B <- tested(find("b_MIX", "C3", f))
  if (!is.null(A) && !is.null(B)) { M <- ma(A, B, 3); procs[[paste("MA", f)]] <- stats::setNames(M$p_ma, M$orthogroup) }
}
for (q in c(0.05, 0.10)) {
  disc <- lapply(procs, function(p) names(p)[stats::p.adjust(p, "BH") < q])
  if (!is.null(disc[["C1 none"]]) && !is.null(disc[["C3 none"]])) disc[["C1 and C3 (each BH) none"]] <- intersect(disc[["C1 none"]], disc[["C3 none"]])
  if (!is.null(disc[["C1 mask"]]) && !is.null(disc[["C3 mask"]])) disc[["C1 and C3 (each BH) mask"]] <- intersect(disc[["C1 mask"]], disc[["C3 mask"]])
  for (r in names(disc)) {
    d <- mix[mix$orthogroup %in% disc[[r]], ]; fd <- sum(!d$alternative); nd <- nrow(d)
    tested_n <- if (!is.null(procs[[r]])) length(procs[[r]]) else NA
    mixrows[[length(mixrows) + 1]] <- data.frame(procedure = r, q = q, tested = tested_n, discoveries = nd, false = fd,
      fdp = if (nd) fd / nd else 0, fdp_lo = cp(fd, nd)[1], fdp_hi = cp(fd, nd)[2],
      power = sum(d$alternative) / sum(mix$alternative), false_corrupted = sum(!d$alternative & d$corrupted),
      false_corrupted_fg = sum(!d$alternative & d$corrupted_fg), stringsAsFactors = FALSE)
  }
}
if (length(mixrows)) wr(do.call(rbind, mixrows), "mixture.tsv")

# --- Filters: positions removed, genes changed.
fl <- do.call(rbind, lapply(sc$scenario, function(s) {
  d <- file.path(root, "sim", s); if (!file.exists(file.path(d, "DONE"))) return(NULL)
  do.call(rbind, lapply(c("mask", "mask_fg", "gapflank"), function(f) {
    lg <- file.path(d, paste0("filtered_", f), "filter_log.tsv"); if (!file.exists(lg)) return(NULL)
    lg <- utils::read.delim(lg, stringsAsFactors = FALSE)
    data.frame(scenario = s, filter = f, genes = nrow(lg), genes_changed = sum(lg$codons_masked > 0),
               codons_masked_mean = mean(lg$codons_masked), fraction_removed = sum(lg$codons_masked) / sum(lg$codons_present))
  }))
}))
wr(fl, "filtering.tsv")

# --- Starts: what the second and third starts changed.
st <- do.call(rbind, lapply(un$unit[grepl(",", un$seeds)], function(u) {
  t <- own[[u]]; if (is.null(t)) return(NULL); t <- t[t$status == "tested", ]; if (!nrow(t)) return(NULL)
  data.frame(unit = u, genes_fitted = nrow(t), changed_1_vs_3 = sum((t$p1 < ALPHA) != (t$p3 < ALPHA)),
             gained = sum(t$p1 >= ALPHA & t$p3 < ALPHA), lost = sum(t$p1 < ALPHA & t$p3 >= ALPHA),
             alt_improved_0.5 = mean(t$lnl_unconstrained_best - t$lnl_unconstrained > 0.5, na.rm = TRUE),
             null_improved_0.5 = mean(t$lnl_null_best - t$lnl_constrained > 0.5, na.rm = TRUE))
}))
if (!is.null(st)) wr(st, "starts.tsv")
cat("summary written to", out, "\n")
