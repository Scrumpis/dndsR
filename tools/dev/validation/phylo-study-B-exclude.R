# Study B, the "warn" and "exclude" classifications (docs/dev/lineage-robustness-study.md,
# "Classification of each error, per method", fixed before running), on the mixture unit:
#  - detection: sensitivity of each signal for the genes the corruption log marks as affected,
#    and its false-flag rate on unaffected genes;
#  - three arms per test: as given; masked (phylo-windows mask); as given with genes flagged by
#    the window profile excluded;
#  - false discovery proportion after the documented correction, power on the selected genes,
#    and each arm's power difference from the unmasked arm, with a 95% interval from 2,000
#    bootstrap resamples of genes; power counts every selected test of the unmasked arm, so an
#    excluded gene counts as not found (the non-inferiority rule: loss under 5 points, lower bound above -10).
#   Rscript --vanilla phylo-study-B-exclude.R <results root> <summary prefix of phylo-study-B-summary.R> <out prefix>
a <- commandArgs(TRUE); root <- a[1]; sprefix <- a[2]; prefix <- a[3]
u <- "B-MIX-T8b"; ud <- file.path(root, "units", u)
genes <- utils::read.delim(file.path(ud, "sim", "genes.tsv"), stringsAsFactors = FALSE)
log <- utils::read.delim(pipe(sprintf("tar -xzOf %s corrupted/corrupt_log.tsv", file.path(ud, "sim", "corrupted.tar.gz"))), stringsAsFactors = FALSE)
affected <- unique(sub("\\.fasta$", "", log$alignment))
dt <- utils::read.delim(paste0(sprefix, "_detection_genes.tsv"), stringsAsFactors = FALSE); dt <- dt[dt$unit == u, ]
sig <- c("dissimilar_window", "high_ds", "few_synonymous", "tip_ds_3x_median", "qc_flag")
det <- data.frame(signal = sig, affected = length(affected),
  sensitivity = vapply(sig, function(s) mean(dt[[s]][dt$gene %in% affected]), 0),
  false_flag_unaffected = vapply(sig, function(s) mean(dt[[s]][!dt$gene %in% affected]), 0))
utils::write.table(det, paste0(prefix, "_mix_detection.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
flagged <- dt$gene[dt$dissimilar_window]
sel <- genes$gene[genes$selected]
# Truth per branch comes from the simulator (sim/branch_truth.tsv), not from the branches table's
# `foreground` column, which is TRUE for every tested branch when aBSREL tests all branches.
bt <- utils::read.delim(file.path(ud, "sim", "branch_truth.tsv"), stringsAsFactors = FALSE)
fg_ids <- bt$branch_id[bt$foreground]

rd <- function(run, what) utils::read.delim(file.path(ud, run, sprintf("%s_%s.tsv", run, what)), stringsAsFactors = FALSE)
get <- function(run, test) {
  r <- rd(run, "results"); t <- r[r$status == "tested", ]
  if (test == "p_two_vs_m0") return(data.frame(og = t$orthogroup, p = t$p_two_vs_m0, true = t$orthogroup %in% sel))
  if (test == "absrel p_value per gene") return(data.frame(og = t$orthogroup, p = t$p_value, true = t$orthogroup %in% sel))
  b <- rd(run, "branches"); b <- b[b$orthogroup %in% t$orthogroup & !is.na(b$p_raw), ]
  data.frame(og = b$orthogroup, p = b$p_raw, true = b$orthogroup %in% sel & b$species_tree_branches %in% fg_ids, branch = b$species_tree_branches)
}
# The documented corrections: BH across genes; for aBSREL branches, BH per species-tree branch (p_bh_lineage).
adjust <- function(x) if (is.null(x$branch)) stats::p.adjust(x$p, "BH") else ave(x$p, x$branch, FUN = function(v) stats::p.adjust(v, "BH"))
# Power over every true test of the unmasked arm: an excluded gene counts as not found.
score <- function(x, ntrue = NULL) { x <- x[!is.na(x$p), ]; d <- adjust(x) <= 0.05
  if (is.null(ntrue)) ntrue <- sum(x$true); c(disc = sum(d), false = sum(d & !x$true), power = sum(d[x$true]) / ntrue, ntrue = ntrue) }
arms <- list(
  absrel = list(tests = c("absrel p_raw per branch (p_bh_lineage)", "absrel p_value per gene"), runs = c("absrel", "absrel_mask")),
  branch_omega = list(tests = "p_two_vs_m0", runs = c("branch_omega", "branch_omega_mask")))
set.seed(20261122)
out <- list()
for (m in names(arms)) for (tn in arms[[m]]$tests) {
  test <- sub(" \\(p_bh_lineage\\)$", "", tn)
  x0 <- get(arms[[m]]$runs[1], test); xm <- get(arms[[m]]$runs[2], test)
  xe <- x0[!x0$og %in% flagged, ]
  sets <- list(as_given = x0, masked = xm, excluded_flagged = xe)
  ogs <- unique(x0$og)
  boot <- replicate(2000, {
    s <- sample(ogs, replace = TRUE); w <- table(s)
    rs <- lapply(sets, function(x) { i <- unlist(lapply(names(w), function(g) rep(which(x$og == g), w[[g]]))); x <- x[i, ]
      x }); nt <- sum(rs$as_given$true[!is.na(rs$as_given$p)])
    vapply(rs, function(x) score(x, nt)[["power"]], 0)
  })
  for (arm in names(sets)) {
    s <- score(sets[[arm]], sum(x0$true[!is.na(x0$p)])); ci <- if (s[["disc"]]) stats::binom.test(s[["false"]], s[["disc"]])$conf.int else c(NA, NA)
    dpow <- boot[arm, ] - boot["as_given", ]
    out[[length(out) + 1]] <- data.frame(method = m, test = tn, arm = arm, genes = length(unique(sets[[arm]]$og)),
      discoveries = s[["disc"]], false = s[["false"]], fdp = if (s[["disc"]]) s[["false"]] / s[["disc"]] else NA, fdp_lo = ci[1], fdp_hi = ci[2],
      power = s[["power"]], true_tests = s[["ntrue"]], power_diff = s[["power"]] - score(x0)[["power"]],
      power_diff_lo = stats::quantile(dpow, 0.025, na.rm = TRUE), power_diff_hi = stats::quantile(dpow, 0.975, na.rm = TRUE))
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL
utils::write.table(res, paste0(prefix, "_mix_arms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("affected genes %d of %d; selected %d; flagged by window profile %d\n", length(affected), nrow(genes), length(sel), length(flagged)))
print(det, row.names = FALSE); print(res, row.names = FALSE, digits = 3)
