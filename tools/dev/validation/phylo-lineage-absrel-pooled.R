# Exploratory E3 (docs/dev/validation-audit.md). Usage: Rscript --vanilla phylo-lineage-absrel-pooled.R <first study root>
root <- commandArgs(TRUE)[1]
set.seed(1)
for (u in c("B0_T5", "B0_T8b", "B0_T8p", "B0_T16")) {
  b <- read.delim(file.path(root, "units", u, "absrel", "absrel_branches.tsv"))
  b <- b[b$in_species_tree & !grepl("+", b$species_tree_branches, fixed = TRUE) & is.finite(b$adaptive_dn + b$adaptive_ds), ]
  per <- do.call(rbind, lapply(split(b, b$species_tree_branches), function(x) {
    est <- sum(x$adaptive_dn) / sum(x$adaptive_ds)
    bo <- replicate(500, { i <- sample.int(nrow(x), replace = TRUE); sum(x$adaptive_dn[i]) / sum(x$adaptive_ds[i]) })
    ci <- quantile(bo, c(.025, .975))
    data.frame(branch = x$species_tree_branches[1], n = nrow(x), est = est, lo = ci[1], hi = ci[2], len = median(x$adaptive_length), base = median(x$baseline_omega), multi = mean(x$rate_classes > 1))
  }))
  cat(sprintf("%s: %d branches; median pooled %.3f (truth 0.22, rel. bias %+.1f%%), range %.3f-%.3f; intervals covering 0.22: %d of %d; median per-gene baseline omega %.3f; genes with >1 class %.0f%%\n",
    u, nrow(per), median(per$est), 100 * (median(per$est) / 0.22 - 1), min(per$est), max(per$est), sum(per$lo <= 0.22 & per$hi >= 0.22), nrow(per), median(per$base), 100 * mean(per$multi)))
  print(round(cor(per$len, per$est / 0.22, method = "spearman"), 2))
}
