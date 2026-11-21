# Summarize a RELAX parametric-bootstrap run (phylo-relax-boot.sb): per gene, the observed LRT,
# how closely the refitted null matches RELAX's own null fit, the bootstrap LRTs, and the
# Monte Carlo p-value (k + 1) / (B + 1). Pooled: bootstrap LRTs against chi-square(1).
#   Rscript --vanilla phylo-relax-boot-summary.R <run dir>
d <- commandArgs(TRUE)[1]
x <- utils::read.delim(file.path(d, "boot.tsv"), stringsAsFactors = FALSE)
tm <- utils::read.delim(file.path(d, "timing.tsv"), header = FALSE, col.names = c("gene", "obs_sim_s", "reps_wall_s", "seed"))
genes <- unique(x$orthogroup); rows <- list(); all <- numeric()
for (g in genes) {
  y <- x[x$orthogroup == g, ]; obs <- y[y$rep == "observed", ]; reps <- y[!y$rep %in% c("observed", "refit_null"), ]
  lrt <- pmax(as.numeric(reps$lrt), 0); lrt <- lrt[!is.na(lrt)]; o <- max(as.numeric(obs$lrt), 0); all <- c(all, lrt)
  rows[[g]] <- data.frame(gene = g, lrt_obs = o, p_chisq = stats::pchisq(o, 1, lower.tail = FALSE), B = length(lrt), failed = sum(is.na(as.numeric(reps$lrt))),
    p_boot = (sum(lrt >= o) + 1) / (length(lrt) + 1), mean_boot = mean(lrt),
    refit_minus_relax_null = as.numeric(y$lrt[y$rep == "refit_null"]) - as.numeric(obs$lnl_null))
}
r <- do.call(rbind, rows); print(r, digits = 4, row.names = FALSE)
p <- c(.5, .75, .9, .95); n <- length(all)
cat(sprintf("\nBootstrap LRTs pooled over %d genes: %d replicates; mean %.3f (SE %.3f; chi-square(1): 1; the simulated null genes of the first study: 1.27)\n", length(genes), n, mean(all), stats::sd(all) / sqrt(n)))
print(rbind(bootstrap = stats::quantile(all, p), chisq1 = stats::qchisq(p, 1)), digits = 3)
cat(sprintf("Share of bootstrap LRTs above 3.84: %.3f (chi-square(1): 0.05)\n", mean(all > stats::qchisq(.95, 1))))
cat(sprintf("Refitted null minus RELAX's null log-likelihood: median %.2g, range %.2g to %.2g\n", stats::median(r$refit_minus_relax_null), min(r$refit_minus_relax_null), max(r$refit_minus_relax_null)))
cat(sprintf("Cost: observed fit + null refit + simulation %.0f s per gene; %d-replicate wall time %.0f s on 16 cores\n", mean(tm$obs_sim_s), max(r$B), mean(tm$reps_wall_s)))
