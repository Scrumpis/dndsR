# Exploratory (docs/dev/validation-audit.md, E2): RELAX's LRT on null genes, against chi-square(1).
#   Rscript --vanilla phylo-lineage-relax-diag.R <relax_results.tsv> [<relax_results.tsv> ...]
f <- commandArgs(TRUE); names(f) <- basename(dirname(dirname(dirname(f))))
x <- do.call(rbind, lapply(names(f), function(n) { r <- read.delim(f[n]); r <- r[r$status == "tested", ]; r$set <- n; r[, c("set", "orthogroup", "lrt", "p_value", "k", "k_at_bound", grep("^(ref|test)_(omega|prop)_", names(r), value = TRUE), "mg94_omega_test", "mg94_omega_reference")] }))
x$lrt <- pmax(x$lrt, 0); x$sig <- x$p_value <= 0.05
om <- as.matrix(x[, grep("^ref_omega_", names(x))]); pr <- as.matrix(x[, grep("^ref_prop_", names(x))])
x$max_omega <- apply(om, 1, max, na.rm = TRUE); x$prop_max <- pr[cbind(seq_len(nrow(x)), apply(om, 1, which.max))]
x$phantom <- x$max_omega > 2 & x$prop_max < 0.05   # a small class far above the simulated maximum (1.5 for C0)
cat("genes:", nrow(x), "; with a small class (prop < 0.05) at omega > 2:", sum(x$phantom), "\n")
for (g in c(TRUE, FALSE)) { y <- x[x$phantom == g, ]; cat(sprintf("phantom=%s: n %d, size %.1f%%, mean LRT %.2f, median %.2f\n", g, nrow(y), 100 * mean(y$sig), mean(y$lrt), median(y$lrt))) }
cat("Spearman LRT ~ log max omega:", round(cor(x$lrt, log(x$max_omega), method = "spearman"), 3), "\n")
cat("k quantiles among significant:\n"); print(round(quantile(x$k[x$sig], c(0, .1, .25, .5, .75, .9, 1)), 3))
cat("k at bound among significant:", sum(x$k_at_bound[x$sig]), "of", sum(x$sig), "\n")
cat("LRT mixture check: share exactly 0:", mean(x$lrt < 1e-6), "; quantiles vs chi2(1):\n")
p <- c(.5, .75, .9, .95, .99); print(rbind(observed = quantile(x$lrt, p), chisq1 = qchisq(p, 1), ratio = quantile(x$lrt, p) / qchisq(p, 1)))
cat("Fit of scaled chi2: a*chi2(1), a = mean:", round(mean(x$lrt), 3), "; a from 95th:", round(quantile(x$lrt, .95) / qchisq(.95, 1), 3), "\n")
cat("By set (size, mean):\n"); print(aggregate(cbind(sig, lrt) ~ set, x, mean))
cat("mg94 test vs ref omega, sig vs not (median ratio):", round(tapply(x$mg94_omega_test / x$mg94_omega_reference, x$sig, median), 3), "\n")
