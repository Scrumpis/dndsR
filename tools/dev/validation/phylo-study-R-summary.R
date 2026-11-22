# Study R, stage 1 summary (docs/dev/relax-bootstrap-study.md, "Stage 1"): RD2-RD7.
#   Rscript --vanilla phylo-study-R-summary.R <results root> <out prefix>
# Reads <root>/units/<unit>/relax/relax_results.tsv and <root>/calib/<unit>/boot/boot.tsv.
# HyPhy's p_value is read unchanged; the two new p-values are computed here, beside it:
#   p_scaled  P(chi-square(1) > LRT / c), c = mean of the dataset's calibration-set bootstrap LRTs
#   p_pool    (1 + #{bootstrap LRT >= LRT}) / (n + 1); beyond the largest bootstrap LRT, p_scaled x n / (n + 1)
a <- commandArgs(TRUE); root <- a[1]; prefix <- a[2]
set.seed(20261124)
cp <- function(x, n) if (n == 0) c(NA, NA) else stats::binom.test(x, n)$conf.int[1:2]
bradley <- function(x, n) { if (n == 0) return("no data"); ci <- cp(x, n); if (ci[1] >= 0.025 && ci[2] <= 0.075) "within tolerance" else if (ci[1] > 0.05) "anticonservative" else if (ci[2] <= 0.05) "calibrated" else "outside tolerance or inconclusive" }
res <- function(u) {
  f <- file.path(root, "units", u, "relax", "relax_results.tsv")
  if (!file.exists(f)) return(NULL)
  r <- utils::read.delim(f, stringsAsFactors = FALSE); r[r$status == "tested", ]
}
boot <- function(u) {
  f <- file.path(root, "calib", u, "boot", "boot.tsv")
  if (!file.exists(f)) return(NULL)
  b <- utils::read.delim(f, stringsAsFactors = FALSE)
  b <- b[grepl("^[0-9]+$", b$rep), ]; b$lrt <- pmax(0, as.numeric(b$lrt)); b[is.finite(b$lrt), ]
}
pvals <- function(lrt, bl) {
  c0 <- mean(bl); n <- length(bl); s <- sort(bl)
  ps <- stats::pchisq(lrt / c0, 1, lower.tail = FALSE)
  ge <- n - findInterval(lrt, s, left.open = TRUE)
  pp <- ifelse(lrt > max(s), ps * n / (n + 1), (1 + ge) / (n + 1))
  list(scale = c0, p_scaled = ps, p_pool = pp)
}
rows <- list(); het <- list()
cells <- c("R-N16", "R-N8", "R-N5", "R-K05", "R-K2")
dat <- list()
for (u in cells) {
  r <- res(u); b <- boot(u)
  if (is.null(r)) next
  dat[[u]] <- list(r = r, b = b)
  p <- if (!is.null(b)) pvals(r$lrt, b$lrt) else list(scale = NA, p_scaled = rep(NA, nrow(r)), p_pool = rep(NA, nrow(r)))
  dat[[u]]$p <- p
  for (k in c("p_value", "p_scaled", "p_pool")) {
    pv <- if (k == "p_value") r$p_value else p[[k]]
    x <- sum(pv <= 0.05, na.rm = TRUE); n <- sum(!is.na(pv)); ci <- cp(x, n)
    right <- if (u == "R-K05") r$k < 1 else if (u == "R-K2") r$k > 1 else rep(NA, nrow(r))
    rows[[length(rows) + 1L]] <- data.frame(unit = u, p = k, scale = p$scale, genes = n, positives = x, rate = x / n, lo = ci[1], hi = ci[2],
      verdict = if (grepl("^R-N", u)) bradley(x, n) else NA, right_direction = if (grepl("^R-K", u)) mean(right[pv <= 0.05], na.rm = TRUE) else NA)
  }
  if (!is.null(b)) { # RD3: per-gene mean bootstrap LRT, against replicate noise (permutation of replicates across genes)
    m <- tapply(b$lrt, b$orthogroup, mean); v <- stats::var(m)
    perm <- replicate(10000, stats::var(tapply(sample(b$lrt), b$orthogroup, mean)))
    het[[u]] <- data.frame(unit = u, genes = length(m), boot_fits = nrow(b), mean_lrt = mean(b$lrt), se = stats::sd(m) / sqrt(length(m)),
      q95 = stats::quantile(b$lrt, 0.95, names = FALSE), share_gt_3.84 = mean(b$lrt > 3.84),
      var_gene_means = v, var_expected = mean(perm), p_heterogeneity = mean(perm >= v))
  }
}
size <- do.call(rbind, rows)
# RD5: size-adjusted chi-square power, at the critical LRT giving 5% in R-N16
if (!is.null(dat[["R-N16"]])) {
  crit <- stats::quantile(dat[["R-N16"]]$r$lrt, 0.95, names = FALSE)
  for (u in c("R-K05", "R-K2")) if (!is.null(dat[[u]])) {
    x <- sum(dat[[u]]$r$lrt > crit); n <- nrow(dat[[u]]$r); ci <- cp(x, n)
    size <- rbind(size, data.frame(unit = u, p = "chi-square, size-adjusted", scale = NA, genes = n, positives = x, rate = x / n, lo = ci[1], hi = ci[2], verdict = NA,
      right_direction = mean((if (u == "R-K05") dat[[u]]$r$k < 1 else dat[[u]]$r$k > 1)[dat[[u]]$r$lrt > crit])))
  }
}
utils::write.table(size, paste0(prefix, "_size_power.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(het)) utils::write.table(do.call(rbind, het), paste0(prefix, "_heterogeneity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# RD2: alignment length
len <- list()
for (u in c("R-N16", "RD-L1000", "RD-L3000")) { r <- res(u); if (is.null(r)) next
  x <- sum(r$p_value <= 0.05); ci <- cp(x, nrow(r))
  len[[u]] <- data.frame(unit = u, codons = r$n_codons[1], genes = nrow(r), mean_lrt = mean(pmax(0, r$lrt)), se = stats::sd(pmax(0, r$lrt)) / sqrt(nrow(r)),
    rate = x / nrow(r), lo = ci[1], hi = ci[2]) }
if (length(len)) utils::write.table(do.call(rbind, len), paste0(prefix, "_length.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# RD6: genome-like mixtures (all R-N16 genes, 30 each of R-K05 and R-K2), one calibration set (R-N16's), BH at 0.05
if (all(c("R-N16", "R-K05", "R-K2") %in% names(dat)) && !is.null(dat[["R-N16"]]$b)) {
  bl <- dat[["R-N16"]]$b$lrt
  mix <- t(replicate(2000, {
    i0 <- seq_len(nrow(dat[["R-N16"]]$r)); i1 <- sample.int(nrow(dat[["R-K05"]]$r), 30); i2 <- sample.int(nrow(dat[["R-K2"]]$r), 30)
    lrt <- c(dat[["R-N16"]]$r$lrt[i0], dat[["R-K05"]]$r$lrt[i1], dat[["R-K2"]]$r$lrt[i2])
    hy <- c(dat[["R-N16"]]$r$p_value[i0], dat[["R-K05"]]$r$p_value[i1], dat[["R-K2"]]$r$p_value[i2])
    alt <- rep(c(FALSE, TRUE), c(length(i0), 60)); p <- pvals(lrt, bl)
    unlist(lapply(list(p_value = hy, p_scaled = p$p_scaled, p_pool = p$p_pool), function(pv) {
      d <- stats::p.adjust(pv, "BH") <= 0.05
      c(fdp = if (any(d)) sum(d & !alt) / sum(d) else 0, power = mean(d[alt]))
    }))
  }))
  m6 <- data.frame(statistic = colnames(mix), mean = colMeans(mix), lo = apply(mix, 2, stats::quantile, 0.025), hi = apply(mix, 2, stats::quantile, 0.975))
  utils::write.table(m6, paste0(prefix, "_mixture.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
# RD1: HyPhy versions, run directly on R-N16's dndsR inputs (phylo-lineage-relax-direct.sb, MODE=dndsr, seed 1)
if (!is.null(dat[["R-N16"]])) {
  ver <- list()
  for (f in Sys.glob(file.path(root, "direct", "*", "direct_dndsr_seed1.tsv"))) {
    v <- basename(dirname(f)); d <- utils::read.delim(f, stringsAsFactors = FALSE)
    m <- merge(dat[["R-N16"]]$r[, c("orthogroup", "lrt")], d[, c("orthogroup", "lrt")], by = "orthogroup", suffixes = c("_dndsr", "_direct"))
    m$lrt_direct <- suppressWarnings(as.numeric(m$lrt_direct)); ok <- is.finite(m$lrt_direct)
    x <- sum(m$lrt_direct[ok] > 3.84); ci <- cp(x, sum(ok))
    ver[[v]] <- data.frame(version = v, genes = nrow(m), failed = sum(!ok), share_gt_3.84 = x / sum(ok), lo = ci[1], hi = ci[2],
      dndsr_share_gt_3.84 = mean(m$lrt_dndsr[ok] > 3.84), mean_lrt = mean(pmax(0, m$lrt_direct[ok])),
      agree_0.01 = mean(abs(m$lrt_direct[ok] - m$lrt_dndsr[ok]) < 0.01), max_abs_diff = max(abs(m$lrt_direct[ok] - m$lrt_dndsr[ok])))
  }
  if (length(ver)) utils::write.table(do.call(rbind, ver), paste0(prefix, "_versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
for (f in Sys.glob(paste0(prefix, "_*.tsv"))) { cat("\n==", basename(f), "\n"); print(utils::read.delim(f), row.names = FALSE) }
