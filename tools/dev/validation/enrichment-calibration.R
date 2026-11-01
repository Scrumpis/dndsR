# Null calibration of Fisher term enrichment, as term_enrichment() does it.
#   Rscript enrichment-calibration.R <dnds_annot.tsv> <real enrichment tsv> <term column> <nperm> <seed> <outdir>
#
# Fast reimplementation (hypergeometric tail on a sparse matrix) so thousands of nulls are
# practical. It must first reproduce the real output exactly, or the script stops.
#
# Three nulls for "which genes are positive", all keeping the number of positives:
#   permute:  unrestricted. Assumes genes are exchangeable.
#   strata:   within length x dS x annotation-richness strata. Keeps those confounders.
#   rotate:   circular shift along the genome. Keeps clustering, such as tandem arrays.
suppressMessages(library(Matrix))
a <- commandArgs(TRUE)
annot <- a[1]; real_file <- a[2]; col <- a[3]; nperm <- as.integer(a[4]); seed <- as.integer(a[5]); out <- a[6]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
min_total <- 2L; min_pos <- 2L; alpha <- 0.05

d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d$dNdS <- suppressWarnings(as.numeric(d$dNdS))
d <- d[is.finite(d$dNdS) & d$dNdS < 10, ]
d <- d[order(d$q_gff_seqname, d$q_gff_start, d$query_id), ]
n <- nrow(d)
terms_of <- strsplit(ifelse(is.na(d[[col]]), "", d[[col]]), ";", fixed = TRUE)
terms_of <- lapply(terms_of, function(x) unique(trimws(x[nzchar(trimws(x))])))
n_terms <- lengths(terms_of)
annotated <- n_terms > 0
ids <- sort(unique(unlist(terms_of)))
m <- sparseMatrix(i = rep(seq_len(n), n_terms), j = match(unlist(terms_of), ids), x = 1,
                  dims = c(n, length(ids)), dimnames = list(NULL, ids))
K <- colSums(m); N <- sum(annotated)

enrich <- function(pos, filter = TRUE) {
  counts <- as.numeric(crossprod(m, as.numeric(pos)))
  npos <- sum(pos & annotated)
  p <- phyper(counts - 1, K, N - K, npos, lower.tail = FALSE)
  keep <- K >= min_total & (if (filter) counts >= min_pos else TRUE)
  adj <- rep(NA_real_, length(p)); adj[keep] <- p.adjust(p[keep], "BH")
  list(a = counts, p = p, keep = keep, adj = adj, npos = npos)
}

# The check: identical to the package at the observed labels.
pos0 <- d$dNdS > 1
obs <- enrich(pos0)
real <- read.table(real_file, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
idx <- match(real$term, ids)
stopifnot(!anyNA(idx), length(idx) == sum(obs$keep), all(obs$keep[idx]),
          all.equal(obs$a[idx], as.numeric(real$pos_count)),
          max(abs(obs$p[idx] - real$p_value)) < 1e-10, max(abs(obs$adj[idx] - real$p_adj)) < 1e-10)
cat(sprintf("Reimplementation equals term_enrichment on %d terms: max |dp| %.2g, max |dpadj| %.2g\n",
            length(idx), max(abs(obs$p[idx] - real$p_value)), max(abs(obs$adj[idx] - real$p_adj))))

ntile <- function(x, k) cut(rank(x, ties.method = "first"), k, labels = FALSE)
rich <- ifelse(n_terms == 0, 0L, ifelse(n_terms <= median(n_terms[annotated]), 1L, 2L))
strata <- interaction(ntile(d$q_len, 5), ntile(d$dS, 5), rich, drop = TRUE)
by_stratum <- split(seq_len(n), strata)
draw <- list(
  permute = function() sample(pos0),
  strata = function() {
    p <- pos0
    for (s in by_stratum) p[s] <- pos0[s[sample.int(length(s))]]
    p
  },
  rotate = function() {
    k <- sample(100:(n - 100), 1)
    pos0[c((k + 1):n, 1:k)]
  }
)

# Filters on term size alone do not depend on the outcome, so BH stays valid after them.
size_filters <- c(5L, 10L)
set.seed(seed)
summ <- list(); tail_count <- list(); min_p <- list()
for (nm in names(draw)) {
  r <- matrix(NA_real_, nperm, 7, dimnames = list(NULL, c("n_tested", "n_sig", "n_sig_nofilter", "frac_p05_tested",
                                                          "n_all", "frac_p05_all", "npos")))
  tc <- numeric(length(ids)); mp <- numeric(nperm)
  size_fp <- matrix(0, nperm, length(size_filters), dimnames = list(NULL, paste0("any_sig_size", size_filters)))
  for (i in seq_len(nperm)) {
    e <- enrich(draw[[nm]]())
    u <- p.adjust(e$p[K >= min_total], "BH")
    r[i, ] <- c(sum(e$keep), sum(e$adj[e$keep] <= alpha), sum(u <= alpha), mean(e$p[e$keep] < 0.05),
                sum(K >= min_total), mean(e$p[K >= min_total] < 0.05), e$npos)
    for (z in seq_along(size_filters)) size_fp[i, z] <- sum(p.adjust(e$p[K >= size_filters[z]], "BH") <= alpha)
    tc <- tc + (e$a >= obs$a)
    mp[i] <- min(e$p[K >= min_total])
  }
  summ[[nm]] <- data.frame(null = nm, nperm = nperm, mean_tested = mean(r[, 1]),
    any_sig = mean(r[, 2] > 0), mean_sig = mean(r[, 2]), p99_sig = quantile(r[, 2], 0.99), max_sig = max(r[, 2]),
    any_sig_nofilter = mean(r[, 3] > 0), mean_sig_nofilter = mean(r[, 3]),
    frac_p05_tested = mean(r[, 4]), frac_p05_all = mean(r[, 6]), mean_npos = mean(r[, 7]),
    t(colMeans(size_fp > 0)), row.names = NULL)
  tail_count[[nm]] <- tc; min_p[[nm]] <- mp
}
calib <- do.call(rbind, summ)
write.table(calib, file.path(out, "null_calibration.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(calib)

# Observed terms under each null. A gene family in a tandem array counts once per cluster.
ord_of <- seq_len(n)
clusters <- function(idx, gap = 5L) if (!length(idx)) 0L else 1L + sum(diff(sort(idx)) > gap)
sig <- which(obs$keep & obs$adj <= alpha)
collapsed_p <- rep(NA_real_, length(ids))
tested <- which(obs$keep)
for (j in tested) {
  g <- which(m[, j] > 0)
  gp <- g[pos0[g]]
  Kc <- clusters(g); ac <- clusters(gp)
  collapsed_p[j] <- phyper(ac - 1, Kc, (N - K[j]) , obs$npos - (obs$a[j] - ac), lower.tail = FALSE)
}
collapsed_adj <- rep(NA_real_, length(ids)); collapsed_adj[tested] <- p.adjust(collapsed_p[tested], "BH")
# The candidate correction: BH over every term of adequate size, whatever its positive count.
adj_all <- function(min_size) { x <- rep(NA_real_, length(ids)); k <- K >= min_size; x[k] <- p.adjust(obs$p[k], "BH"); x }
fixed2 <- adj_all(2L); fixed5 <- adj_all(5L); fixed10 <- adj_all(10L)
res <- data.frame(term = ids[sig], pos = obs$a[sig], total = K[sig], p = obs$p[sig], p_adj = obs$adj[sig],
  p_adj_all = fixed2[sig], p_adj_size5 = fixed5[sig], p_adj_size10 = fixed10[sig],
  pos_clusters = vapply(sig, function(j) { g <- which(m[, j] > 0); clusters(g[pos0[g]]) }, 1L),
  total_clusters = vapply(sig, function(j) clusters(which(m[, j] > 0)), 1L),
  collapsed_adj = collapsed_adj[sig],
  stringsAsFactors = FALSE)
for (nm in names(draw)) {
  res[[paste0("emp_p_", nm)]] <- (1 + tail_count[[nm]][sig]) / (nperm + 1)
  res[[paste0("wy_", nm)]] <- vapply(obs$p[sig], function(x) (1 + sum(min_p[[nm]] <= x)) / (nperm + 1), 0)
}
res <- res[order(res$p), ]
write.table(res, file.path(out, "observed_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Observed: %d tested, %d significant at BH %.2f; survive: collapsed %d; Westfall-Young permute %d, strata %d, rotate %d\n",
            sum(obs$keep), nrow(res), alpha, sum(res$collapsed_adj <= alpha), sum(res$wy_permute <= alpha),
            sum(res$wy_strata <= alpha), sum(res$wy_rotate <= alpha)))

cat(sprintf("Terms by size: at least 2, %d; 5, %d; 10, %d. Significant with BH before filtering: size 2, %d; size 5, %d; size 10, %d\n",
            sum(K >= 2), sum(K >= 5), sum(K >= 10), sum(fixed2 <= alpha, na.rm = TRUE), sum(fixed5 <= alpha, na.rm = TRUE), sum(fixed10 <= alpha, na.rm = TRUE)))

# Confounders of being positive.
q5 <- function(x) ntile(x, 5)
conf <- rbind(
  data.frame(variable = "length", quintile = 1:5, rate = as.numeric(tapply(pos0, q5(d$q_len), mean)),
             median = as.numeric(tapply(d$q_len, q5(d$q_len), median))),
  data.frame(variable = "dS", quintile = 1:5, rate = as.numeric(tapply(pos0, q5(d$dS), mean)),
             median = as.numeric(tapply(d$dS, q5(d$dS), median))))
write.table(conf, file.path(out, "confounder_rates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(conf)
fit <- glm(pos0 ~ log10(d$q_len) + log10(d$dS) + annotated + log1p(n_terms), family = binomial)
print(round(summary(fit)$coefficients, 3))
cat(sprintf("Positive rate: annotated %.4f (%d of %d), unannotated %.4f (%d of %d)\n", mean(pos0[annotated]), sum(pos0[annotated]), N,
            mean(pos0[!annotated]), sum(pos0[!annotated]), sum(!annotated)))
cat(sprintf("Median length aa: positive %g, other %g. Median dS: positive %.3f, other %.3f\n",
            median(d$q_len[pos0]), median(d$q_len[!pos0]), median(d$dS[pos0]), median(d$dS[!pos0])))

# Candidate selection: how much do the conclusions depend on the rule?
base <- ids[sig]
rule <- function(name, p, keep = rep(TRUE, n)) {
  mk <- m[keep, , drop = FALSE]; p <- p[keep]; ann <- annotated[keep]
  cnt <- as.numeric(crossprod(mk, as.numeric(p))); Kk <- colSums(mk); Nk <- sum(ann); np <- sum(p & ann)
  pv <- phyper(cnt - 1, Kk, Nk - Kk, np, lower.tail = FALSE)
  ok <- Kk >= min_total & cnt >= min_pos
  s <- ids[ok][p.adjust(pv[ok], "BH") <= alpha]
  data.frame(rule = name, genes = sum(keep), positive = sum(p), annotated_positive = np, significant = length(s),
             shared_with_default = length(intersect(s, base)),
             jaccard = if (length(union(s, base))) length(intersect(s, base)) / length(union(s, base)) else NA)
}
topn <- rank(-d$dNdS, ties.method = "first") <= sum(pos0)
rules <- rbind(
  rule("dNdS>1 (default)", pos0),
  rule("dNdS>0.5", d$dNdS > 0.5), rule("dNdS>0.75", d$dNdS > 0.75), rule("dNdS>1.5", d$dNdS > 1.5),
  rule("top N by ratio, same N", topn),
  rule("dNdS>1, length>=150 aa", pos0, d$q_len >= 150),
  rule("dNdS>1, dS>=0.05", pos0, d$dS >= 0.05),
  rule("dNdS>1, dS<=1", pos0, d$dS <= 1),
  rule("dN>=0.05 and dNdS>1", pos0 & d$dN >= 0.05))
write.table(rules, file.path(out, "candidate_rules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(rules)

# Threshold-free comparison: rank test per term.
rk <- rank(d$dNdS)[annotated]; mA <- m[annotated, , drop = FALSE]; ok <- which(K >= 5)
sum_rank <- as.numeric(crossprod(mA[, ok, drop = FALSE], rk)); kk <- K[ok]
zs <- (sum_rank - kk * (N + 1) / 2) / sqrt(kk * (N - kk) * (N + 1) / 12)
pw <- pnorm(zs, lower.tail = FALSE); adjw <- p.adjust(pw, "BH")
wil <- data.frame(term = ids[ok], total = kk, z = zs, p = pw, p_adj = adjw)
write.table(wil[order(wil$p), ], file.path(out, "rank_test.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Rank test: %d terms, %d significant, %d of the %d Fisher terms with at least 5 genes among them\n",
            length(ok), sum(adjw <= alpha), sum(ids[ok][adjw <= alpha] %in% base), sum(base %in% ids[ok])))
