# Enrichment-correction study (docs/dev/enrichment-correction-study.md): every candidate test on one
# labelled table. Base R, splines and mgcv; GOseq runs through a helper in a separate R (envs/comparators).
#   Rscript enrich-corr-methods.R <dnds_annot.tsv> <go_ancestors.tsv> <custom_sets.tsv> <labels.tsv.gz> <thresholds> <out.tsv> [goseq Rscript]
# labels: an enrich-corr-labels.R (exact) or sim-genome.R (estimator) table. thresholds: comma list of
# "1", "q0.9", "q0.95"; exact-label tables use only "1" (their rate is set by the generator).
# Families: GO (BP, MF, CC, ancestors propagated, one BH family across them, as go_enrichment's
# default scope), IPR, CUSTOM (custom_sets.tsv; min_total 5 so sparse sets are tested; others 10).
# Writes one "(family)" row per family, threshold and method (n_tested, positives, the domain
# diagnostic's p-value), and one row per term with BH p_adj <= 0.05 (raw p, p_adj, size, positives, expected).
suppressMessages({ library(splines); library(mgcv) })
a <- commandArgs(TRUE)
annot <- a[1]; ancestors <- a[2]; custom <- a[3]; labels <- a[4]; thresholds <- strsplit(a[5], ",")[[1]]; out <- a[6]
goseq_rscript <- if (length(a) >= 7) a[7] else ""
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))

d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d$dNdS <- suppressWarnings(as.numeric(d$dNdS)); d$dS <- suppressWarnings(as.numeric(d$dS))
d <- d[!duplicated(d$query_id), c("query_id", "dNdS", "dS", "q_len", "q_go", "q_ipr")]
s <- read.table(gzfile(labels), sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
head_line <- readLines(gzfile(labels), n = 1L)
kind <- if (grepl("kind=exact", head_line)) "exact" else "estimator"
d <- d[match(s$query_id, d$query_id), ]; ok <- !is.na(d$query_id); d <- d[ok, ]; s <- s[ok, ]
d$dNdS <- s$dNdS; d$dS <- s$dS
if (!is.null(s$q_len)) d$q_len <- s$q_len
p_true <- if (!is.null(s$p_true)) s$p_true else rep(NA_real_, nrow(d))
keep <- is.finite(d$dNdS) & d$dNdS < 10 & is.finite(d$dS) & d$dS > 0 & is.finite(d$q_len) & d$q_len > 0
d <- d[keep, ]; p_true <- p_true[keep]
if (kind == "exact") thresholds <- "1"

split_terms <- function(x) lapply(strsplit(ifelse(is.na(x), "", x), "[;,|]"), function(t) unique(t[nzchar(t) & t != "NA"]))
anc <- read.table(ancestors, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
direct <- split_terms(d$q_go)
go_terms <- lapply(c("BP", "MF", "CC"), function(o) { up <- split(anc$ancestor[anc$ontology == o], anc$term[anc$ontology == o])
  lapply(direct, function(t) unique(unlist(up[t], use.names = FALSE))) })
cs <- read.table(custom, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cs_by_gene <- split(cs$set, cs$gene)
custom_terms <- lapply(d$query_id, function(g) { x <- cs_by_gene[[g]]; if (is.null(x)) character() else x })
richness <- lengths(direct) + lengths(split_terms(d$q_ipr))

pb_tail <- function(p, k) { if (k <= 0L) return(1); if (k > length(p)) return(0)
  f <- c(1, numeric(k)); for (q in p) { top <- f[k + 1L] + f[k] * q; f[-1L] <- f[-1L] * (1 - q) + f[-(k + 1L)] * q; f[1L] <- f[1L] * (1 - q); f[k + 1L] <- top }
  min(1, f[k + 1L]) }
ntile <- function(x, k) cut(rank(x, ties.method = "first"), k, labels = FALSE)
stratified <- function(y, idx, bins) {
  obs <- sum(y[idx]); if (obs == 0L) return(1)
  sn <- table(bins); sk <- tapply(y, bins, sum); m <- table(bins[idx]); tail <- c(1, numeric(obs))
  for (st in names(m)[m > 0L]) {
    lo <- max(0L, m[[st]] + sk[[st]] - sn[[st]]); hi <- min(m[[st]], sk[[st]]); if (hi == 0L) next
    h <- stats::dhyper(lo:hi, sk[[st]], sn[[st]] - sk[[st]], m[[st]]); new <- numeric(obs + 1L)
    for (v in seq_along(h)) { j <- pmin(seq_len(obs + 1L) - 1L + (lo + v - 1L), obs) + 1L; add <- tail * h[v]; for (u in seq_along(j)) new[j[u]] <- new[j[u]] + add[u] }
    tail <- new }
  min(1, tail[obs + 1L]) }

methods <- c("fisher", "pb_floor", "pb_floor05", "pb_floor20", "pb_gam", "pb_rich", "pb_len", "goglm", "cmh_10x10", "goseq_len", "goseq_both", "oracle")

# One universe: genes with at least one term of the family; returns a matrix of raw p-values per term.
test_universe <- function(terms, min_total, thr) {
  has <- lengths(terms) > 0L
  dn <- d$dNdS[has]; cut <- if (thr == "1") 1 else stats::quantile(dn, as.numeric(sub("^q", "", thr)), names = FALSE, type = 7)
  y <- dn > cut; ln <- log(d$q_len[has]); lds <- log(d$dS[has]); lr <- log1p(richness[has]); pt <- p_true[has]; terms <- terms[has]
  members <- split(rep(seq_along(terms), lengths(terms)), unlist(terms, use.names = FALSE)); members <- members[lengths(members) >= min_total]
  N <- length(y); K <- sum(y); if (!length(members) || K == 0L) return(NULL)
  fitp <- function(f) stats::fitted(suppressWarnings(stats::glm(f, family = stats::binomial())))
  p_add3 <- fitp(y ~ ns(ln, 3) + ns(lds, 3)); rate <- K / N
  p_gam <- tryCatch(stats::fitted(mgcv::gam(y ~ te(ln, lds), family = stats::binomial(), method = "REML")), error = function(e) p_add3)
  p_rich <- fitp(y ~ ns(ln, 3) + ns(lds, 3) + ns(lr, 2)); p_len <- fitp(y ~ ns(ln, 3))
  probs <- list(pb_floor = pmax(p_add3, 0.1 * rate), pb_floor05 = pmax(p_add3, 0.05 * rate), pb_floor20 = pmax(p_add3, 0.2 * rate),
                pb_gam = pmax(p_gam, 0.1 * rate), pb_rich = pmax(p_rich, 0.1 * rate), pb_len = pmax(p_len, 0.1 * rate))
  # Domain diagnostic: does the call depend on length and dS? (likelihood-ratio test against a constant)
  null0 <- suppressWarnings(stats::glm(y ~ 1, family = stats::binomial())); full <- suppressWarnings(stats::glm(y ~ ns(ln, 3) + ns(lds, 3), family = stats::binomial()))
  diag_p <- stats::pchisq(null0$deviance - full$deviance, null0$df.residual - full$df.residual, lower.tail = FALSE)
  # GOglm: one-sided logistic score test of the term indicator, given ns(log length, 3) + ns(log dS, 3)
  X0 <- stats::model.matrix(~ ns(ln, 3) + ns(lds, 3)); W <- p_add3 * (1 - p_add3); r <- y - p_add3
  XtWX_inv <- solve(crossprod(X0 * sqrt(W)))
  bins <- interaction(ntile(ln, 10), ntile(lds, 10), drop = TRUE)
  m <- do.call(rbind, lapply(members, function(i) {
    k <- sum(y[i]); xw <- colSums(X0[i, , drop = FALSE] * W[i]); info <- sum(W[i]) - drop(t(xw) %*% XtWX_inv %*% xw)
    c(fisher = stats::phyper(k - 1L, K, N - K, length(i), lower.tail = FALSE),
      vapply(probs, function(pp) pb_tail(pp[i], k), 0),
      goglm = if (info > 1e-12) stats::pnorm(sum(r[i]) / sqrt(info), lower.tail = FALSE) else 1,
      cmh_10x10 = stratified(y, i, bins),
      goseq_len = NA, goseq_both = NA,
      oracle = if (all(is.finite(pt[i]))) pb_tail(pt[i], k) else NA,
      size = length(i), positive = k, expected = sum(probs$pb_floor[i]))
  }))
  rownames(m) <- names(members)
  if (nzchar(goseq_rscript)) {
    tmp <- tempfile("gs"); dir.create(tmp)
    utils::write.table(data.frame(gene = seq_len(N), y = as.integer(y), len = exp(ln), pfit = probs$pb_floor), file.path(tmp, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(data.frame(gene = unlist(members), term = rep(names(members), lengths(members))), file.path(tmp, "g2c.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    st <- system2(goseq_rscript, c("--vanilla", file.path(here, "enrich-corr-goseq.R"), tmp), stdout = FALSE, stderr = file.path(tmp, "err.txt"), env = c("R_LIBS_USER=", "R_LIBS=", "R_LIBS_SITE="))
    if (st != 0) message("GOseq helper failed: ", paste(readLines(file.path(tmp, "err.txt")), collapse = " "))
    if (st == 0 && file.exists(file.path(tmp, "out.tsv"))) {
      g <- utils::read.delim(file.path(tmp, "out.tsv"), stringsAsFactors = FALSE)
      m[, "goseq_len"] <- g$p_len[match(rownames(m), g$term)]; m[, "goseq_both"] <- g$p_both[match(rownames(m), g$term)]
    }
    unlink(tmp, recursive = TRUE)
  }
  list(m = m, universe = N, positive = K, diag_p = diag_p, cut = cut)
}

rows <- list()
emit <- function(family, thr, parts) {
  parts <- Filter(Negate(is.null), parts); if (!length(parts)) return()
  m <- do.call(rbind, lapply(parts, `[[`, "m"))
  for (method in methods) {
    pv <- m[, method]; if (all(is.na(pv))) next
    q <- stats::p.adjust(pv, "BH"); hit <- which(q <= 0.05)
    rows[[length(rows) + 1L]] <<- data.frame(kind = kind, family = family, threshold = thr, method = method, term = "(family)", n_tested = sum(!is.na(pv)),
      universe = paste(vapply(parts, `[[`, 1, "universe"), collapse = ","), positive = paste(vapply(parts, `[[`, 1, "positive"), collapse = ","),
      diag_p = paste(signif(vapply(parts, `[[`, 1, "diag_p"), 4), collapse = ","), cut = signif(parts[[1]]$cut, 5),
      size = NA, term_positive = NA, expected = NA, raw_p = NA, p_adj = NA)
    if (length(hit)) rows[[length(rows) + 1L]] <<- data.frame(kind = kind, family = family, threshold = thr, method = method, term = rownames(m)[hit], n_tested = NA,
      universe = NA, positive = NA, diag_p = NA, cut = NA, size = m[hit, "size"], term_positive = m[hit, "positive"], expected = signif(m[hit, "expected"], 4),
      raw_p = signif(pv[hit], 4), p_adj = signif(q[hit], 4))
  }
}
for (thr in thresholds) {
  emit("GO", thr, lapply(go_terms, test_universe, min_total = 10L, thr = thr))
  emit("IPR", thr, list(test_universe(split_terms(d$q_ipr), 10L, thr)))
  emit("CUSTOM", thr, list(test_universe(custom_terms, 5L, thr)))
}
res <- do.call(rbind, rows)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
