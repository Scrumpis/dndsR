# Enrichment-correction study (docs/dev/enrichment-correction-study.md): labels drawn with known
# probabilities on a real gene table, so every term's truth is exact and no estimator is involved.
# Extends enrichment-direct-labels.R (study N, unchanged).
#   Rscript enrich-corr-labels.R <dnds_annot.tsv> <scenario> <rate> <seed> <out.tsv.gz>
# rate: "obs" (the table's own share above dN/dS 1), or a share such as 0.05 or 0.10 (the exact-label
# counterpart of the quantile thresholds q0.95 and q0.9).
# Scenarios. eta is a smooth logit in log length and log dS (quadratics and their product) fitted to the
# observed labels, a basis unlike the tests'; every p is rescaled to the target rate on the logit scale.
#   N0          constant p
#   N1_<s>      logit p = c + s * eta, s in 0.5, 1, 2 (graded strength of confounding)
#   N1x         eta plus 0.5 * z(log length) * z(log dS) (a stronger interaction)
#   N1r         eta plus 0.5 * z(log number of GO and InterPro terms) (annotation density)
#   N1c         eta, with dS replaced by a rank-preserving transform that raises its correlation with
#               length (Spearman about 0.6); the tests see the transformed dS
#   TL          eta, and the genes of 5 InterPro terms (20-300 genes) get dS halved before eta is
#               evaluated; the tests see it. Their call rate rises only through dS: null terms
#   NOISY       eta on the true covariates; the tests see length x exp(N(0, 0.2)) and dS x exp(N(0, 0.3))
#   <base>s_<F> as <base>, then a share F of the genes of 5 InterPro terms (20-300 genes) are made positive
#   N1sds_<F>   N1_1 with the 5 spiked terms' genes' dS halved first (signal together with low dS)
# Output columns: query_id, dNdS (2 if positive, else 0.5), dS and q_len as the tests must see them,
# spiked, p_true (the generating probability).
a <- commandArgs(TRUE)
annot <- a[1]; scenario <- a[2]; rate_arg <- a[3]; seed <- as.integer(a[4]); out <- a[5]
RNGversion("4.0.0"); set.seed(seed)
d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d$dNdS <- suppressWarnings(as.numeric(d$dNdS)); d$dS <- suppressWarnings(as.numeric(d$dS))
d <- d[!duplicated(d$query_id), ]
d <- d[is.finite(d$dNdS) & d$dNdS < 10 & is.finite(d$dS) & d$dS > 0 & is.finite(d$q_len), ]
y_obs <- d$dNdS > 1
rate <- if (rate_arg == "obs") mean(y_obs) else as.numeric(rate_arg)
z <- function(x) (x - mean(x)) / stats::sd(x)
split_terms <- function(x) strsplit(ifelse(is.na(x), "", x), "[;,|]")
ln <- log(d$q_len); lds <- log(d$dS)
fit <- suppressWarnings(stats::glm(y_obs ~ ln + I(ln^2) + lds + I(lds^2) + ln:lds, family = stats::binomial()))
eta_of <- function(ln, lds) stats::predict(fit, newdata = data.frame(ln = ln, lds = lds), type = "link")
to_rate <- function(e) { c0 <- stats::uniroot(function(c) mean(stats::plogis(c + e)) - rate, c(-50, 50))$root; stats::plogis(c0 + e) }
terms_ipr <- split_terms(d$q_ipr)
members <- split(rep(seq_along(terms_ipr), lengths(terms_ipr)), unlist(terms_ipr, use.names = FALSE))
pool <- sort(names(lengths(members))[lengths(members) >= 20 & lengths(members) <= 300])
sds <- grepl("^N1sds_[0-9.]+$", scenario)
if (sds) {
  base <- "N1"; spike_share <- as.numeric(sub("^N1sds_", "", scenario))
} else if (grepl("^(N0|N1)s_[0-9.]+$", scenario)) {
  base <- sub("s_.*$", "", scenario); spike_share <- as.numeric(sub("^.*s_", "", scenario))
} else {
  base <- scenario; spike_share <- 0
}
chosen <- character(); ds_seen <- d$dS; len_seen <- d$q_len
if (base == "TL" || sds) {
  chosen <- sample(pool, 5); halve <- unique(unlist(members[chosen])); ds_seen[halve] <- ds_seen[halve] * 0.5
}
if (base == "N1c") {
  r <- 0.6 * rank(ln) + 0.4 * rank(lds) + stats::rnorm(length(ln), 0, 0.15 * length(ln))
  ds_seen <- sort(d$dS)[rank(r, ties.method = "first")]
}
lds_seen <- log(ds_seen)
p <- switch(sub("_[0-9.]+$", "", base),
  N0 = rep(rate, nrow(d)),
  N1 = { s <- if (grepl("^N1_", base)) as.numeric(sub("^N1_", "", base)) else 1; to_rate(s * eta_of(ln, lds_seen)) },
  N1x = to_rate(eta_of(ln, lds) + 0.5 * z(ln) * z(lds)),
  N1r = to_rate(eta_of(ln, lds) + 0.5 * z(log1p(lengths(split_terms(d$q_go)) + lengths(terms_ipr)))),
  N1c = to_rate(eta_of(ln, lds_seen)),
  TL = to_rate(eta_of(ln, lds_seen)),
  NOISY = to_rate(eta_of(ln, lds)),
  stop("unknown scenario ", scenario))
if (base == "NOISY") { len_seen <- d$q_len * exp(stats::rnorm(nrow(d), 0, 0.2)); ds_seen <- d$dS * exp(stats::rnorm(nrow(d), 0, 0.3)) }
y <- stats::runif(nrow(d)) < p
spiked <- rep(FALSE, nrow(d)); spiked_terms <- character()
if (spike_share > 0) {
  spiked_terms <- if (length(chosen)) chosen else sample(pool, 5)
  for (t in spiked_terms) { i <- members[[t]]; j <- i[sample.int(length(i), max(1L, round(spike_share * length(i))))]; y[j] <- TRUE; spiked[j] <- TRUE }
}
res <- data.frame(query_id = d$query_id, dNdS = ifelse(y, 2, 0.5), dS = signif(ds_seen, 6), q_len = signif(len_seen, 6),
                  spiked = spiked, p_true = signif(p, 6))
con <- gzfile(out, "w")
writeLines(paste0("# kind=exact scenario=", scenario, " rate=", signif(rate, 4), " seed=", seed, " spiked_terms=", paste(spiked_terms, collapse = ","),
                  " null_shifted_terms=", if (base == "TL") paste(chosen, collapse = ",") else ""), con)
utils::write.table(res, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
