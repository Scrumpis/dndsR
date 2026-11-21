# Study N (docs/dev/enrichment-error-study.md): labels drawn directly, with known probabilities,
# on a real gene table (lengths, dS, terms kept), so that the truth of every term is known exactly
# and no estimator is involved. Writes a table in sim-genome.R's format, so enrichment-methods.R
# and enrichment-methods-summary.py read it unchanged.
#   Rscript enrichment-direct-labels.R <dnds_annot.tsv> <scenario> <seed> <out.tsv.gz>
# Scenarios (P(positive) for gene i is p_i; the positive rate matches the observed table's):
#   N0            p_i constant: labels exchangeable, every term null in every sense
#   N1            logit p_i a smooth function of log length and log dS, fitted to the observed labels
#                 with a basis unlike covariate_enrichment's (quadratics and their product), no term effect
#   N0s_F, N1s_F  as N0 or N1, then a share F of the genes of 5 InterPro terms (20-300 genes, chosen
#                 by the seed) become positive, as sim-genome.R's spikes do
# Columns: query_id, dNdS (2 if positive, else 0.5), dS and dN (observed, unchanged), spiked, p_true.
a <- commandArgs(TRUE)
annot <- a[1]; scenario <- a[2]; seed <- as.integer(a[3]); out <- a[4]
RNGversion("4.0.0"); set.seed(seed)
d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d$dNdS <- suppressWarnings(as.numeric(d$dNdS)); d$dS <- suppressWarnings(as.numeric(d$dS))
d <- d[!duplicated(d$query_id), ]
d <- d[is.finite(d$dNdS) & d$dNdS < 10 & is.finite(d$dS) & d$dS > 0 & is.finite(d$q_len), ]
y_obs <- d$dNdS > 1
base <- sub("s_.*$", "", scenario)
if (base == "N0") {
  p <- rep(mean(y_obs), nrow(d))
} else if (base == "N1") {
  ln <- log(d$q_len); lds <- log(d$dS)
  fit <- suppressWarnings(stats::glm(y_obs ~ ln + I(ln^2) + lds + I(lds^2) + ln:lds, family = stats::binomial()))
  p <- pmin(pmax(stats::fitted(fit), 1e-4), 0.5)
  p <- p * mean(y_obs) / mean(p) # the same expected number of positives as N0
} else stop("unknown scenario ", scenario)
y <- stats::runif(nrow(d)) < p
spiked_terms <- character(); spiked <- rep(FALSE, nrow(d))
if (grepl("s_", scenario, fixed = TRUE)) {
  share <- as.numeric(sub("^.*s_", "", scenario))
  terms <- strsplit(ifelse(is.na(d$q_ipr), "", d$q_ipr), "[;,|]")
  members <- split(rep(seq_along(terms), lengths(terms)), unlist(terms, use.names = FALSE))
  sizes <- lengths(members)
  pool <- sort(names(sizes)[sizes >= 20 & sizes <= 300])
  spiked_terms <- sample(pool, 5)
  for (t in spiked_terms) {
    i <- members[[t]]; k <- max(1L, round(share * length(i)))
    j <- i[sample.int(length(i), k)]; y[j] <- TRUE; spiked[j] <- TRUE
  }
}
res <- data.frame(query_id = d$query_id, dNdS = ifelse(y, 2, 0.5), dS = d$dS, dN = d$dS * ifelse(y, 2, 0.5),
                  spiked = spiked, p_true = signif(p, 6))
con <- gzfile(out, "w")
writeLines(paste0("# scenario=", scenario, " seed=", seed, " spiked_terms=", paste(spiked_terms, collapse = ",")), con)
utils::write.table(res, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
