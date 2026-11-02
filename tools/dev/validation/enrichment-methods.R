# Candidate enrichment tests, compared on the same labels. Base R plus splines, except the
# optional topGO default below.
#   Rscript enrichment-methods.R <dnds_annot.tsv> <go_ancestors.tsv> <labels> <out.tsv>
# <labels> is one of:
#   observed                 the table's own dN/dS
#   <file>.tsv.gz            a sim-genome.R table
#   permute:S, rotate:S, strata5:S, strata10:S    a reassignment of the table's dN/dS with seed S;
#                            strataK shuffles within K x K bins of length and dS, by annotation richness
# Writes a row per method and family with family size, and a row per significant term.
# With DNDSR_GO_DEFAULT=1 also runs the real go_enrichment() defaults on these labels.
#
# Methods. y is dN/dS > 1. Each GO ontology has its own universe of annotated genes, as in topGO, and
# BH runs once across a family (GO: all ontologies; IPR: all entries), which is what the defaults do.
# Terms are tested if at least min_total genes in the universe carry them, independent of y.
#   fisher     Fisher's test, on all ancestors: go_enrichment(algorithm = "classic") and term_enrichment
#   pb_len     Poisson-binomial tail, with P(y = 1) from logistic regression on spline(log length)
#   pb_len_ds  The same, on spline(log length) + spline(log dS)
#   pb_add3    The same with 3 degrees of freedom each. An interaction was tried and rejected: with about
#              100 positive genes the fit separates, giving probabilities near 1e-14 and p-values near 1e-260
#   pb_floor   pb_add3, with probabilities no lower than a tenth of the positive rate, to limit extrapolation
#   pb_rich    pb_add3, and spline(log number of terms)
#   cmh_5x5, cmh_10x10  Exact stratified test: a sum of hypergeometrics over bins of length and dS
#   filter     Fisher's test on genes of at least 150 amino acids and dS >= 0.05
suppressMessages(library(splines))
args <- commandArgs(TRUE)
annot <- args[1]; ancestors <- args[2]; labels <- args[3]; out <- args[4]
min_total <- 10L

d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d$dNdS <- suppressWarnings(as.numeric(d$dNdS)); d$dS <- suppressWarnings(as.numeric(d$dS))
d <- d[!duplicated(d$query_id), c("query_id", "subject_id", "dN", "dNdS", "dS", "q_len", "q_go", "q_ipr", "q_gff_seqname", "q_gff_start")]
d <- d[is.finite(d$dNdS) & d$dNdS < 10 & is.finite(d$dS) & d$dS > 0 & is.finite(d$q_len), ]
ntile <- function(x, k) cut(rank(x, ties.method = "first"), k, labels = FALSE)
richness <- function(x) lengths(strsplit(ifelse(is.na(x), "", x), "[;,|]"))
if (grepl("\\.gz$", labels)) {
  s <- read.table(gzfile(labels), sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  s <- s[is.finite(s$dNdS) & s$dNdS < 10 & is.finite(s$dS) & s$dS > 0, ]
  d <- d[match(s$query_id, d$query_id), ]
  d$dNdS <- s$dNdS; d$dS <- s$dS; d$dN <- s$dN
  d <- d[!is.na(d$query_id), ]
} else if (labels != "observed") {
  p <- strsplit(labels, ":", fixed = TRUE)[[1]]
  RNGversion("4.0.0"); set.seed(as.integer(p[2]))
  d <- d[order(d$q_gff_seqname, d$q_gff_start, d$query_id), ]
  n <- nrow(d)
  if (p[1] == "permute") d$dNdS <- sample(d$dNdS)
  if (p[1] == "rotate") { k <- sample(100:(n - 100), 1); d$dNdS <- d$dNdS[c((k + 1):n, 1:k)] }
  if (startsWith(p[1], "strata")) {
    k <- as.integer(sub("strata", "", p[1])); r <- richness(d$q_go)
    rich <- ifelse(r == 0, 0L, ifelse(r <= stats::median(r[r > 0]), 1L, 2L))
    for (s in split(seq_len(n), interaction(ntile(d$q_len, k), ntile(d$dS, k), rich, drop = TRUE))) {
      d$dNdS[s] <- d$dNdS[s[sample.int(length(s))]]
    }
  }
}

split_terms <- function(x) strsplit(ifelse(is.na(x), "", x), "[;,|]")
anc <- read.table(ancestors, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
direct <- split_terms(d$q_go)
go_by <- lapply(c(BP = "BP", MF = "MF", CC = "CC"), function(o) {
  up <- split(anc$ancestor[anc$ontology == o], anc$term[anc$ontology == o])
  lapply(direct, function(t) unique(unlist(up[t], use.names = FALSE)))
})

poisson_binomial_tail <- function(p, k) {
  if (k <= 0L) return(1)
  f <- c(1, numeric(k)) # f[j + 1] = P(S = j); everything at k or above is kept in f[k + 1]
  for (q in p) {
    top <- f[k + 1L] + f[k] * q
    f[-1L] <- f[-1L] * (1 - q) + f[-(k + 1L)] * q
    f[1L] <- f[1L] * (1 - q)
    f[k + 1L] <- top
  }
  min(1, f[k + 1L])
}

methods <- c("fisher", "pb_len", "pb_len_ds", "pb_add3", "pb_floor", "pb_rich", "cmh_5x5", "cmh_10x10", "filter")
test_universe <- function(terms, richness_terms = terms) {
  has <- lengths(terms) > 0L
  y <- d$dNdS[has] > 1; ln <- log(d$q_len[has]); lds <- log(d$dS[has]); lr <- log(lengths(richness_terms)[has]); terms <- terms[has]
  members <- split(rep(seq_along(terms), lengths(terms)), unlist(terms, use.names = FALSE))
  members <- members[lengths(members) >= min_total]
  N <- length(y); K <- sum(y)
  fit <- function(f) stats::fitted(suppressWarnings(stats::glm(f, family = stats::binomial())))
  p_len <- fit(y ~ ns(ln, 4)); p_both <- fit(y ~ ns(ln, 4) + ns(lds, 4)); p_add3 <- fit(y ~ ns(ln, 3) + ns(lds, 3))
  p_floor <- pmax(p_add3, 0.1 * K / N)
  p_rich <- fit(y ~ ns(ln, 3) + ns(lds, 3) + ns(lr, 2))
  bin5 <- interaction(ntile(ln, 5), ntile(lds, 5), drop = TRUE)
  bin10 <- interaction(ntile(ln, 10), ntile(lds, 10), drop = TRUE)
  keep <- exp(ln) >= 150 & exp(lds) >= 0.05
  N2 <- sum(keep); K2 <- sum(y[keep])
  stratified <- function(idx, bins) {
    obs <- sum(y[idx]); if (obs == 0L) return(1)
    stratum_n <- table(bins); stratum_k <- tapply(y, bins, sum)
    m <- table(bins[idx]); tail <- c(1, numeric(obs)) # as above, absorbing at obs
    for (s in names(m)[m > 0L]) {
      lo <- max(0L, m[[s]] + stratum_k[[s]] - stratum_n[[s]]); hi <- min(m[[s]], stratum_k[[s]])
      if (hi == 0L) next
      h <- stats::dhyper(lo:hi, stratum_k[[s]], stratum_n[[s]] - stratum_k[[s]], m[[s]])
      new <- numeric(obs + 1L)
      for (v in seq_along(h)) {
        j <- pmin(seq_len(obs + 1L) - 1L + (lo + v - 1L), obs) + 1L
        add <- tail * h[v]
        for (u in seq_along(j)) new[j[u]] <- new[j[u]] + add[u]
      }
      tail <- new
    }
    min(1, tail[obs + 1L])
  }
  m <- do.call(rbind, lapply(names(members), function(t) {
    i <- members[[t]]; k <- sum(y[i]); j <- i[keep[i]]
    c(fisher = stats::phyper(k - 1L, K, N - K, length(i), lower.tail = FALSE),
      pb_len = poisson_binomial_tail(p_len[i], k), pb_len_ds = poisson_binomial_tail(p_both[i], k),
      pb_add3 = poisson_binomial_tail(p_add3[i], k), pb_floor = poisson_binomial_tail(p_floor[i], k),
      pb_rich = poisson_binomial_tail(p_rich[i], k), cmh_5x5 = stratified(i, bin5), cmh_10x10 = stratified(i, bin10),
      filter = if (length(j) == 0L) 1 else stats::phyper(sum(y[j]) - 1L, K2, N2 - K2, length(j), lower.tail = FALSE),
      size = length(i), positive = k, expected = sum(p_add3[i]))
  }))
  rownames(m) <- names(members)
  list(m = m, universe = N, positive = K)
}

summarize <- function(name, parts) {
  m <- do.call(rbind, lapply(parts, `[[`, "m"))
  res <- list()
  for (method in methods) {
    q <- stats::p.adjust(m[, method], "BH"); hit <- which(q <= 0.05)
    res[[length(res) + 1L]] <- data.frame(family = name, method = method, term = "(family)", n_tested = nrow(m),
      universe = paste(vapply(parts, `[[`, 1, "universe"), collapse = ","), positive = paste(vapply(parts, `[[`, 1, "positive"), collapse = ","),
      size = NA, term_positive = NA, expected = NA, raw_p = NA, p_adj = NA, n_raw05 = sum(m[, method] < 0.05))
    if (length(hit)) res[[length(res) + 1L]] <- data.frame(family = name, method = method, term = rownames(m)[hit], n_tested = nrow(m),
      universe = NA, positive = NA, size = m[hit, "size"], term_positive = m[hit, "positive"], expected = signif(m[hit, "expected"], 4),
      raw_p = signif(m[hit, method], 4), p_adj = signif(q[hit], 4), n_raw05 = NA)
  }
  do.call(rbind, res)
}
out_rows <- rbind(summarize("GO", lapply(go_by, test_universe)), summarize("IPR", list(test_universe(split_terms(d$q_ipr)))))

if (nzchar(Sys.getenv("DNDSR_GO_DEFAULT"))) {
  suppressMessages(library(dndsR))
  w <- tempfile("go"); dir.create(w); f <- file.path(w, "N_dnds_annot.tsv")
  write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
  files <- unlist(go_enrichment(dnds_annot_file = f, sides = "query", ontologies = c("BP", "MF", "CC"), make_plots = FALSE, threads = 1))
  tabs <- do.call(rbind, lapply(files[grepl("\\.tsv$", files)], function(p) read.table(p, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)))
  hit <- tabs[tabs$p_adj <= 0.05, ]
  out_rows <- rbind(out_rows, data.frame(family = "GO", method = "go_default", term = "(family)", n_tested = nrow(tabs), universe = NA, positive = NA,
                                         size = NA, term_positive = NA, expected = NA, raw_p = NA, p_adj = NA, n_raw05 = sum(tabs$raw_p < 0.05)))
  if (nrow(hit)) out_rows <- rbind(out_rows, data.frame(family = "GO", method = "go_default", term = hit$GO_ID, n_tested = nrow(tabs), universe = NA, positive = NA,
                                                        size = hit$annotated, term_positive = hit$significant, expected = NA, raw_p = hit$raw_p, p_adj = hit$p_adj, n_raw05 = NA))
  unlink(w, recursive = TRUE)
}
rownames(out_rows) <- NULL
write.table(out_rows, out, sep = "\t", quote = FALSE, row.names = FALSE)
