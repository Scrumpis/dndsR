# Enrichment-correction study summary (docs/dev/enrichment-correction-study.md; rules fixed there).
#   Rscript --vanilla enrich-corr-summary.R <runs dir> <annot dir map tsv: pair annot custom> <out prefix> [phase: D|C]
# For every cell (pair, kind, scenario, rate or threshold, divergence) x family x method:
#   complete nulls (no spiked term): share of runs with any significant term (every term is null), Wilson 95%;
#   planted cells: per run, a significant term is TRUE if spiked, RELATED if its genes overlap the spiked
#   genes beyond chance (hypergeometric p < 0.001), otherwise FALSE; FDP = false / significant (0 if none),
#   mean FDP with a bootstrap 95% interval over runs; power = share of the 5 spiked terms found (IPR family).
# Phase D applies the selection rule; phase C the confirmatory criteria (C1-C3) and the domain diagnostic.
a <- commandArgs(TRUE); runs <- a[1]; map <- utils::read.delim(a[2], stringsAsFactors = FALSE); prefix <- a[3]; phase <- if (length(a) >= 4) a[4] else "D"
wilson <- function(k, n, z = 1.96) { if (n == 0) return(c(NA, NA)); p <- k / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(max(0, c0 - h), min(1, c0 + h)) }
split_terms <- function(x) strsplit(ifelse(is.na(x), "", x), "[;,|]")
members_cache <- list()
members_of <- function(pair) {
  if (!is.null(members_cache[[pair]])) return(members_cache[[pair]])
  m <- map[map$pair == pair, ]
  d <- utils::read.table(m$annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  d <- d[!duplicated(d$query_id), ]
  ipr <- split_terms(d$q_ipr); go <- split_terms(d$q_go)
  cs <- utils::read.delim(m$custom, stringsAsFactors = FALSE)
  anc <- utils::read.delim(m$go_ancestors, stringsAsFactors = FALSE); up <- split(anc$ancestor, anc$term)
  goall <- lapply(go, function(t) unique(unlist(up[t], use.names = FALSE)))
  res <- list(IPR = split(rep(d$query_id, lengths(ipr)), unlist(ipr)), GO = split(rep(d$query_id, lengths(goall)), unlist(goall)),
              CUSTOM = split(cs$gene, cs$set), n = nrow(d))
  members_cache[[pair]] <<- res; res
}
cells <- list(); planted <- list()
for (cd in list.dirs(runs, recursive = TRUE)) {
  files <- list.files(cd, pattern = "^terms_[0-9]+\\.tsv$", full.names = TRUE); if (!length(files)) next
  pair <- basename(dirname(cd)); cell <- basename(cd)
  for (f in files) {
    seed <- sub("^terms_([0-9]+)\\.tsv$", "\\1", basename(f))
    lab <- file.path(cd, sprintf("labels_%s.tsv.gz", seed)); h <- readLines(gzfile(lab), n = 1L)
    spiked_terms <- strsplit(sub("^.*spiked_terms=([^ ]*).*$", "\\1", h), ",")[[1]]; spiked_terms <- spiked_terms[nzchar(spiked_terms)]
    r <- utils::read.delim(f, stringsAsFactors = FALSE)
    fam <- r[r$term == "(family)", ]; hits <- r[r$term != "(family)", ]
    if (length(spiked_terms)) {
      lb <- utils::read.table(gzfile(lab), sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
      sg <- lb$query_id[lb$spiked]; mem <- members_of(pair); N <- mem$n
    }
    for (i in seq_len(nrow(fam))) {
      x <- fam[i, ]; h2 <- hits[hits$family == x$family & hits$threshold == x$threshold & hits$method == x$method, ]
      key <- paste(pair, cell, x$family, x$threshold, x$method, sep = "|")
      if (!length(spiked_terms)) {
        cells[[length(cells) + 1L]] <- data.frame(key = key, pair = pair, cell = cell, family = x$family, threshold = x$threshold, method = x$method,
          any_false = nrow(h2) > 0, diag_p = x$diag_p)
      } else {
        status <- vapply(h2$term, function(t) {
          if (t %in% spiked_terms) return("true")
          g <- mem[[x$family]][[t]]; if (is.null(g)) return("false")
          k <- sum(g %in% sg); if (stats::phyper(k - 1, length(sg), N - length(sg), length(g), lower.tail = FALSE) < 0.001) "related" else "false"
        }, "")
        planted[[length(planted) + 1L]] <- data.frame(key = key, pair = pair, cell = cell, family = x$family, threshold = x$threshold, method = x$method,
          sig = nrow(h2), false = sum(status == "false"), fdp = if (nrow(h2)) sum(status == "false") / nrow(h2) else 0,
          power = if (x$family == "IPR") mean(spiked_terms %in% h2$term) else NA)
      }
    }
  }
}
w <- function(x, what) utils::write.table(x, paste0(prefix, "_", what, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
nulls <- if (length(cells)) do.call(rbind, cells) else NULL
pl <- if (length(planted)) do.call(rbind, planted) else NULL
if (!is.null(nulls)) {
  ns <- do.call(rbind, lapply(split(nulls, nulls$key), function(x) { ci <- wilson(sum(x$any_false), nrow(x))
    data.frame(x[1, c("pair", "cell", "family", "threshold", "method")], runs = nrow(x), any_false = sum(x$any_false), rate = mean(x$any_false), lo = ci[1], hi = ci[2]) }))
  w(ns, "nulls")
}
if (!is.null(pl)) {
  ps <- do.call(rbind, lapply(split(pl, pl$key), function(x) { set.seed(1); bs <- replicate(2000, mean(sample(x$fdp, replace = TRUE)))
    data.frame(x[1, c("pair", "cell", "family", "threshold", "method")], runs = nrow(x), mean_fdp = mean(x$fdp), fdp_lo = stats::quantile(bs, 0.025, names = FALSE),
      fdp_hi = stats::quantile(bs, 0.975, names = FALSE), power = mean(x$power)) }))
  w(ps, "planted")
}
if (phase == "D" && !is.null(nulls) && !is.null(pl)) {
  # Selection rule: complete-null rate <= 5% (point) in every null cell at dN/dS > 1 (exact cells at the
  # observed rate, estimator cells at threshold 1); then highest median power relative to the oracle in
  # the planted exact-label cells; ties within 2 points to the simpler model.
  at1 <- ns[(grepl("^exact-.*-obs-", ns$cell) | (grepl("^estimator-", ns$cell) & ns$threshold == "1")), ]
  cand <- setdiff(unique(ns$method), "oracle")
  eligible <- cand[vapply(cand, function(m) all(at1$rate[at1$method == m] <= 0.05), TRUE)]
  pe <- ps[grepl("^exact-", ps$cell) & ps$family == "IPR", ]
  rel <- vapply(cand, function(m) { x <- merge(pe[pe$method == m, c("cell", "power")], pe[pe$method == "oracle", c("cell", "power")], by = "cell")
    stats::median(x$power.x / pmax(x$power.y, 1e-9)) }, 0)
  simplicity <- c("pb_floor", "pb_len", "pb_rich", "pb_floor05", "pb_floor20", "pb_gam", "goglm", "goseq_len", "goseq_both", "cmh_10x10", "fisher")
  pick <- if (length(eligible)) { best <- max(rel[eligible]); near <- eligible[rel[eligible] >= best - 0.02]; intersect(simplicity, near)[1] } else NA
  sel <- data.frame(method = cand, eligible = cand %in% eligible, worst_null_rate_at_1 = vapply(cand, function(m) max(at1$rate[at1$method == m]), 0),
                    median_power_vs_oracle = rel, chosen = cand == pick)
  w(sel, "selection"); print(sel, row.names = FALSE, digits = 3)
  # The same rule at q0.95 (exact cells at rate 0.05, estimator cells at q0.95): Q3
  atq <- ns[(grepl("^exact-.*-0.05-", ns$cell) | (grepl("^estimator-", ns$cell) & ns$threshold == "q0.95")), ]
  q3 <- data.frame(method = cand, worst_null_rate_q095 = vapply(cand, function(m) if (any(atq$method == m)) max(atq$rate[atq$method == m]) else NA, 0))
  w(q3, "q095"); print(q3, row.names = FALSE, digits = 3)
}
