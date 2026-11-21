# Study N summary (docs/dev/enrichment-error-study.md).
#   Rscript enrichment-direct-summary.R <runs dir> <dnds_annot.tsv> <go_ancestors.tsv> <out prefix>
# Per scenario, family, and method:
#   runs, runs_any (runs with at least one significant term) and its Wilson interval. Under a
#     complete null (N0, N1) this is the family-wise error rate, and also the FDR, since every
#     discovery is false.
#   mean_fdp: false discovery proportion averaged over runs (0 in runs without discoveries), with a
#     normal interval. A discovery is false unless the term is spiked or overlaps the spiked genes
#     beyond chance (hypergeometric p < 0.001), the rule of enrichment-methods-summary.py.
#   power: share of spiked InterPro terms found (IPR family only).
#   For the test's own hypothesis (independence of label and term): rr_median, the median over
#   false discoveries of the term's true relative risk, mean p_true in the term over mean p_true in
#   the rest of its universe, and share_rr_gt_1.1. Under N0 every term has relative risk 1.
suppressMessages(library(stats))
a <- commandArgs(TRUE); runs_dir <- a[1]; annot <- a[2]; ancestors <- a[3]; prefix <- a[4]
d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d <- d[!duplicated(d$query_id), ]
split_terms <- function(x) strsplit(ifelse(is.na(x), "", x), "[;,|]")
anc <- read.table(ancestors, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
direct <- split_terms(d$q_go); names(direct) <- d$query_id
memb <- list()
for (o in c("BP", "MF", "CC")) {
  up <- split(anc$ancestor[anc$ontology == o], anc$term[anc$ontology == o])
  g <- lapply(direct, function(t) unique(unlist(up[t], use.names = FALSE)))
  g <- g[lengths(g) > 0]
  memb[[o]] <- list(universe = names(g), m = split(rep(names(g), lengths(g)), unlist(g, use.names = FALSE)))
}
ipr <- split_terms(d$q_ipr); names(ipr) <- d$query_id; ipr <- ipr[lengths(ipr) > 0]
memb$IPR <- list(universe = names(ipr), m = split(rep(names(ipr), lengths(ipr)), unlist(ipr, use.names = FALSE)))
which_family <- function(term, fam) {
  if (fam == "IPR") return("IPR")
  for (o in c("BP", "MF", "CC")) if (!is.null(memb[[o]]$m[[term]])) return(o)
  NA_character_
}
wilson <- function(k, n, z = 1.959964) {
  if (n == 0) return(c(NA, NA)); p <- k / n; den <- 1 + z^2 / n
  c(p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)), p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den
}
rows <- list(); disc <- list()
for (folder in sort(list.dirs(runs_dir, recursive = FALSE))) {
  scen <- basename(folder)
  files <- list.files(folder, "^terms_[0-9]+\\.tsv$", full.names = TRUE)
  acc <- list()
  for (f in files) {
    s <- sub("^terms_([0-9]+)\\.tsv$", "\\1", basename(f))
    tab <- read.delim(f, stringsAsFactors = FALSE)
    lab_f <- file.path(folder, paste0("sim_", s, ".tsv.gz"))
    first <- readLines(gzfile(lab_f), n = 1)
    spiked_terms <- strsplit(sub("^.*spiked_terms=", "", first), ",", fixed = TRUE)[[1]]
    lab <- read.table(gzfile(lab_f), sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
    p_true <- stats::setNames(lab$p_true, lab$query_id); spiked_genes <- lab$query_id[lab$spiked]
    for (key in unique(paste(tab$family, tab$method, sep = "|"))) {
      fam <- sub("\\|.*", "", key); meth <- sub(".*\\|", "", key)
      hits <- tab[tab$family == fam & tab$method == meth & tab$term != "(family)", ]
      n_false <- 0; found <- NA
      for (t in hits$term) {
        o <- which_family(t, fam); m <- if (is.na(o)) character() else memb[[o]]$m[[t]]
        u <- if (is.na(o)) character() else memb[[o]]$universe
        is_false <- !(t %in% spiked_terms)
        if (is_false && length(spiked_genes) && length(m)) {
          k <- length(intersect(m, spiked_genes)); sp <- length(intersect(u, spiked_genes))
          if (k > 0 && phyper(k - 1, sp, length(u) - sp, length(m), lower.tail = FALSE) < 0.001) is_false <- FALSE
        }
        if (is_false) {
          n_false <- n_false + 1
          rr <- if (length(m)) mean(p_true[m], na.rm = TRUE) / mean(p_true[setdiff(u, m)], na.rm = TRUE) else NA
          disc[[length(disc) + 1L]] <- data.frame(scenario = scen, family = fam, method = meth, seed = s, term = t, size = length(m), rr_true = rr)
        }
      }
      if (fam == "IPR" && length(spiked_terms) && nzchar(spiked_terms[1])) found <- mean(spiked_terms %in% hits$term)
      acc[[length(acc) + 1L]] <- data.frame(family = fam, method = meth, any = nrow(hits) > 0, any_false = n_false > 0,
                                            fdp = if (nrow(hits)) n_false / nrow(hits) else 0, found = found)
    }
  }
  if (!length(acc)) next
  acc <- do.call(rbind, acc)
  for (key in unique(paste(acc$family, acc$method, sep = "|"))) {
    x <- acc[paste(acc$family, acc$method, sep = "|") == key, ]
    n <- nrow(x); k <- sum(x$any_false); w <- wilson(k, n)
    se <- stats::sd(x$fdp) / sqrt(n)
    rows[[length(rows) + 1L]] <- data.frame(scenario = scen, family = x$family[1], method = x$method[1], runs = n,
      runs_any = sum(x$any), runs_any_false = k, pct_any_false = 100 * k / n, lo = 100 * w[1], hi = 100 * w[2],
      mean_fdp = 100 * mean(x$fdp), fdp_lo = 100 * (mean(x$fdp) - 1.96 * se), fdp_hi = 100 * (mean(x$fdp) + 1.96 * se),
      power = if (all(is.na(x$found))) NA else 100 * mean(x$found, na.rm = TRUE))
  }
}
res <- do.call(rbind, rows)
dd <- if (length(disc)) do.call(rbind, disc) else data.frame(scenario = character(), family = character(), method = character(), rr_true = numeric())
rr <- if (nrow(dd)) aggregate(rr_true ~ scenario + family + method, dd, function(v) c(median = median(v), gt = mean(v > 1.1))) else NULL
if (!is.null(rr)) {
  rr <- data.frame(rr[, 1:3], rr_median = rr$rr_true[, "median"], share_rr_gt_1.1 = rr$rr_true[, "gt"])
  res <- merge(res, rr, all.x = TRUE)
}
res <- res[order(res$scenario, res$family, res$method), ]
num <- vapply(res, is.numeric, logical(1)); res[num] <- lapply(res[num], function(v) round(v, 3))
write.table(res, paste0(prefix, "_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dd, paste0(prefix, "_false_discoveries.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
