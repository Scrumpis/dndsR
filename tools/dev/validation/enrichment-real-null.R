# Confirms enrichment-calibration.R with the package's own functions, not the reimplementation.
#   Rscript enrichment-real-null.R <dnds_annot.tsv> <runs> <seed> <outdir>
# Reassigns dN/dS at random, so every significant term in a run is a false positive, and
# counts those from term_enrichment() and ipr_enrichment() at their defaults.
suppressMessages(library(dndsR))
a <- commandArgs(TRUE)
annot <- a[1]; runs <- as.integer(a[2]); seed <- as.integer(a[3]); out <- a[4]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)
d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d <- d[, c("query_id", "subject_id", "dN", "dS", "dNdS", "q_ipr", "q_gff_seqname", "q_gff_start")]
d$dNdS <- suppressWarnings(as.numeric(d$dNdS))
d <- d[is.finite(d$dNdS) & d$dNdS < 10, ]
rows <- list()
count <- function(files, fun, i) {
  for (p in unlist(files)) {
    r <- read.table(p, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <<- data.frame(fun = fun, run = i, file = sub("^N_", "", basename(p)), tested = nrow(r),
                                             n_sig = sum(r$p_adj <= 0.05, na.rm = TRUE))
  }
}
for (i in seq_len(runs)) {
  w <- tempfile("real_", tmpdir = out); dir.create(w)
  x <- d; x$dNdS <- sample(d$dNdS)
  f <- file.path(w, "N_dnds_annot.tsv")
  write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE)
  count(term_enrichment(dnds_annot_file = f, terms = "ipr", exclude_terms = NULL, sides = "query",
                        make_plots = FALSE, threads = 1, verbose = FALSE), "term_enrichment", i)
  count(ipr_enrichment(dnds_annot_file = f, sides = "query", make_plots = FALSE, threads = 1,
                       entries_source = "local", entries_path = system.file("extdata", "interpro_entry.list.tsv", package = "dndsR"),
                       tree_source = "none"), "ipr_enrichment", i)
  unlink(w, recursive = TRUE)
}
r <- do.call(rbind, rows)
write.table(r, file.path(out, sprintf("real_null_%d.tsv", seed)), sep = "\t", quote = FALSE, row.names = FALSE)
per_run <- aggregate(n_sig ~ fun + run, r, sum)
print(aggregate(n_sig ~ fun, per_run, function(v) c(runs = length(v), any_sig = mean(v > 0), mean_sig = mean(v))))
