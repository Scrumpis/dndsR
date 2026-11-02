# Confirms enrichment-calibration.R with the package's own functions, not the reimplementation.
#   Rscript enrichment-real-null.R <dnds_annot.tsv> <runs per null> <seed> <outdir> [InterPro tree]
# Reassigns dN/dS by each of three nulls, so every significant term in a run is a false
# positive, and counts those at the functions' defaults. With a tree, tests only the
# parent-child method. IPR_ADJUST_SCOPE sets adjust_scope for ipr_enrichment; when set,
# term_enrichment is left out. Results are written after each run, in case a job is stopped.
suppressMessages(library(dndsR))
a <- commandArgs(TRUE)
annot <- a[1]; runs <- as.integer(a[2]); seed <- as.integer(a[3]); out <- a[4]
tree <- if (length(a) > 4L) a[5] else NULL
scope <- Sys.getenv("IPR_ADJUST_SCOPE")
scope_args <- if (nzchar(scope)) list(adjust_scope = scope) else list()
fun_suffix <- if (nzchar(scope)) paste0("_", scope) else ""
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)
d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d <- d[, c("query_id", "subject_id", "dN", "dS", "dNdS", "q_len", "q_ipr", "q_gff_seqname", "q_gff_start")]
d$dNdS <- suppressWarnings(as.numeric(d$dNdS))
d <- d[is.finite(d$dNdS) & d$dNdS < 10, ]
d <- d[order(d$q_gff_seqname, d$q_gff_start, d$query_id), ]
n <- nrow(d)
ntile <- function(x, k) cut(rank(x, ties.method = "first"), k, labels = FALSE)
terms_n <- lengths(strsplit(ifelse(is.na(d$q_ipr), "", d$q_ipr), ";", fixed = TRUE))
rich <- ifelse(terms_n == 0, 0L, ifelse(terms_n <= median(terms_n[terms_n > 0]), 1L, 2L))
by_stratum <- split(seq_len(n), interaction(ntile(d$q_len, 5), ntile(d$dS, 5), rich, drop = TRUE))
draw <- list(
  permute = function(x) sample(x),
  strata = function(x) { for (s in by_stratum) x[s] <- x[s[sample.int(length(s))]]; x },
  rotate = function(x) { k <- sample(100:(n - 100), 1); x[c((k + 1):n, 1:k)] }
)
entries <- system.file("extdata", "interpro_entry.list.tsv", package = "dndsR")
results <- file.path(out, sprintf("real_null_%d.tsv", seed))
rows <- list()
count <- function(files, fun, null, i) {
  for (p in unlist(files)) {
    r <- read.table(p, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <<- data.frame(
      fun = fun, null = null, run = i, file = sub("^N_", "", basename(p)), reported = nrow(r),
      n_tested = if ("n_tested" %in% names(r)) r$n_tested[1] else NA_integer_, n_sig = sum(r$p_adj <= 0.05, na.rm = TRUE))
  }
  # A run with nothing to report still counts, as a run with no significant term.
  rows[[length(rows) + 1L]] <<- data.frame(fun = fun, null = null, run = i, file = "(run)", reported = 0L, n_tested = NA_integer_, n_sig = 0L)
}
for (null in names(draw)) for (i in seq_len(runs)) {
  w <- tempfile("real_", tmpdir = out); dir.create(w)
  x <- d; x$dNdS <- draw[[null]](d$dNdS)
  f <- file.path(w, "N_dnds_annot.tsv")
  write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE)
  ipr_args <- c(list(dnds_annot_file = f, sides = "query", make_plots = FALSE, threads = 1, auto_detect_release = FALSE,
                     entries_source = "local", entries_path = entries), scope_args)
  if (is.null(tree)) {
    if (!nzchar(scope)) {
      count(term_enrichment(dnds_annot_file = f, terms = "ipr", exclude_terms = NULL, sides = "query",
                            make_plots = FALSE, threads = 1, verbose = FALSE), "term_enrichment", null, i)
    }
    count(do.call(ipr_enrichment, c(ipr_args, list(tree_source = "none"))), paste0("ipr_enrichment", fun_suffix), null, i)
  } else {
    r <- tryCatch(do.call(ipr_enrichment, c(ipr_args, list(method = "parent_child", tree_source = "local", tree_path = tree))),
                  error = function(e) e)
    # An error must not count as a run with no significant term, so it is recorded, and such
    # runs are left out when rates are computed.
    if (inherits(r, "error")) {
      rows[[length(rows) + 1L]] <- data.frame(fun = paste0("ipr_parent_child", fun_suffix), null = null, run = i, file = "(error)",
                                              reported = 0L, n_tested = NA_integer_, n_sig = 0L)
      message(conditionMessage(r))
    } else count(r, paste0("ipr_parent_child", fun_suffix), null, i)
  }
  unlink(w, recursive = TRUE)
  write.table(do.call(rbind, rows), results, sep = "\t", quote = FALSE, row.names = FALSE)
}
