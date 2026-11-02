# Runs the real go_enrichment() with its defaults on null assignments of dN/dS to genes.
#   Rscript enrichment-go-null.R <dnds_annot.tsv> <task> <runs per null> <outdir>
# GO_ALGORITHM sets another topGO algorithm than the default. n_padj05 counts terms at
# p_adj <= 0.05 as written (the function's default correction), and n_bh05 is BH on that
# ontology's raw_p alone.
# Each run reassigns the dN/dS column, by the same three nulls as enrichment-calibration.R,
# so the only thing separating a positive gene from another is the null.
suppressMessages(library(dndsR))
a <- commandArgs(TRUE)
annot <- a[1]; task <- as.integer(a[2]); runs <- as.integer(a[3]); out <- a[4]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(1000L + task)

d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
keep <- c("query_id", "subject_id", "dN", "dS", "dNdS", "q_len", "q_go", "q_gff_seqname", "q_gff_start")
d <- d[, keep]
d$dNdS <- suppressWarnings(as.numeric(d$dNdS))
d <- d[is.finite(d$dNdS) & d$dNdS < 10, ]
d <- d[order(d$q_gff_seqname, d$q_gff_start, d$query_id), ]
n <- nrow(d)
ntile <- function(x, k) cut(rank(x, ties.method = "first"), k, labels = FALSE)
go_n <- lengths(strsplit(ifelse(is.na(d$q_go), "", d$q_go), ";", fixed = TRUE))
rich <- ifelse(go_n == 0, 0L, ifelse(go_n <= median(go_n[go_n > 0]), 1L, 2L))
by_stratum <- split(seq_len(n), interaction(ntile(d$q_len, 5), ntile(d$dS, 5), rich, drop = TRUE))
draw <- list(
  observed = function(x) x,
  permute = function(x) sample(x),
  strata = function(x) { for (s in by_stratum) x[s] <- x[s[sample.int(length(s))]]; x },
  rotate = function(x) { k <- sample(100:(n - 100), 1); x[c((k + 1):n, 1:k)] }
)

summary_rows <- list(); term_rows <- list()
for (nm in names(draw)) {
  # The observed labels need one run, from the first task only.
  iters <- if (nm == "observed") as.integer(task == 1L) else runs
  for (i in seq_len(iters)) {
    w <- tempfile(paste0("go_", task, "_"), tmpdir = out)
    dir.create(w)
    x <- d; x$dNdS <- draw[[nm]](d$dNdS)
    f <- file.path(w, "N_dnds_annot.tsv")
    write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE)
    args <- list(dnds_annot_file = f, sides = "query", make_plots = FALSE, threads = 1)
    if (nzchar(Sys.getenv("GO_ALGORITHM"))) args$algorithm <- Sys.getenv("GO_ALGORITHM")
    files <- suppressMessages(do.call(go_enrichment, args))
    for (p in unlist(files)) {
      r <- read.table(p, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
      pv <- suppressWarnings(as.numeric(r$raw_p))
      ont <- r$ontology[1]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        null = nm, run = i, task = task, ontology = ont, tested = sum(!is.na(pv)),
        n_p05 = sum(pv < 0.05, na.rm = TRUE), n_p01 = sum(pv < 0.01, na.rm = TRUE),
        n_p001 = sum(pv < 0.001, na.rm = TRUE), min_p = min(pv, na.rm = TRUE),
        n_bh05 = sum(p.adjust(pv, "BH") <= 0.05, na.rm = TRUE),
        n_padj05 = if ("p_adj" %in% names(r)) sum(r$p_adj <= 0.05, na.rm = TRUE) else NA_integer_,
        n_tested = if ("n_tested" %in% names(r)) r$n_tested[1] else NA_integer_,
        method = if ("p_adjust_method" %in% names(r)) r$p_adjust_method[1] else NA_character_)
      low <- which(pv < 0.05)
      if (length(low)) term_rows[[length(term_rows) + 1L]] <- data.frame(
        null = nm, run = i, task = task, ontology = ont, term = r$GO_ID[low], p = pv[low], p_adj = r$p_adj[low],
        annotated = r$annotated[low], significant = r$significant[low])
    }
    unlink(w, recursive = TRUE)
  }
}
write.table(do.call(rbind, summary_rows), file.path(out, sprintf("summary_%02d.tsv", task)), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, term_rows), file.path(out, sprintf("terms_%02d.tsv", task)), sep = "\t", quote = FALSE, row.names = FALSE)
