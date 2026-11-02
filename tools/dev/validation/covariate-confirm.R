# Runs the installed covariate_enrichment() on the labels enrichment-methods.R uses, to confirm its
# calibration with the shipped code and seeds not used to choose the method.
#   Rscript covariate-confirm.R <dnds_annot.tsv> <labels> <out.tsv>
# <labels> is a sim-genome.R table, or permute:S, rotate:S, strata5:S, strata10:S.
# Writes the format enrichment-methods-summary.py reads. "fisher" is BH on the table's fisher_p.
suppressMessages(library(dndsR))
args <- commandArgs(TRUE)
annot <- args[1]; labels <- args[2]; out <- args[3]
d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d <- d[!duplicated(d$query_id), c("query_id", "subject_id", "dN", "dNdS", "dS", "q_len", "q_go", "q_ipr", "q_gff_seqname", "q_gff_start")]
d$dNdS <- suppressWarnings(as.numeric(d$dNdS)); d$dS <- suppressWarnings(as.numeric(d$dS))
d <- d[is.finite(d$dNdS) & d$dNdS < 10 & is.finite(d$dS) & d$dS > 0 & is.finite(d$q_len), ]
ntile <- function(x, k) cut(rank(x, ties.method = "first"), k, labels = FALSE)
if (grepl("\\.gz$", labels)) {
  s <- read.table(gzfile(labels), sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  s <- s[is.finite(s$dNdS) & s$dNdS < 10 & is.finite(s$dS) & s$dS > 0, ]
  d <- d[match(s$query_id, d$query_id), ]
  d$dNdS <- s$dNdS; d$dS <- s$dS; d$dN <- s$dN
  d <- d[!is.na(d$query_id), ]
} else {
  p <- strsplit(labels, ":", fixed = TRUE)[[1]]
  RNGversion("4.0.0"); set.seed(as.integer(p[2]))
  d <- d[order(d$q_gff_seqname, d$q_gff_start, d$query_id), ]
  n <- nrow(d)
  if (p[1] == "permute") d$dNdS <- sample(d$dNdS)
  if (p[1] == "rotate") { k <- sample(100:(n - 100), 1); d$dNdS <- d$dNdS[c((k + 1):n, 1:k)] }
  if (startsWith(p[1], "strata")) {
    k <- as.integer(sub("strata", "", p[1]))
    r <- lengths(strsplit(ifelse(is.na(d$q_go), "", d$q_go), "[;,|]"))
    rich <- ifelse(r == 0, 0L, ifelse(r <= stats::median(r[r > 0]), 1L, 2L))
    for (s in split(seq_len(n), interaction(ntile(d$q_len, k), ntile(d$dS, k), rich, drop = TRUE))) {
      d$dNdS[s] <- d$dNdS[s[sample.int(length(s))]]
    }
  }
}
w <- tempfile("cov"); dir.create(w); f <- file.path(w, "N_dnds_annot.tsv")
write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
files <- suppressMessages(covariate_enrichment(dnds_annot_file = f, sides = "query", min_pos = 0, verbose = FALSE))
rows <- list()
for (path in files) {
  t <- read.delim(path, quote = "", stringsAsFactors = FALSE)
  family <- t$term_type[1]
  for (method in c("covariate", "fisher")) {
    raw <- if (method == "fisher") t$fisher_p else t$p_value
    q <- p.adjust(raw, "BH"); hit <- which(q <= 0.05)
    rows[[length(rows) + 1L]] <- data.frame(family = family, method = method, term = "(family)", n_tested = nrow(t), positive = paste(unique(t$pos_total), collapse = ","),
                                            size = NA, term_positive = NA, expected = NA, raw_p = NA, p_adj = NA, n_raw05 = sum(raw < 0.05))
    if (length(hit)) rows[[length(rows) + 1L]] <- data.frame(family = family, method = method, term = t$term[hit], n_tested = nrow(t), positive = NA, size = t$total_count[hit],
                                                             term_positive = t$pos_count[hit], expected = signif(t$expected[hit], 4), raw_p = signif(raw[hit], 4), p_adj = signif(q[hit], 4), n_raw05 = NA)
  }
}
write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE, row.names = FALSE)
unlink(w, recursive = TRUE)
