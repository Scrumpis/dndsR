# How much does isoform choice change dN/dS? For gene pairs where at least one gene has several
# isoforms, compares every isoform combination on one independently built alignment.
#   Rscript isoform-effects.R <work dir> <comparison> <query FASTA> <query GFF> <subject FASTA> <subject GFF> <gene pairs> <output dir>
# Needs dndsR's dN/dS table and CDS QC tables in <work dir>/<comparison>.
suppressPackageStartupMessages({
  library(Biostrings); library(pwalign); library(GenomicFeatures); library(txdbmaker); library(Rsamtools)
})
a <- commandArgs(TRUE)
work <- a[1]; cmp <- a[2]; qfa <- a[3]; qgff <- a[4]; sfa <- a[5]; sgff <- a[6]; n_pairs <- as.integer(a[7]); out <- a[8]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
RNGversion("4.0.0"); set.seed(7)
tmp <- tempfile("iso"); dir.create(tmp)

load_isoforms <- function(fa, gff) {
  txdb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(gff, format = "gff3")))
  cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
  map <- suppressMessages(AnnotationDbi::select(txdb, keys = names(cds), keytype = "TXNAME", columns = "GENEID"))
  seqs <- extractTranscriptSeqs(FaFile(fa), cds)
  list(seq = setNames(as.character(seqs), names(cds)), gene = map$GENEID[match(names(cds), map$TXNAME)], tx = names(cds))
}
q <- load_isoforms(qfa, qgff); s <- load_isoforms(sfa, sgff)
d <- read.delim(file.path(work, cmp, paste0(cmp, "_dnds.tsv")))[, c("query_id", "subject_id", "dN", "dS", "dNdS")]

valid <- function(x) {
  if (nchar(x) %% 3L != 0L || nchar(x) < 30L) return(FALSE)
  p <- as.character(suppressWarnings(translate(DNAString(x), if.fuzzy.codon = "X")))
  !grepl("*", sub("\\*$", "", p), fixed = TRUE)
}
strip_stop <- function(x) {
  p <- as.character(suppressWarnings(translate(DNAString(x), if.fuzzy.codon = "X")))
  if (endsWith(p, "*")) substr(x, 1, nchar(x) - 3) else x
}
pair_rate <- function(x, y) {
  x <- strip_stop(x); y <- strip_stop(y)
  al <- pwalign::pairwiseAlignment(AAString(as.character(translate(DNAString(x), if.fuzzy.codon = "X"))),
                                   AAString(as.character(translate(DNAString(y), if.fuzzy.codon = "X"))),
                                   substitutionMatrix = "BLOSUM62", type = "global")
  back <- function(aln, cds) {
    cod <- substring(cds, seq(1, nchar(cds) - 2, 3), seq(3, nchar(cds), 3))
    p <- strsplit(aln, "")[[1]]; res <- character(length(p)); j <- 1L
    for (k in seq_along(p)) if (p[k] == "-") res[k] <- "---" else { res[k] <- cod[j]; j <- j + 1L }
    paste(res, collapse = "")
  }
  f <- file.path(tmp, "p.fa")
  writeLines(c(">q", back(as.character(pwalign::alignedPattern(al)), x), ">s", back(as.character(pwalign::alignedSubject(al)), y)), f)
  r <- orthologr::substitutionrate(f, est.method = "Comeron", quiet = TRUE)
  c(dN = r$dN, dS = r$dS, ratio = r$dNdS)
}

n_iso <- function(g, set) sum(set$gene == g, na.rm = TRUE)
d$q_iso <- vapply(d$query_id, n_iso, 1L, set = q); d$s_iso <- vapply(d$subject_id, n_iso, 1L, set = s)
multi <- d[d$q_iso > 1L | d$s_iso > 1L, ]
cat("Pairs:", nrow(d), " with an isoform choice:", nrow(multi), "\n")
take <- multi[sample.int(nrow(multi), min(n_pairs, nrow(multi))), ]
rows <- list()
for (i in seq_len(nrow(take))) {
  g <- take[i, ]
  qt <- q$tx[which(q$gene == g$query_id)]; st <- s$tx[which(s$gene == g$subject_id)]
  qt <- qt[vapply(q$seq[qt], valid, TRUE)]; st <- st[vapply(s$seq[st], valid, TRUE)]
  if (!length(qt) || !length(st)) next
  # Bound the combinations in genes with many isoforms.
  grid <- expand.grid(q = qt, s = st, stringsAsFactors = FALSE)
  if (nrow(grid) > 16L) grid <- grid[sample.int(nrow(grid), 16L), ]
  for (k in seq_len(nrow(grid))) {
    r <- tryCatch(pair_rate(q$seq[[grid$q[k]]], s$seq[[grid$s[k]]]), error = function(e) NULL)
    if (is.null(r)) next
    rows[[length(rows) + 1L]] <- data.frame(pair = i, q = grid$q[k], s = grid$s[k], qlen = nchar(q$seq[[grid$q[k]]]), slen = nchar(s$seq[[grid$s[k]]]), t(r), reported = g$dNdS)
  }
}
x <- do.call(rbind, rows)
utils::write.table(x, file.path(out, "isoform_combinations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
per <- do.call(rbind, lapply(split(x, x$pair), function(p) {
  ok <- p[is.finite(p$ratio), ]
  if (nrow(ok) < 2L) return(NULL)
  longest <- ok[which.max(ok$qlen + ok$slen), ]
  data.frame(pair = p$pair[1], combos = nrow(ok), min = min(ok$ratio), max = max(ok$ratio), longest = longest$ratio,
             median_other = stats::median(ok$ratio), range = max(ok$ratio) - min(ok$ratio),
             crosses_1 = min(ok$ratio) <= 1 & max(ok$ratio) > 1)
}))
utils::write.table(per, file.path(out, "isoform_per_pair.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summ <- data.frame(
  pairs_total = nrow(d), pairs_with_isoform_choice = nrow(multi), sampled_with_two_or_more_defined = nrow(per),
  median_ratio_range = stats::median(per$range), p90_ratio_range = as.numeric(stats::quantile(per$range, 0.9)),
  frac_range_over_0.1 = mean(per$range > 0.1), frac_twofold = mean(per$max > 2 * pmax(per$min, 1e-4)),
  frac_crossing_1 = mean(per$crosses_1), reported_above_1 = sum(take$dNdS > 1, na.rm = TRUE)
)
utils::write.table(summ, file.path(out, "isoform_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(t(summ))
