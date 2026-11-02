# On real gene pairs: compares dndsR's reported values with PAML yn00 and KaKs_Calculator,
# on one independently built codon alignment, and with Comeron on that same alignment.
#   Rscript real-estimators.R <work dir> <tools bin> <pairs per comparison>
# Writes <work>/summary/real_estimators*.tsv. Sampling is stratified so that genes with
# dN/dS above 1 are included. Alignment: global protein alignment, back-translated.
suppressPackageStartupMessages({ library(Biostrings); library(pwalign) })
args <- commandArgs(trailingOnly = TRUE)
work <- args[1]; tools_dir <- args[2]; n_pairs <- as.integer(args[3])
RNGversion("4.0.0"); set.seed(99)
comparisons <- read.table(file.path(work, "comparison.txt"), stringsAsFactors = FALSE,
                          col.names = c("name", "qfa", "qgff", "sfa", "sgff"))
stem <- function(p) tools::file_path_sans_ext(basename(p))
tmp <- tempfile("est"); dir.create(tmp)

backtranslate <- function(aligned_aa, cds) {
  codons <- substring(cds, seq(1, nchar(cds) - 2, 3), seq(3, nchar(cds), 3))
  aa <- strsplit(aligned_aa, "")[[1]]
  out <- character(length(aa)); j <- 1L
  for (k in seq_along(aa)) {
    if (aa[k] == "-") out[k] <- "---" else { out[k] <- codons[j]; j <- j + 1L }
  }
  out
}

align_pair <- function(q, s) {
  # Remove terminal stops, then align proteins.
  drop_stop <- function(x) { p <- as.character(translate(DNAString(x), if.fuzzy.codon = "X")); if (endsWith(p, "*")) substr(x, 1, nchar(x) - 3) else x }
  q <- drop_stop(q); s <- drop_stop(s)
  pq <- as.character(translate(DNAString(q), if.fuzzy.codon = "X"))
  ps <- as.character(translate(DNAString(s), if.fuzzy.codon = "X"))
  a <- pwalign::pairwiseAlignment(AAString(pq), AAString(ps), substitutionMatrix = "BLOSUM62", type = "global")
  cq <- backtranslate(as.character(pwalign::alignedPattern(a)), q)
  cs <- backtranslate(as.character(pwalign::alignedSubject(a)), s)
  keep <- cq != "---" & cs != "---" & !grepl("N", cq) & !grepl("N", cs)
  list(q = cq, s = cs, ungapped_q = paste(cq[keep], collapse = ""), ungapped_s = paste(cs[keep], collapse = ""),
       gapped_q = paste(cq, collapse = ""), gapped_s = paste(cs, collapse = ""))
}

comeron <- function(al) {
  f <- file.path(tmp, "pair.fa")
  writeLines(c(">q", al$gapped_q, ">s", al$gapped_s), f)
  r <- orthologr::substitutionrate(f, est.method = "Comeron", quiet = TRUE)
  c(comeron_dN = r$dN, comeron_dS = r$dS)
}
yn00 <- function(al) {
  old <- setwd(tmp); on.exit(setwd(old))
  writeLines(c(sprintf(" 2 %d", nchar(al$ungapped_q)), "q", al$ungapped_q, "s", al$ungapped_s), "pair.phy")
  writeLines(c("seqfile = pair.phy", "outfile = yn.out", "verbose = 0", "icode = 0", "weighting = 0", "commonf3x4 = 0"), "yn00.ctl")
  system2(file.path(tools_dir, "yn00"), "yn00.ctl", stdout = FALSE, stderr = FALSE)
  res <- readLines("yn.out")
  row <- grep("^\\s*2\\s+1\\s+\\d", res, value = TRUE)[1]
  v <- as.numeric(strsplit(trimws(row), "\\s+")[[1]])
  ng <- grep("^s\\s", grep("\\(", res, value = TRUE), value = TRUE)[1]
  w <- as.numeric(regmatches(ng, gregexpr("-?[0-9]+\\.[0-9]+", ng))[[1]])
  c(yn00_dN = v[8], yn00_dS = v[11], ng86_dN = w[2], ng86_dS = w[3])
}
kaks <- function(al) {
  old <- setwd(tmp); on.exit(setwd(old))
  writeLines(c("q-s", al$ungapped_q, al$ungapped_s, ""), "pair.axt")
  system2(file.path(tools_dir, "KaKs_Calculator"), c("-i", "pair.axt", "-o", "pair.kaks", "-m", "MLWL", "-m", "MYN"), stdout = FALSE, stderr = FALSE)
  r <- utils::read.delim("pair.kaks")
  out <- as.numeric(c(rbind(r$Ka, r$Ks)))
  names(out) <- c(rbind(paste0(r$Method, "_dN"), paste0(r$Method, "_dS")))
  out
}

rows <- list()
for (i in seq_len(nrow(comparisons))) {
  cmp <- comparisons[i, ]
  d <- read.delim(file.path(work, "cur", cmp$name, paste0(cmp$name, "_dnds.tsv")))[, c("query_id", "subject_id", "dN", "dS", "dNdS")]
  q <- as.character(readDNAStringSet(file.path(work, "cur", cmp$name, paste0(stem(cmp$qfa), "_CDS.fasta"))))
  s <- as.character(readDNAStringSet(file.path(work, "cur", cmp$name, paste0(stem(cmp$sfa), "_CDS.fasta"))))
  d <- d[!is.na(d$dNdS), ]
  high <- d[d$dNdS > 1, ]
  low <- d[d$dNdS <= 1, ]
  take <- rbind(high[sample.int(nrow(high), min(nrow(high), n_pairs %/% 4L)), ],
                low[sample.int(nrow(low), n_pairs - min(nrow(high), n_pairs %/% 4L)), ])
  for (k in seq_len(nrow(take))) {
    r <- take[k, ]
    res <- tryCatch({
      al <- align_pair(q[[r$query_id]], s[[r$subject_id]])
      c(comeron(al), yn00(al), kaks(al), aligned_codons = nchar(al$ungapped_q) / 3)
    }, error = function(e) NULL)
    if (!is.null(res)) rows[[length(rows) + 1L]] <- data.frame(comparison = cmp$name, r[, c("dN", "dS", "dNdS")], t(res))
  }
}
x <- do.call(rbind, rows)
utils::write.table(x, file.path(work, "summary", "real_estimators_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ratio <- function(a, b) ifelse(b > 0, a / b, NA)
methods <- list(comeron_same_alignment = c("comeron_dN", "comeron_dS"), yn00 = c("yn00_dN", "yn00_dS"),
                ng86 = c("ng86_dN", "ng86_dS"), MLWL = c("MLWL_dN", "MLWL_dS"), MYN = c("MYN_dN", "MYN_dS"))
summary_rows <- list()
for (cmp in unique(x$comparison)) for (m in names(methods)) {
  d <- x[x$comparison == cmp, ]
  dn <- d[[methods[[m]][1]]]; ds <- d[[methods[[m]][2]]]
  r <- ratio(dn, ds); base <- d$dNdS
  ok <- is.finite(r) & is.finite(base) & r < 10 & base < 10
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    comparison = cmp, method = m, pairs = sum(ok),
    spearman_ratio = stats::cor(r[ok], base[ok], method = "spearman"),
    spearman_dN = stats::cor(dn[ok], d$dN[ok], method = "spearman"),
    spearman_dS = stats::cor(ds[ok], d$dS[ok], method = "spearman"),
    median_ratio_method = stats::median(r[ok]), median_ratio_dndsR = stats::median(base[ok]),
    median_ds_method = stats::median(ds[ok]), median_ds_dndsR = stats::median(d$dS[ok]),
    agree_side_of_1 = mean((r[ok] > 1) == (base[ok] > 1))
  )
}
utils::write.table(do.call(rbind, summary_rows), file.path(work, "summary", "real_estimators.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Done")
