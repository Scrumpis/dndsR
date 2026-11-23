# Custom gene-set family for the enrichment-correction study (docs/dev/enrichment-correction-study.md):
# 500 random gene sets (sizes log-uniform 5-500, so sparse sets are included) and 100 overlapping and
# nested sets, each the union of 2 or 3 InterPro entries' gene sets (entries of 10-300 genes).
#   Rscript enrich-corr-custom-sets.R <dnds_annot.tsv> <seed> <out.tsv>   (gene, set)
a <- commandArgs(TRUE); RNGversion("4.0.0"); set.seed(as.integer(a[2]))
d <- read.table(a[1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d <- d[!duplicated(d$query_id), ]
genes <- d$query_id
sizes <- round(exp(stats::runif(500, log(5), log(500))))
rnd <- lapply(seq_along(sizes), function(i) data.frame(gene = sample(genes, sizes[i]), set = sprintf("RAND%03d", i)))
ipr <- strsplit(ifelse(is.na(d$q_ipr), "", d$q_ipr), "[;,|]")
members <- split(rep(genes, lengths(ipr)), unlist(ipr)); members <- members[nzchar(names(members))]
pool <- names(members)[lengths(members) >= 10 & lengths(members) <= 300]
uni <- lapply(seq_len(100), function(i) { k <- sample(2:3, 1); t <- sample(pool, k)
  data.frame(gene = unique(unlist(members[t])), set = sprintf("UNION%03d", i)) })
out <- do.call(rbind, c(rnd, uni)); out <- out[!duplicated(out), ]
utils::write.table(out, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
