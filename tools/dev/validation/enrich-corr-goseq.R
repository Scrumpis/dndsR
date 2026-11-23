# GOseq comparator for enrich-corr-methods.R (docs/dev/enrichment-correction-study.md). Runs in the
# separate comparators environment (bioconductor-goseq 1.58.0, r-biasedurn 2.0.12), not a dndsR dependency.
#   Rscript enrich-corr-goseq.R <dir with genes.tsv (gene, y, len, pfit) and g2c.tsv (gene, term)>
# Writes out.tsv: term, p_len (PWF on length, GOseq's documented use) and p_both (PWF on the fitted
# length-and-dS probability, the adaptation to dN/dS-based selection), Wallenius over-representation p.
suppressMessages(library(goseq))
dir <- commandArgs(TRUE)[1]
g <- utils::read.delim(file.path(dir, "genes.tsv")); m <- utils::read.delim(file.path(dir, "g2c.tsv"), stringsAsFactors = FALSE)
de <- stats::setNames(g$y, as.character(g$gene)); g2c <- split(m$term, as.character(m$gene))
run <- function(bias) {
  pwf <- suppressMessages(suppressWarnings(goseq::nullp(de, bias.data = bias, plot.fit = FALSE)))
  r <- suppressMessages(goseq::goseq(pwf, gene2cat = g2c, method = "Wallenius", use_genes_without_cat = FALSE))
  stats::setNames(r$over_represented_pvalue, r$category)
}
# A variant whose weighting function cannot be fitted is reported as missing, not as a p-value.
safe <- function(b) tryCatch(run(b), error = function(e) stats::setNames(numeric(0), character(0)))
pl <- safe(g$len); pb <- safe(g$pfit)
terms <- unique(m$term)
utils::write.table(data.frame(term = terms, p_len = unname(pl[terms]), p_both = unname(pb[terms])), file.path(dir, "out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
