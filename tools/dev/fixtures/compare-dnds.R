# Summarises how two dN/dS tables differ: pairs lost or new, and per column the number of
# changed values and the largest absolute change. Base R, run in the image.
#   Rscript --vanilla compare-dnds.R <expected _dnds.tsv> <observed _dnds.tsv>
a <- commandArgs(TRUE)
read <- function(p) { d <- utils::read.delim(p, stringsAsFactors = FALSE); rownames(d) <- paste(d$query_id, d$subject_id); d }
e <- read(a[1]); o <- read(a[2])
lost <- setdiff(rownames(e), rownames(o)); new <- setdiff(rownames(o), rownames(e))
cat(sprintf("pairs: expected %d, observed %d, lost %d, new %d\n", nrow(e), nrow(o), length(lost), length(new)))
for (k in utils::head(lost, 10)) cat("  lost", k, "\n")
for (k in utils::head(new, 10)) cat("  new ", k, "\n")
shared <- intersect(rownames(e), rownames(o))
for (col in c("dN", "dS", "dNdS")) {
  x <- e[shared, col]; y <- o[shared, col]
  changed <- xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
  worst <- suppressWarnings(max(abs(x - y), na.rm = TRUE))
  cat(sprintf("%s: %d changed among shared pairs, largest absolute change %g\n", col, sum(changed), if (is.finite(worst)) worst else 0))
}
