# Writes each GO term's ancestors from GO.db, including itself, for enrichment-verify.py go.
#   Rscript enrichment-go-ancestors.R <out.tsv>
suppressMessages(library(GO.db))
out <- commandArgs(TRUE)[1]
maps <- list(BP = GOBPANCESTOR, MF = GOMFANCESTOR, CC = GOCCANCESTOR)
rows <- lapply(names(maps), function(o) {
  a <- as.list(maps[[o]])
  a <- lapply(names(a), function(t) data.frame(ontology = o, term = t, ancestor = unique(c(t, setdiff(a[[t]], "all")))))
  do.call(rbind, a)
})
write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE, row.names = FALSE)
