#!/usr/bin/env Rscript
# Fails when R/ uses pkg::fun from a package DESCRIPTION does not declare.
#
# `R CMD check` reports this as "checking dependencies in R code ... WARNING", and it
# is the failure that reached CI on 2026-11-14 for GenomicRanges and IRanges: the
# local gate (tests, skill, safeguards, tutorial) never runs R CMD check, so nothing
# caught it before the push. This costs a fraction of a second.
#
#   Rscript --vanilla tools/dev/check-imports.R [repo root]
#
# Exit status: 0 when every package used with `::` is declared, 1 otherwise.
root <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
root <- normalizePath(root, mustWork = TRUE)

desc <- read.dcf(file.path(root, "DESCRIPTION"))
declared <- character(0)
for (field in c("Depends", "Imports", "Suggests", "LinkingTo")) {
  if (!field %in% colnames(desc)) next
  entries <- trimws(strsplit(desc[1, field], ",")[[1]])
  entries <- sub("\\s*\\(.*$", "", entries)
  declared <- c(declared, entries[nzchar(entries)])
}
declared <- setdiff(unique(declared), "R")

files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
if (!length(files)) {
  cat("check-imports: no R/ files found\n")
  quit(status = 0)
}

# pkg::fun and pkg:::fun, ignoring comments and string literals (messages such as
# "[dndsR::go_enrichment]" and XPath "./*[not(self::defs)]" are not code).
base_packages <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "tools", "parallel")
used <- character(0)
for (f in files) {
  lines <- readLines(f, warn = FALSE)
  lines <- sub("#.*$", "", lines)
  lines <- gsub('"[^"]*"', '""', lines)
  lines <- gsub("'[^']*'", "''", lines)
  hits <- regmatches(lines, gregexpr("[A-Za-z][A-Za-z0-9._]*:::{0,2}[A-Za-z._][A-Za-z0-9._]*", lines))
  used <- c(used, sub(":::?.*$", "", unlist(hits)))
}
used <- sort(setdiff(unique(used[nzchar(used)]), base_packages))
missing <- setdiff(used, declared)

if (length(missing)) {
  cat("check-imports: used in R/ but not declared in DESCRIPTION:\n")
  for (p in missing) {
    where <- vapply(files, function(f) {
      any(grepl(paste0("\\b", gsub(".", "\\.", p, fixed = TRUE), "::"), readLines(f, warn = FALSE)))
    }, logical(1))
    cat("  ", p, " (", paste(basename(files[where]), collapse = ", "), ")\n", sep = "")
  }
  cat("Add it to Imports (or Suggests) in DESCRIPTION, or use an existing dependency.\n")
  quit(status = 1)
}

cat("check-imports: ", length(used), " packages used with ::, all declared\n", sep = "")
