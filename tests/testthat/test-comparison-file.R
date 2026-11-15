# How calculate_dnds reads a comparison file. Each call stops at the first missing CDS
# FASTA, before any alignment, and that error names the comparison and the files it
# looked for, which shows how the file was parsed. aligner = "blast" skips the DIAMOND
# check, so these tests need no external software.

read_as <- function(lines) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  f <- file.path(dir, "comparisons.txt")
  writeLines(lines, f)
  tryCatch(
    calculate_dnds(comparison_file = f, output_dir = dir, aligner = "blast"),
    error = function(e) conditionMessage(e)
  )
}

test_that("a headerless five-column file is read by position, as before", {
  msg <- read_as("c1 a.fa a.gff b.fa b.gff")
  expect_match(msg, "Sequence files not found for c1 ")
  expect_match(msg, "c1/a_CDS.fasta", fixed = TRUE)
  expect_match(msg, "c1/b_CDS.fasta", fixed = TRUE)
})

test_that("a header row is read as a header, not as a comparison", {
  msg <- read_as(c("comparison_name query_fasta query_gff subject_fasta subject_gff",
                   "c1 a.fa a.gff b.fa b.gff"))
  expect_match(msg, "Sequence files not found for c1 ")
  expect_no_match(msg, "comparison_name")
})

test_that("header names, not positions, decide the columns", {
  msg <- read_as(c("comparison_name subject_fasta subject_gff query_fasta query_gff",
                   "c1 b.fa b.gff a.fa a.gff"))
  expect_match(msg, "query_file:   .*c1/a_CDS.fasta")
  expect_match(msg, "subject_file: .*c1/b_CDS.fasta")
})

test_that("an empty line before the header is skipped", {
  msg <- read_as(c("", "comparison_name query_fasta query_gff subject_fasta subject_gff",
                   "c1 a.fa a.gff b.fa b.gff"))
  expect_match(msg, "Sequence files not found for c1 ")
})

test_that("a header with empty MK columns passes the MK check", {
  msg <- read_as(c("comparison_name query_fasta query_gff subject_fasta subject_gff snp_file snp_side",
                   "c1 a.fa a.gff b.fa b.gff NA NA"))
  expect_match(msg, "Sequence files not found for c1 ")
})

test_that("a header with a bad MK column fails before any work, without the headerless hint", {
  msg <- read_as(c("comparison_name query_fasta query_gff subject_fasta subject_gff snp_file snp_side",
                   "c1 a.fa a.gff b.fa b.gff snps.vcf sideways"))
  expect_match(msg, "snp_side must be 'query' or 'subject', not 'sideways' (comparison 'c1')", fixed = TRUE)
  expect_no_match(msg, "has no header")
})

test_that("the headerless MK hint still appears for a headerless file", {
  msg <- read_as(c("c1 a.fa a.gff b.fa b.gff extra", "c2 a.fa a.gff b.fa b.gff extra"))
  expect_match(msg, "snp_side is missing")
  expect_match(msg, "has no header, so column 6 is read as snp_file")
})

test_that("a headered FASTA list is read as a header", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "seqs.txt")
  writeLines(c("comparison_name query_seq subject_seq", "c1 q.fa s.fa"), f)
  msg <- tryCatch(suppressWarnings(calculate_dnds(seq_list_file = f, output_dir = dir, aligner = "blast")),
                  error = function(e) conditionMessage(e))
  expect_no_match(msg, "comparison_name")
  expect_match(msg, "not found")
})
