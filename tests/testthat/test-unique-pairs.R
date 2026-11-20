# Tests for the unordered-pair identity and calculate_dnds(unique_pairs = TRUE).
# A multi-species comparison file may list both A-vs-B and B-vs-A; the option runs
# each unordered pair once. Expectations are the file's own names, nothing computed.

pair_id <- function(...) get(".dnds_pair_identity", envir = asNamespace("dndsR"))(list(...))

test_that("the identity of a genome pair does not depend on order", {
  a <- file.path(tempdir(), "spA.fasta")
  b <- file.path(tempdir(), "spB.fasta")
  expect_identical(pair_id(query_fasta = a, subject_fasta = b),
                   pair_id(query_fasta = b, subject_fasta = a))
  c <- file.path(tempdir(), "spC.fasta")
  expect_false(identical(pair_id(query_fasta = a, subject_fasta = b),
                         pair_id(query_fasta = a, subject_fasta = c)))
})

test_that("relative and absolute spellings of the same file compare equal", {
  rel <- "spA.fasta"; abs <- file.path(getwd(), "spA.fasta")
  expect_identical(pair_id(query_fasta = rel, subject_fasta = "spB.fasta"),
                   pair_id(query_fasta = abs, subject_fasta = file.path(getwd(), "spB.fasta")))
})

test_that("a row without a FASTA pair has no identity", {
  expect_true(is.na(pair_id(query_seq = "a", subject_seq = "b")))
  expect_true(is.na(pair_id(query_fasta = "", subject_fasta = "b.fasta")))
})

test_that("unique_pairs skips a reversed duplicate before any dN/dS work", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "pairs.txt")
  # The files do not exist: the run must report the skip, then fail on the first
  # row's missing sequence files, never on the skipped one.
  writeLines(c(
    "A_v_B a.fasta a.gff3 b.fasta b.gff3",
    "B_v_A b.fasta b.gff3 a.fasta a.gff3",
    "A_v_C a.fasta a.gff3 c.fasta c.gff3"
  ), f)
  expect_message(
    err <- tryCatch(dndsR::calculate_dnds(comparison_file = f, output_dir = dir,
                                          unique_pairs = TRUE, aligner = "blast"),
                    error = function(e) conditionMessage(e)),
    "skipping 1 comparison\\(s\\).*B_v_A"
  )
  expect_match(err, "Sequence files not found for A_v_B")
})

test_that("without unique_pairs every row is attempted as written", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "pairs.txt")
  writeLines(c("A_v_B a.fasta a.gff3 b.fasta b.gff3",
               "B_v_A b.fasta b.gff3 a.fasta a.gff3"), f)
  expect_error(dndsR::calculate_dnds(comparison_file = f, output_dir = dir,
                                     aligner = "blast"),
               "Sequence files not found for A_v_B")
})

# ---- The skip rule shared by every batch command (.unique_pairs_skip) ----
# Expectations are the row names the rule must keep or skip, written by hand.

skipped <- function(lines, header = NULL) {
  f <- withr::local_tempfile(fileext = ".txt", .local_envir = parent.frame())
  writeLines(c(header, lines), f)
  suppressMessages(get(".unique_pairs_skip", envir = asNamespace("dndsR"))(f, TRUE))
}

test_that("reversed and repeated rows without MK inputs are skipped", {
  expect_identical(skipped(c("A_v_B a.fa a.gff b.fa b.gff",
                             "B_v_A b.fa b.gff a.fa a.gff",
                             "A_v_B2 a.fa a.gff b.fa b.gff",
                             "A_v_C a.fa a.gff c.fa c.gff")),
                   c("B_v_A", "A_v_B2"))
})

test_that("nothing is skipped unless asked", {
  f <- withr::local_tempfile(fileext = ".txt")
  writeLines(c("A_v_B a.fa a.gff b.fa b.gff", "B_v_A b.fa b.gff a.fa a.gff"), f)
  expect_identical(get(".unique_pairs_skip", envir = asNamespace("dndsR"))(f, FALSE), character(0))
})

test_that("a reversed row with its own MK inputs is kept, and says why", {
  hdr <- "comparison_name query_fasta query_gff subject_fasta subject_gff snp_file snp_side"
  lines <- c("A_v_B a.fa a.gff b.fa b.gff a.vcf query",
             "B_v_A b.fa b.gff a.fa a.gff b.vcf query")
  f <- withr::local_tempfile(fileext = ".txt"); writeLines(c(hdr, lines), f)
  expect_message(get(".unique_pairs_skip", envir = asNamespace("dndsR"))(f, TRUE, "test"),
                 "keeping B_v_A.*own McDonald-Kreitman inputs")
  expect_identical(skipped(lines, hdr), character(0))
})

test_that("the same SNP file on the same genome, from the other side, is a duplicate", {
  hdr <- "comparison_name query_fasta query_gff subject_fasta subject_gff snp_file snp_side"
  # a.vcf is on genome a in both rows: query of A_v_B, subject of B_v_A.
  expect_identical(skipped(c("A_v_B a.fa a.gff b.fa b.gff a.vcf query",
                             "B_v_A b.fa b.gff a.fa a.gff a.vcf subject"), hdr), "B_v_A")
  # The same file named on the other genome is a different MK input, so it is kept.
  expect_identical(skipped(c("A_v_B a.fa a.gff b.fa b.gff a.vcf query",
                             "B_v_A b.fa b.gff a.fa a.gff a.vcf query"), hdr), character(0))
})

test_that("a later row with MK inputs is kept after an earlier row without them", {
  hdr <- "comparison_name query_fasta query_gff subject_fasta subject_gff snp_file snp_side"
  expect_identical(skipped(c("A_v_B a.fa a.gff b.fa b.gff NA NA",
                             "B_v_A b.fa b.gff a.fa a.gff b.vcf query",
                             "B_v_A2 b.fa b.gff a.fa a.gff b.vcf query"), hdr), "B_v_A2")
})

test_that("headerless files with MK columns follow the same rule", {
  expect_identical(skipped(c("A_v_B a.fa a.gff b.fa b.gff a.vcf query",
                             "B_v_A b.fa b.gff a.fa a.gff b.vcf query",
                             "B_v_A2 b.fa b.gff a.fa a.gff NA NA")), "B_v_A2")
})

test_that("calculate_dnds keeps a reversed row that carries its own MK inputs", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "pairs.txt")
  writeLines(c("A_v_B a.fasta a.gff3 b.fasta b.gff3 a.vcf query",
               "B_v_A b.fasta b.gff3 a.fasta a.gff3 b.vcf query"), f)
  expect_message(
    tryCatch(dndsR::calculate_dnds(comparison_file = f, output_dir = dir, unique_pairs = TRUE, aligner = "blast"),
             error = function(e) NULL),
    "keeping B_v_A")
})

test_that("every batch command that reads a comparison file accepts unique_pairs", {
  for (fn in c("extract_cds", "calculate_dnds", "append_annotations", "dnds_summary", "dnds_contrast",
               "dnds_state_contrast", "go_enrichment", "term_enrichment", "ipr_enrichment",
               "covariate_enrichment", "dnds_ideogram")) {
    f <- formals(getExportedValue("dndsR", fn))
    expect_true("unique_pairs" %in% names(f), info = fn)
    expect_false(isTRUE(f$unique_pairs), info = fn)
  }
})
