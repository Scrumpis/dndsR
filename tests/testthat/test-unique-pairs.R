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
