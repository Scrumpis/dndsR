# The documented workflow with a reversed row: every step runs with --unique-pairs,
# and without it a later step reports the missing output instead of skipping it.
# Runs orthologr on the bundled example, so it is slow; kept in its own file.
skip_if_not_installed("orthologr")
skip_if_not_installed("txdbmaker")
skip_if(Sys.which("diamond") == "", "DIAMOND is not on PATH")
skip_if(fast_tests(), "slow: runs orthologr")

use_example_copy()
writeLines(c("example_AvB spA.fasta spA.gff3 spB.fasta spB.gff3",
             "example_BvA spB.fasta spB.gff3 spA.fasta spA.gff3"), "both.txt")

test_that("extract, calculate, annotate, and summarise all skip the reversed row", {
  expect_message(extract_cds(comparison_file = "both.txt", output_dir = ".", unique_pairs = TRUE),
                 "skipping 1 comparison\\(s\\).*example_BvA")
  suppressWarnings(expect_message(
    calculate_dnds(comparison_file = "both.txt", output_dir = ".", unique_pairs = TRUE, comp_cores = test_cores()),
    "skipping 1 comparison\\(s\\).*example_BvA"))
  expect_false(file.exists("example_BvA/example_BvA_dnds.tsv"))
  out <- suppressMessages(append_annotations(comparison_file = "both.txt", output_dir = ".", unique_pairs = TRUE))
  expect_true(file.exists("example_AvB/example_AvB_dnds_annot.tsv"))
  expect_false(file.exists("example_BvA/example_BvA_dnds_annot.tsv"))
  s <- suppressMessages(dnds_summary(comparison_file = "both.txt", output_dir = ".", regions_bed = "regions.bed",
                                     make_plots = FALSE, unique_pairs = TRUE))
  expect_false(any(grepl("example_BvA", s)))
})

test_that("the kept pair's results equal a run of that row alone", {
  writeLines("example_AvB spA.fasta spA.gff3 spB.fasta spB.gff3", "single.txt")
  dir.create("alone")
  suppressMessages(extract_cds(comparison_file = "single.txt", output_dir = "alone"))
  suppressMessages(suppressWarnings(calculate_dnds(comparison_file = "single.txt", output_dir = "alone",
                                                   comp_cores = test_cores())))
  a <- utils::read.delim("alone/example_AvB/example_AvB_dnds.tsv")
  b <- utils::read.delim("example_AvB/example_AvB_dnds.tsv")
  expect_identical(a, b)
})

test_that("without unique_pairs a later step reports the missing output", {
  expect_error(suppressMessages(append_annotations(comparison_file = "both.txt", output_dir = ".")),
               "dN/dS file not found: .*example_BvA")
})
