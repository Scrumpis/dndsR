# Results must not depend on how many threads produced them.
skip_if_not_installed("orthologr")
skip_if_not_installed("txdbmaker")
skip_if(Sys.which("diamond") == "", "DIAMOND is not on PATH")
skip_on_cran()
skip_if(fast_tests(), "fast mode")
cores <- min(4L, parallel::detectCores())
skip_if(cores < 2L, "needs two cores")

use_example_copy()

run <- function(threads, cores) {
  folder <- sprintf("t%d", threads)
  dir.create(folder)
  file.copy(list.files(".", "^sp|^comparison_file", full.names = TRUE), folder)
  withr::local_dir(folder)
  suppressMessages(suppressWarnings({
    extract_cds(comparison_file = "comparison_file.txt", output_dir = ".", threads = threads)
    calculate_dnds(comparison_file = "comparison_file.txt", output_dir = ".", comp_cores = cores)
    append_annotations(comparison_file = "comparison_file.txt", output_dir = ".", threads = threads)
    term_enrichment(comparison_file = "comparison_file.txt", output_dir = ".", make_plots = FALSE, verbose = FALSE, threads = threads)
    ipr_enrichment(comparison_file = "comparison_file.txt", output_dir = ".", make_plots = FALSE, threads = threads,
                   entries_source = "local", tree_source = "none",
                   entries_path = system.file("extdata", "interpro_entry.list.tsv", package = "dndsR"))
  }))
  out <- list.files("example_AvB", "\\.(tsv|fasta)$", full.names = TRUE)
  stats::setNames(tools::md5sum(out), basename(out))
}

test_that("one worker and several give identical files", {
  one <- run(1L, 2L)
  several <- run(cores, cores)
  expect_gt(length(one), 3L)
  expect_true(any(grepl("_dnds\\.tsv$", names(one))))
  expect_true(any(grepl("enrichment", names(one))))
  expect_setequal(names(several), names(one))
  expect_identical(several[names(one)], one)
})
