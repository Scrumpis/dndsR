# Separate from the other example tests because converting SVG to PNG takes over a minute.
skip_if_not_installed("orthologr")
skip_if_not_installed("txdbmaker")
skip_if(Sys.which("diamond") == "", "DIAMOND is not on PATH")
skip_if_not_installed("RIdeogram")

use_example_copy()
suppressMessages(suppressWarnings({
  extract_cds(comparison_file = "comparison_file.txt", output_dir = ".")
  calculate_dnds(comparison_file = "comparison_file.txt", output_dir = ".", comp_cores = test_cores())
  append_annotations(comparison_file = "comparison_file.txt", output_dir = ".")
}))

test_that("ideograms are written", {
  suppressMessages(dnds_ideogram(comparison_file = "comparison_file.txt", output_dir = ".", sides = "query", make_png = FALSE))
  expect_true(length(list.files("example_AvB", "ideogram\\.svg$")) > 0L)
})

test_that("ideograms are converted to PNG", {
  skip_if(fast_tests(), "fast mode")
  suppressMessages(dnds_ideogram(comparison_file = "comparison_file.txt", output_dir = ".", sides = "query", overwrite = TRUE))
  expect_true(length(list.files("example_AvB", "ideogram\\.png$")) > 0L)
})
