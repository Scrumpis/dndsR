test_that("terminal stops are allowed; internal stops and bad lengths are found", {
  skip_if_not_installed("Biostrings")
  cds <- Biostrings::DNAStringSet(c(ok = "ATGAAATAA", stop_inside = "ATGTAAAAATAA", short = "ATGAAAT"))
  got <- .cds_quality(cds)
  expect_equal(got$multiple_of_3, c(TRUE, TRUE, FALSE))
  expect_equal(got$internal_stop, c(FALSE, TRUE, NA))
})

test_that("the genetic code is respected", {
  skip_if_not_installed("Biostrings")
  # TGA is a stop in the standard code and tryptophan in code 4.
  cds <- Biostrings::DNAStringSet("ATGTGAAAATAA")
  expect_true(.cds_quality(cds, 1L)$internal_stop)
  expect_false(.cds_quality(cds, 4L)$internal_stop)
})
