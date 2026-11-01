# Set by tools/dev/run-package-tests.R. The defaults suit R CMD check, which allows two cores.
test_cores <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("DNDSR_TEST_CORES", "2")))
  if (is.na(n) || n < 1L) 2L else n
}

# Fast mode leaves out the slowest tests, for use while developing. Run all of them before committing.
fast_tests <- function() nzchar(Sys.getenv("DNDSR_TEST_FAST"))

# Works in a private copy of the bundled example until the test file ends.
use_example_copy <- function() {
  example <- system.file("extdata", "example", package = "dndsR")
  work <- withr::local_tempdir(.local_envir = teardown_env())
  file.copy(list.files(example, full.names = TRUE), work, recursive = TRUE)
  old <- setwd(work)
  withr::defer(setwd(old), teardown_env())
  invisible(work)
}
