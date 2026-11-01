# Runs the package tests against the installed dndsR, and fails if any test fails.
#   Rscript tools/dev/run-package-tests.R <path to tests/testthat>
r <- testthat::test_dir(commandArgs(TRUE)[1], package = "dndsR", load_package = "installed",
                        reporter = "summary", stop_on_failure = FALSE)
d <- as.data.frame(r)
cat(sprintf("\nRESULT tests=%d expectations=%d failed=%d errors=%d skipped=%d\n",
            nrow(d), sum(d$nb), sum(d$failed), sum(d$error), sum(d$skipped)))
if (sum(d$failed) + sum(d$error) > 0L) {
  bad <- d[d$failed > 0L | d$error, c("file", "test")]
  if (nzchar(Sys.getenv("GITHUB_ACTIONS"))) {
    cat(sprintf("::error title=Package tests failed::%s\n", paste(bad$test, collapse = "%0A")))
  }
  quit(status = 1L, save = "no")
}
