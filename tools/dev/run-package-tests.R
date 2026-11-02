# Runs the package tests against the installed dndsR, and fails if any test fails.
#   Rscript tools/dev/run-package-tests.R <path to tests/testthat> [filter]
# DNDSR_TEST_CPUS sets the CPUs to use. Test files run side by side, and the tests that
# calculate dN/dS share what is left. DNDSR_TEST_FAST=1 leaves out the slowest tests.
a <- commandArgs(TRUE)
filter <- if (length(a) > 1L) a[2] else NULL
files <- list.files(a[1], "^test-.*[.]R$")
if (!is.null(filter)) files <- grep(filter, sub("^test-(.*)[.]R$", "\\1", files), value = TRUE)
cpus <- suppressWarnings(as.integer(Sys.getenv("DNDSR_TEST_CPUS", as.character(parallel::detectCores()))))
if (is.na(cpus) || cpus < 1L) cpus <- 1L
workers <- max(1L, min(cpus, length(files)))
Sys.setenv(NOT_CRAN = "true", TESTTHAT_PARALLEL = if (workers > 1L) "true" else "false", TESTTHAT_CPUS = workers,
           DNDSR_TEST_CORES = max(2L, min(8L, cpus %/% 3L)))
cat(sprintf("Running %d test files on %d workers, with %d CPUs\n", length(files), workers, cpus))
wall <- system.time(
  r <- testthat::test_dir(a[1], filter = filter, package = "dndsR", load_package = "installed",
                          reporter = "summary", stop_on_failure = FALSE)
)[["elapsed"]]
d <- as.data.frame(r)
cat(sprintf("\nRESULT tests=%d expectations=%d failed=%d errors=%d skipped=%d wall=%.0fs\n",
            nrow(d), sum(d$nb), sum(d$failed), sum(d$error), sum(d$skipped), wall))
if (sum(d$failed) + sum(d$error) > 0L) {
  bad <- d[d$failed > 0L | d$error, c("file", "test")]
  if (nzchar(Sys.getenv("GITHUB_ACTIONS"))) {
    cat(sprintf("::error title=Package tests failed::%s\n", paste(bad$test, collapse = "%0A")))
  }
  quit(status = 1L, save = "no")
}
