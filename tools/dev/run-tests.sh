#!/usr/bin/env bash
# Installs the working tree into the container library and runs the package tests.
# Requires DNDSR_SIF and DNDSR_DEV_LIB (see in-container.sh). Exits nonzero on failure.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
"$here/in-container.sh" Rscript -e '
  r <- testthat::test_dir(file.path(Sys.getenv("DNDSR_REPO"), "tests/testthat"), package = "dndsR",
    load_package = "installed", reporter = "summary", stop_on_failure = FALSE)
  d <- as.data.frame(r)
  cat(sprintf("\nRESULT tests=%d expectations=%d failed=%d errors=%d skipped=%d\n",
    nrow(d), sum(d$nb), sum(d$failed), sum(d$error), sum(d$skipped)))
  if (sum(d$failed) + sum(d$error) > 0) quit(status = 1, save = "no")' 2>&1 | grep -v 'WARNING: While'
exit "${PIPESTATUS[0]}"
