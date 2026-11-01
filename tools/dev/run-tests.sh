#!/usr/bin/env bash
# Installs the working tree into the container library and runs the package tests.
#   tools/dev/run-tests.sh [--fast] [filter]
# --fast leaves out the slowest tests, for use while developing. A filter runs matching files, such as "enrichment".
# Requires DNDSR_SIF and DNDSR_DEV_LIB (see in-container.sh). Exits nonzero on failure.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
fast=
if [ "${1:-}" = --fast ]; then fast=1; shift; fi
cpus=${DNDSR_TEST_CPUS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}
"$here/in-container.sh" env DNDSR_TEST_CPUS="$cpus" DNDSR_TEST_FAST="$fast" Rscript --vanilla \
  "$here/run-package-tests.R" "$here/../../tests/testthat" "$@" 2>&1 | grep -v 'WARNING: While'
exit "${PIPESTATUS[0]}"
