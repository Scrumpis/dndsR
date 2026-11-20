#!/usr/bin/env bash
# Runs the development checks side by side, and reports each one's result and time.
#   DNDSR_SIF=/path/image.sif DNDSR_DEV_LIB=/path/rlib tools/dev/check.sh [--fast] [--serial]
# --fast takes under a minute with several CPUs: fast package tests, skill, and safeguards.
# Without it, also the slow tests and the tutorial, which is checked out from HEAD.
# --serial runs the checks one at a time, for limited memory; combine with
# DNDSR_TEST_CPUS=4 to cap the test workers too.
# --rcmd-check adds R CMD check on the built tarball (tools/dev/rcmd-check.sh, about 95 s),
# which catches what only check sees, such as a namespace note or a documentation mismatch.
# --regression adds the real-genome fixture regression (tools/dev/regression.sh, about 2-4 min):
# the full pipeline on ~540 public Arabidopsis genes, compared with the committed baseline.
# Always includes an import check: every pkg::fun in R/ must be declared in DESCRIPTION,
# which is what R CMD check reports and what CI caught on 2026-11-14.
# Exits nonzero if any check fails, and prints the end of its output.
set -u
here=$(cd "$(dirname "$0")" && pwd)
fast=; serial=; rcmd=; regress=
for arg in "$@"; do
  case "$arg" in
    --fast) fast=--fast ;;
    --serial) serial=1 ;;
    --rcmd-check) rcmd=1 ;;
    --regression) regress=1 ;;
    *) echo "Unknown option: $arg" >&2; exit 2 ;;
  esac
done
log=$(mktemp -d); trap 'rm -rf "$log"' EXIT
"$here/in-container.sh" true || exit 1
export DNDSR_NO_INSTALL=1

start() { # name, command...
  local name=$1; shift
  ( s=$SECONDS; "$@" >"$log/$name.log" 2>&1; rc=$?; echo "$rc $((SECONDS - s))" >"$log/$name.result" ) &
}
check_tests() { start tests "$here/run-tests.sh" $fast; }
check_imports() { start imports "$here/in-container.sh" Rscript --vanilla "$here/check-imports.R" "$here/../.."; }
check_skill() { start skill "$here/in-container.sh" Rscript --vanilla "$here/validate-skill.R" "$here/../.."; }
check_safeguards() { start safeguards "$here/test-safeguards.sh"; }
check_tutorial() { start tutorial "$here/test-tutorial.sh"; }
check_rcmd() { start rcmdcheck "$here/rcmd-check.sh"; }
check_regression() { start regression "$here/regression.sh" --threads "${DNDSR_TEST_CPUS:-${SLURM_CPUS_PER_TASK:-4}}"; }

if [ -n "$serial" ]; then
  # One check at a time: the low-memory route (limited RAM, WSL, a small laptop).
  check_imports; wait
  check_tests; wait
  check_skill; wait
  check_safeguards; wait
  [ -n "$fast" ] || { check_tutorial; wait; }
  [ -z "$rcmd" ] || { check_rcmd; wait; }
  [ -z "$regress" ] || { check_regression; wait; }
else
  check_imports; check_tests; check_skill; check_safeguards
  [ -n "$fast" ] || check_tutorial
  [ -z "$rcmd" ] || check_rcmd
  [ -z "$regress" ] || check_regression
  wait
fi

failed=0
for f in "$log"/*.result; do
  name=$(basename "$f" .result); read -r rc secs <"$f"
  if [ "$rc" = 0 ]; then printf 'ok    %-11s %4ss\n' "$name" "$secs"
  else printf 'FAIL  %-11s %4ss\n' "$name" "$secs"; tail -25 "$log/$name.log"; failed=1; fi
done
grep -hE '^(SKIP|RESULT) ' "$log/tests.log" || true
[ -z "$regress" ] || grep -hE '^(pipeline|regression):' "$log/regression.log" || true
exit "$failed"
