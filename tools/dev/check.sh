#!/usr/bin/env bash
# Runs the development checks side by side, and reports each one's result and time.
#   DNDSR_SIF=/path/image.sif DNDSR_DEV_LIB=/path/rlib tools/dev/check.sh [--fast]
# --fast takes under a minute with several CPUs: fast package tests, skill, and safeguards.
# Without it, also the slow tests and the tutorial, which is checked out from HEAD.
# Exits nonzero if any check fails, and prints the end of its output.
set -u
here=$(cd "$(dirname "$0")" && pwd)
fast=; [ "${1:-}" = --fast ] && fast=--fast
log=$(mktemp -d); trap 'rm -rf "$log"' EXIT
"$here/in-container.sh" true || exit 1
export DNDSR_NO_INSTALL=1

start() { # name, command...
  local name=$1; shift
  ( s=$SECONDS; "$@" >"$log/$name.log" 2>&1; rc=$?; echo "$rc $((SECONDS - s))" >"$log/$name.result" ) &
}
start tests "$here/run-tests.sh" $fast
start skill "$here/in-container.sh" Rscript --vanilla "$here/validate-skill.R" "$here/../.."
start safeguards "$here/test-safeguards.sh"
[ -n "$fast" ] || start tutorial "$here/test-tutorial.sh"
wait

failed=0
for f in "$log"/*.result; do
  name=$(basename "$f" .result); read -r rc secs <"$f"
  if [ "$rc" = 0 ]; then printf 'ok    %-11s %4ss\n' "$name" "$secs"
  else printf 'FAIL  %-11s %4ss\n' "$name" "$secs"; tail -25 "$log/$name.log"; failed=1; fi
done
grep -h '^RESULT' "$log/tests.log" || true
exit "$failed"
