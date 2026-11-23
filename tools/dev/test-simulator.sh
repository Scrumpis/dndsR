#!/usr/bin/env bash
# Regression test of the validation simulator (tools/dev/validation/phylo-lineage-sim2.R), run in the
# image (check.sh does this through in-container.sh). For every cell the studies use:
#  1. generates a small data set with warnings as errors, and R's partial-matching warnings on
#     ($, arguments, attributes), so an ambiguous field access fails instead of silently matching;
#  2. checks it against the cell's design (sim-check.R, expected values derived from the design
#     documents' numbers);
#  3. compares an MD5 of every output file with tools/dev/validation/sim-golden.tsv, so any change to
#     what a seed generates is caught. UPDATE=1 rewrites the golden file (only after the frozen-study
#     reproduction in docs/dev/simulator-integrity.md has been rerun).
#   tools/dev/test-simulator.sh            (from the repository root, inside the image)
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd); cd "$root"
v=tools/dev/validation; golden=$v/sim-golden.tsv
params=$v/sim-test-params.json
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
fail=0; out=()
while IFS=$'\t' read -r cell tree n seed codons; do
  [ "$cell" = cell ] && continue
  d="$work/${cell}_$tree"
  if ! Rscript -e 'options(warn = 2, warnPartialMatchDollar = TRUE, warnPartialMatchArgs = TRUE, warnPartialMatchAttr = TRUE); source("tools/dev/validation/phylo-lineage-sim2.R")' \
      "$params" "$cell" "$tree" "$n" "$seed" "$d" "$codons" > "$work/log" 2>&1; then
    echo "FAIL $cell $tree: the simulator stopped (warnings are errors here)"; cat "$work/log"; fail=1; continue
  fi
  if [ "$cell" != C0 ]; then Rscript --vanilla "$v/sim-check.R" "$d" "$cell" > "$work/check" 2>&1 || { cat "$work/check"; fail=1; }; { grep -c '^ok' "$work/check" || true; } | sed "s/^/ok   $cell $tree: /; s/$/ design checks/"; fi
  sum=$( (cd "$d" && find . -type f | LC_ALL=C sort | xargs md5sum) | md5sum | cut -c1-32)
  out+=("$cell	$tree	$n	$seed	$codons	$sum")
done < "$v/sim-test-cells.tsv"
if [ "${UPDATE:-0}" = 1 ]; then
  { printf 'cell\ttree\tn\tseed\tcodons\tmd5\n'; printf '%s\n' "${out[@]}"; } > "$golden"; echo "wrote $golden"
else
  diff <(printf 'cell\ttree\tn\tseed\tcodons\tmd5\n'; printf '%s\n' "${out[@]}") "$golden" > "$work/diff" && echo "ok   every cell's output equals the golden checksums" ||
    { echo "FAIL output differs from $golden:"; cat "$work/diff"; fail=1; }
fi
exit $fail
