#!/usr/bin/env bash
# Frozen-study reproduction (docs/dev/simulator-integrity.md): for every cell a study used, generates N
# genes with that study's frozen simulator (<root>/src) and with the current one, from the same seed
# and the study's own params.json, and compares every output file byte for byte. Run in the image.
#   bash sim-reproduce.sh <runtime dir> <current repo> <table> <out tsv>
# <table>: root, cell, tree, seed, codons, n (tab-separated, header).
set -euo pipefail
rt=$1; repo=$2; table=$3; outf=$4
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
printf 'root\tcommit\tcell\ttree\tseed\tn\tfiles\tidentical\n' > "$outf"
while IFS=$'\t' read -r root cell tree seed codons n; do
  [ "$root" = root ] && continue
  r="$rt/$root"; c=$(cut -c1-7 "$r/lib/COMMIT")
  (cd "$r/src" && Rscript --vanilla tools/dev/validation/phylo-lineage-sim2.R "$r/params.json" "$cell" "$tree" "$n" "$seed" "$work/frozen" "$codons" > /dev/null 2>&1)
  (cd "$repo" && Rscript --vanilla tools/dev/validation/phylo-lineage-sim2.R "$r/params.json" "$cell" "$tree" "$n" "$seed" "$work/current" "$codons" > /dev/null 2>&1)
  files=$(cd "$work/frozen" && find . -type f | wc -l)
  if diff -r -q "$work/frozen" "$work/current" > /dev/null; then same=TRUE; else same=FALSE; fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$root" "$c" "$cell" "$tree" "$seed" "$n" "$files" "$same" >> "$outf"
  rm -rf "$work/frozen" "$work/current"
done < "$table"
cat "$outf"
