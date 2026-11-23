#!/usr/bin/env bash
# Checks the simulated inputs each study actually used (<root>/units/<unit>/sim) against their cells'
# designs with sim-check.R (docs/dev/simulator-integrity.md). Cells come from the units files in each
# root's frozen src. Run in the image.  bash sim-stored-check.sh <runtime dir> <current repo> <out tsv>
set -uo pipefail
rt=$1; repo=$2; outf=$3
printf 'root\tunit\tcell\tchecks_ok\tchecks_failed\n' > "$outf"
for root in study-L study-B2 study-B3 study-Q aln-pilot study-R2 aln-confirm study-Q2; do
  for uf in "$rt/$root"/src/tools/dev/validation/*units*.tsv; do
    awk -F'\t' 'NR > 1 {print $1 "\t" $2}' "$uf"
  done | sort -u | while IFS=$'\t' read -r unit cell; do
    d="$rt/$root/units/$unit/sim"; [ -f "$d/gene_branch_truth.tsv" ] || continue
    res=$(Rscript --vanilla "$repo/tools/dev/validation/sim-check.R" "$d" "$cell" 2>&1)
    printf '%s\t%s\t%s\t%s\t%s\n' "$root" "$unit" "$cell" "$(grep -c '^ok' <<<"$res")" "$(grep -c '^FAIL' <<<"$res")" >> "$outf"
    grep '^FAIL' <<<"$res" | sed "s/^/  $root $unit: /"
  done
done
cat "$outf"
