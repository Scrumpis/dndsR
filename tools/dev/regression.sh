#!/usr/bin/env bash
# Real-genome regression on the small public fixture (tools/dev/fixtures/atal/): runs
# extract_cds, calculate_dnds, and append_annotations through the CLI in the image, and
# compares every output with the committed baseline (expected/outputs.sha256). A
# difference fails the run and prints which pairs and values changed.
#   DNDSR_SIF=... DNDSR_DEV_LIB=... tools/dev/regression.sh [--threads N] [--keep DIR] [--update]
# --update rewrites the baseline from this run; commit it only with a recorded reason
# (docs/dev/TESTING.md). Uses the image's dependencies, so compare versions in one image.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
fix="$here/fixtures/atal"
threads=4; keep=; update=
while [ $# -gt 0 ]; do
  case "$1" in
    --threads) threads=$2; shift ;;
    --keep) keep=$2; shift ;;
    --update) update=1 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac; shift
done
work=${keep:-$(mktemp -d)}
[ -n "$keep" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"; work=$(cd "$work" && pwd)
for f in Athaliana.fasta Athaliana.gff3 Alyrata.fasta Alyrata.gff3; do gzip -dc "$fix/$f.gz" > "$work/$f"; done
echo "AtAl Athaliana.fasta Athaliana.gff3 Alyrata.fasta Alyrata.gff3" > "$work/comparison.txt"

run() { "$here/in-container.sh" bash -c 'cd "$1" && shift && Rscript -e "dndsR::cli_main()" -- "$@"' _ "$work" "$@"; }
s=$SECONDS
run extract_cds -c comparison.txt -o out -t "$threads" >"$work/extract.log" 2>&1
DNDSR_NO_INSTALL=1 run calculate_dnds -c comparison.txt -o out -t "$threads" >"$work/calculate.log" 2>&1
DNDSR_NO_INSTALL=1 run append_annotations -c comparison.txt -o out -t "$threads" >"$work/annotate.log" 2>&1
echo "pipeline: $((SECONDS - s)) s with $threads threads"

# Outputs compared: every file except logs and provenance, which hold times and versions.
(cd "$work/out" && find . -type f ! -name '*.log' ! -name '*provenance.json' ! -name '*.fai' | sort | xargs sha256sum) > "$work/outputs.sha256"
if [ -n "$update" ]; then
  mkdir -p "$fix/expected"
  cp "$work/outputs.sha256" "$fix/expected/outputs.sha256"
  for t in AtAl_dnds.tsv; do cp "$work/out/AtAl/$t" "$fix/expected/$t"; done
  echo "baseline updated: $fix/expected"; exit 0
fi
if diff "$fix/expected/outputs.sha256" "$work/outputs.sha256" >"$work/sha.diff"; then
  echo "regression: PASS, $(wc -l < "$work/outputs.sha256") outputs byte-identical to the baseline"
else
  echo "regression: FAIL, outputs differ from the baseline:"; cat "$work/sha.diff"
  DNDSR_NO_INSTALL=1 "$here/in-container.sh" Rscript --vanilla "$here/fixtures/compare-dnds.R" "$fix/expected/AtAl_dnds.tsv" "$work/out/AtAl/AtAl_dnds.tsv" || true
  echo "Work kept in: $work"; trap - EXIT; exit 1
fi
