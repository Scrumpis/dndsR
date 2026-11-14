#!/usr/bin/env bash
# Submits stage 2B (docs/dev/phylo-calibration-stage2b.md) from a frozen results root.
#   tools/dev/validation/phylo-stage2b-submit.sh <results root> <K>
#       simulation of every 2B scenario, the essential units except B2, and the summary
#   tools/dev/validation/phylo-stage2b-submit.sh <results root> <K> B2
#       B2 with the start count in <root>/b2_decision.txt ("1start" or "3start"), which must
#       equal the D4 audit result in <root>/summary/audit_decision.txt; then the summary again
# Resubmitting is safe: finished scenarios and units are skipped.
set -euo pipefail
root=$(realpath "$1"); K=${2:-16}; mode=${3:-main}
[ -f "$root/freeze.json" ] || { echo "freeze first (phylo-stage2-freeze.sh)" >&2; exit 1; }
v=$root/src/tools/dev/validation
mkdir -p "$root/logs"
sub() { sbatch --parsable --export=ALL "$@"; }
sif=$(python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))['image'])" "$root/freeze.json")
lib=$(python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))['library'])" "$root/freeze.json")
summary() {  # $1: dependency; remaining: unit tables
  local dep=$1; shift
  sub --dependency=afterany:"$dep" --constraint=agg --cpus-per-task=4 --mem=32G --time=6:00:00 --output="$root/logs/summary-%j.out" \
      --wrap "singularity exec --cleanenv --bind $root --env R_LIBS_USER=$lib $sif Rscript --vanilla $v/phylo-stage2b-summary.R $root $root/summary $*"
}
if [ "$mode" = main ]; then
  ns=$(($(wc -l < "$v/phylo-stage2b-scenarios.tsv") - 1))
  a=$(sub --array=1-"$ns" --output="$root/logs/sim-%A_%a.out" "$v/phylo-stage2-sim.sb" "$root" phylo-stage2b-scenarios.tsv)
  b=$(sub --dependency=afterok:"$a" --array=0-$((K - 1)) --output="$root/logs/run-%A_%a.out" "$v/phylo-stage2-run.sb" "$root" "$K" phylo-stage2b-units.tsv)
  c=$(summary "$b" phylo-stage2b-units.tsv)
  printf 'sim %s\nrun %s\nsummary %s\n' "$a" "$b" "$c" | tee -a "$root/jobs.txt"
elif [ "$mode" = B2 ]; then
  want=$(cat "$root/b2_decision.txt"); got=$(cat "$root/summary/audit_decision.txt")
  [ "$want" = "$got" ] || { echo "b2_decision.txt ($want) differs from the audit ($got)" >&2; exit 1; }
  t=phylo-stage2b-units-B2-$want.tsv
  b=$(sub --array=0-$((K - 1)) --output="$root/logs/runB2-%A_%a.out" "$v/phylo-stage2-run.sb" "$root" "$K" "$t")
  c=$(summary "$b" phylo-stage2b-units.tsv "$t")
  printf 'runB2 %s (%s)\nsummary %s\n' "$b" "$want" "$c" | tee -a "$root/jobs.txt"
else
  echo "mode: main or B2" >&2; exit 1
fi
