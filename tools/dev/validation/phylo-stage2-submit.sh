#!/usr/bin/env bash
# Submits stage 2A (docs/dev/phylo-calibration-stage2.md) from a frozen results root:
# simulation (one array task per scenario), analyses (K array tasks), and the summary.
# Resubmitting is safe: finished scenarios and units are skipped.
#   tools/dev/validation/phylo-stage2-submit.sh <results root> [K]
set -euo pipefail
root=$(realpath "$1"); K=${2:-12}
[ -f "$root/freeze.json" ] || { echo "freeze first (phylo-stage2-freeze.sh)" >&2; exit 1; }
v=$root/src/tools/dev/validation
mkdir -p "$root/logs"
ns=$(($(wc -l < "$v/phylo-stage2-scenarios.tsv") - 1))
sub() { sbatch --parsable --export=ALL "$@"; }
a=$(sub --array=1-"$ns" --output="$root/logs/sim-%A_%a.out" "$v/phylo-stage2-sim.sb" "$root")
b=$(sub --dependency=afterok:"$a" --array=0-$((K - 1)) --output="$root/logs/run-%A_%a.out" "$v/phylo-stage2-run.sb" "$root" "$K")
sif=$(python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))['image'])" "$root/freeze.json")
lib=$(python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))['library'])" "$root/freeze.json")
c=$(sub --dependency=afterany:"$b" --constraint=agg --cpus-per-task=4 --mem=16G --time=2:00:00 --output="$root/logs/summary-%j.out" \
      --wrap "singularity exec --cleanenv --bind $root --env R_LIBS_USER=$lib $sif Rscript --vanilla $v/phylo-stage2-summary.R $root")
printf 'sim %s\nrun %s\nsummary %s\n' "$a" "$b" "$c" | tee -a "$root/jobs.txt"
