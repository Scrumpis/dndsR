#!/usr/bin/env bash
# Submits stage 2B (docs/dev/phylo-calibration-stage2b.md) from a frozen results root.
#   tools/dev/validation/phylo-stage2b-submit.sh <results root> <K>
#       simulation of every 2B scenario, the essential units except B2, and the summary
#   tools/dev/validation/phylo-stage2b-submit.sh <results root> <K> B2
#       B2 with the start count in <root>/b2_decision.txt ("1start" or "3start"), which must
#       equal the D4 audit result in <root>/summary/audit_decision.txt; then the summary again
# Resubmitting is safe: finished scenarios and units are skipped. Submission is refused while
# any job in <root>/jobs.txt is still queued or running, so a study cannot be submitted twice.
# Full-fit units run on 48 cores; units that fit only the alignments a filter changed (often a
# few genes) run in a separate array on 8 cores, so that cores are not held idle. The thread
# count only sets how many orthogroups run at once; each HyPhy fit uses one CPU.
set -euo pipefail
root=$(realpath "$1"); K=${2:-16}; mode=${3:-main}
[ -f "$root/freeze.json" ] || { echo "freeze first (phylo-stage2-freeze.sh)" >&2; exit 1; }
v=$root/src/tools/dev/validation
mkdir -p "$root/logs"
exec 8>"$root/submit.lock"; flock -n 8 || { echo "another submission of $root is in progress" >&2; exit 1; }
if [ -f "$root/jobs.txt" ]; then
  live=$(squeue -h -u "$USER" -o '%A %F' | tr ' ' '\n' | sort -u) || { echo "squeue failed; not submitting" >&2; exit 1; }
  active=$(grep -o '[0-9]\{6,\}' "$root/jobs.txt" | sort -u | comm -12 - <(printf '%s\n' "$live"))
  [ -z "$active" ] || { echo "jobs from $root/jobs.txt are still queued or running; not submitting again:" $active >&2; exit 1; }
fi
rm -rf "$root/running"   # no job of this root is active, so any claim marker is stale
CHANGED_CPUS=8
sub() { sbatch --parsable --export=ALL "$@"; }
subu() { local units=$1; shift; sbatch --parsable --export=ALL,P3S2_UNITS="$units" "$@"; }
sif=$(python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))['image'])" "$root/freeze.json")
lib=$(python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))['library'])" "$root/freeze.json")
summary() {  # $1: dependency; remaining: unit tables
  local dep=$1; shift
  sub --dependency=afterany:"$dep" --constraint=agg --cpus-per-task=4 --mem=32G --time=6:00:00 --output="$root/logs/summary-%j.out" \
      --wrap "singularity exec --cleanenv --bind $root --env R_LIBS_USER=$lib $sif Rscript --vanilla $v/phylo-stage2b-summary.R $root $root/summary $*"
}
# Ceiling approved by the maintainer for the whole essential workflow, and the share kept for B2.
CEILING=1806; B2_ONE_START=120
projected() {  # core-hours used so far, and a table's estimate scaled by measured/estimated cost
  python3 -I - "$root" "$v/phylo-stage2b-estimates.tsv" "$v/$1" <<'PY'
import csv, glob, os, sys
root, estf, table = sys.argv[1:]
est = {r["unit"]: float(r["core_hours"]) for r in csv.DictReader(open(estf), delimiter="\t")}
used = 0.0; e_done = 0.0
for m in glob.glob(os.path.join(root, "units", "*", "meta.tsv")) + glob.glob(os.path.join(root, "ledger", "*.tsv")):
    r = next(csv.DictReader(open(m), delimiter="\t"))
    used += float(r["seconds"]) * float(r["threads"]) / 3600
    if "/units/" in m: e_done += est.get(r["unit"], 0.0)
ratio = used / e_done if e_done > 0 else 1.0
todo = [r["unit"] for r in csv.DictReader(open(table), delimiter="\t")
        if not os.path.exists(os.path.join(root, "units", r["unit"], "DONE"))]
proj = sum(est.get(u, 0.0) for u in todo) * max(ratio, 1.0)
print("%.1f %.1f %.2f" % (used, proj, ratio))
PY
}
if [ "$mode" = main ]; then
  read -r used proj ratio < <(projected phylo-stage2b-units.tsv)
  cap=$((CEILING - B2_ONE_START))
  echo "used $used, main blocks projected $proj (ratio $ratio), cap for the main blocks $cap"
  awk -v u="$used" -v p="$proj" -v c="$cap" 'BEGIN { exit !(u + p <= c) }' || { echo "projected total exceeds the cap; ask the maintainer" >&2; exit 1; }
  echo "$cap" > "$root/budget_cap.txt"
  ns=$(($(wc -l < "$v/phylo-stage2b-scenarios.tsv") - 1))
  a=$(sub --array=1-"$ns" --output="$root/logs/sim-%A_%a.out" "$v/phylo-stage2-sim.sb" "$root" phylo-stage2b-scenarios.tsv)
  b=$(subu full --dependency=afterok:"$a" --array=0-$((K - 1)) --output="$root/logs/run-%A_%a.out" \
      "$v/phylo-stage2-run.sb" "$root" "$K" phylo-stage2b-units.tsv)
  b2=$(subu changed --dependency=afterok:"$a" --array=0-$((K - 1)) --cpus-per-task=$CHANGED_CPUS --mem=24G \
       --output="$root/logs/runchg-%A_%a.out" "$v/phylo-stage2-run.sb" "$root" "$K" phylo-stage2b-units.tsv)
  c=$(summary "$b:$b2" phylo-stage2b-units.tsv)
  printf 'sim %s\nrun %s\nrunchanged %s\nsummary %s\n' "$a" "$b" "$b2" "$c" | tee -a "$root/jobs.txt"
elif [ "$mode" = B2 ]; then
  missing=$(tail -n +2 "$v/phylo-stage2b-units.tsv" | cut -f1 | while read -r u; do [ -f "$root/units/$u/DONE" ] || echo "$u"; done)
  [ -z "$missing" ] || { echo "main units not finished; B2 waits:" >&2; echo "$missing" >&2; exit 1; }
  [ ! -s "$root/BUDGET_STOP" ] || { echo "$root/BUDGET_STOP is not empty; the maintainer decides" >&2; exit 1; }
  want=$(cat "$root/b2_decision.txt"); got=$(cat "$root/summary/audit_decision.txt")
  [ "$want" = "$got" ] || { echo "b2_decision.txt ($want) differs from the audit ($got)" >&2; exit 1; }
  t=phylo-stage2b-units-B2-$want.tsv
  read -r used proj ratio < <(projected "$t")
  echo "used $used, B2 ($want) projected $proj (measured/estimated $ratio), ceiling $CEILING"
  awk -v u="$used" -v p="$proj" -v c="$CEILING" 'BEGIN { exit !(u + p <= c) }' || { echo "projected total exceeds the ceiling; ask the maintainer" >&2; exit 1; }
  echo "$CEILING" > "$root/budget_cap.txt"; rm -f "$root/BUDGET_STOP"
  b=$(sub --array=0-$((K - 1)) --output="$root/logs/runB2-%A_%a.out" "$v/phylo-stage2-run.sb" "$root" "$K" "$t")
  c=$(summary "$b" phylo-stage2b-units.tsv "$t")
  printf 'runB2 %s (%s)\nsummary %s\n' "$b" "$want" "$c" | tee -a "$root/jobs.txt"
else
  echo "mode: main or B2" >&2; exit 1
fi
