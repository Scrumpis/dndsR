#!/usr/bin/env bash
# Queues P3 as dependent Slurm jobs, to run unattended: real-genome benchmark, simulation
# parameters and unit scale, bounded pilot (AVX-512 nodes), the processor block (AVX2
# node), and the pilot summary. Run from the repository root.
#   DNDSR_SIF=<image> DNDSR_DEV_LIB=<lib> tools/dev/validation/phylo-p3-chain.sh <results root>
set -euo pipefail
out=$(realpath -m "$1"); mkdir -p "$out/logs"
v=tools/dev/validation
tools/dev/in-container.sh true   # install this checkout into DNDSR_DEV_LIB before jobs start
sub() { sbatch --parsable --export=ALL "$@"; }
a=$(sub --output="$out/logs/realgenome-%j.out" $v/phylo-realgenome.sb "$out/p3a")
b=$(sub --dependency=afterok:$a --output="$out/logs/scale-%j.out" $v/phylo-calibration-scale.sb "$out/p3a" "$out/params.json")
c=$(sub --dependency=afterok:$b --output="$out/logs/pilot-%j.out" $v/phylo-calibration-pilot.sb "$out/params.json" "$out/pilot")
d=$(sub --dependency=afterok:$b --constraint=amr --output="$out/logs/pilot-avx2-%j.out" $v/phylo-calibration-pilot.sb "$out/params.json" "$out/pilot" processor)
e=$(sub --dependency=afterany:$c:$d --cpus-per-task=2 --mem=8G --time=1:00:00 --output="$out/logs/summary-%j.out" \
      --wrap "singularity exec --cleanenv --bind $PWD --bind $out $DNDSR_SIF Rscript --vanilla $v/phylo-calibration-summary.R $out/pilot $out/pilot_summary")
printf 'realgenome %s\nscale %s\npilot %s\npilot_avx2 %s\nsummary %s\n' "$a" "$b" "$c" "$d" "$e" | tee "$out/jobs.txt"
