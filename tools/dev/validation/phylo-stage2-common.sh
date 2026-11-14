# Shared by the stage 2 jobs (docs/dev/phylo-calibration-stage2.md). Sourced, not run.
# Needs $root, the results root written by phylo-stage2-freeze.sh.
set -euo pipefail
freeze=$root/freeze.json
[ -f "$freeze" ] || { echo "no $freeze: run phylo-stage2-freeze.sh first" >&2; exit 1; }
jget() { python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])" "$freeze" "$1"; }
commit=$(jget commit); sif=$(jget image); sif_sha=$(jget image_sha256); lib=$(jget library); src=$(jget source)
# The image and the library must be exactly what was frozen.
[ "$(sha256sum "$sif" | cut -d' ' -f1)" = "$sif_sha" ] || { echo "image digest differs from freeze.json" >&2; exit 1; }
[ ! -w "$lib" ] && [ ! -w "$lib/dndsR" ] || { echo "library $lib is writable; it must be read-only" >&2; exit 1; }
v=$src/tools/dev/validation
threads=${SLURM_CPUS_PER_TASK:-8}
run() { singularity exec --cleanenv --bind "$root" --bind "${TMPDIR:-/tmp}" --pwd "$PWD" \
          --env "R_LIBS_USER=$lib" "$sif" "$@"; }
dndsr() { run Rscript --vanilla -e 'dndsR::cli_main()' -- "$@"; }
cpu_model() { grep -m1 "model name" /proc/cpuinfo | sed 's/.*: //'; }
# A unit's or scenario's specification hash: its table row and the frozen commit.
spec_hash() { printf '%s\n%s\n' "$1" "$commit" | md5sum | cut -d' ' -f1; }
done_ok() { [ -f "$1/DONE" ] && [ "$(cat "$1/DONE")" = "$2" ]; }
