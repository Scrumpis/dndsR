#!/usr/bin/env bash
# Runs a command in the pinned dndsR runtime image, with this clone and a separate
# R library mounted. Installs the working tree first unless DNDSR_NO_INSTALL=1.
#   DNDSR_SIF=/path/image.sif DNDSR_DEV_LIB=/path/rlib tools/dev/in-container.sh Rscript x.R
set -eu
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
: "${DNDSR_SIF:?Set DNDSR_SIF to the runtime image (see docs/dev/RUNBOOK.md)}"
: "${DNDSR_DEV_LIB:?Set DNDSR_DEV_LIB to a writable R library outside the clone}"
mkdir -p "$DNDSR_DEV_LIB"
engine=$(command -v apptainer || command -v singularity) || { echo "Needs Apptainer or Singularity" >&2; exit 1; }
run() {
  "$engine" exec --cleanenv --bind "$root" --bind "$DNDSR_DEV_LIB" --bind "$PWD" \
    --env "R_LIBS_USER=$DNDSR_DEV_LIB" --env "DNDSR_REPO=$root" "$DNDSR_SIF" "$@"
}
if [ "${DNDSR_NO_INSTALL:-0}" != 1 ]; then
  run R CMD INSTALL --no-docs --no-test-load --library="$DNDSR_DEV_LIB" "$root" >/dev/null 2>"$DNDSR_DEV_LIB/.install.log" ||
    { cat "$DNDSR_DEV_LIB/.install.log" >&2; exit 1; }
fi
run "$@"
