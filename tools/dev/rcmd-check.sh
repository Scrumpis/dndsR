#!/usr/bin/env bash
# Builds the package and runs R CMD check --no-manual on the tarball in the pinned image,
# as CI does, and fails on an ERROR or WARNING, as CI does. About 95 s, serial.
# Extra arguments go to R CMD check, for example --as-cran before a release.
#   DNDSR_SIF=/path/image.sif DNDSR_DEV_LIB=/path/rlib tools/dev/rcmd-check.sh [--as-cran]
# Also run by tools/dev/check.sh --rcmd-check.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
root=$(git -C "$here" rev-parse --show-toplevel)
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
cd "$work"
export DNDSR_NO_INSTALL=1
"$here/in-container.sh" R CMD build --no-manual --no-build-vignettes "$root"
tarball=$(ls dndsR_*.tar.gz)
"$here/in-container.sh" env _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual "$@" "$tarball" ||
  { tail -40 dndsR.Rcheck/00check.log 2>/dev/null; exit 1; }
status=$(grep '^Status:' dndsR.Rcheck/00check.log)
echo "$status"
case "$status" in *ERROR*|*WARNING*) exit 1 ;; esac
