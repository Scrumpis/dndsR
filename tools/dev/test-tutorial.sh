#!/usr/bin/env bash
# Runs every command block in the getting-started tutorial, in order, in a clean HOME.
# The clone and image download are replaced with local copies of what they would fetch.
#   DNDSR_SIF=/path/image.sif [DNDSR_ENGINE=apptainer] tools/dev/test-tutorial.sh
# Exits nonzero if any block fails or the tutorial's expected final message is absent.
set -eu
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
: "${DNDSR_SIF:?Set DNDSR_SIF to the runtime image}"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
export HOME=$work/home
mkdir -p "$HOME" "$work/run"

# Test what would be committed, not the working tree's untracked files.
git -C "$root" stash list >/dev/null
git -C "$root" archive --format=tar HEAD | { mkdir -p "$work/source" && tar -x -C "$work/source"; }
git -C "$work/source" init -q
script=$work/tutorial.sh
{
  echo 'set -euo pipefail'
  awk '/^```$/ { inside = !inside; next } inside' "$root/docs/tutorials/getting-started.md"
} |
  sed -e "s#^git clone .*#cp -r '$work/source' dndsR#" \
      -e "s#^singularity pull dndsr.sif .*#cp '$DNDSR_SIF' dndsr.sif#" \
      -e "s#--engine singularity#--engine ${DNDSR_ENGINE:-singularity}#" \
      >"$script"

cd "$work/run"
bash -x "$script" >"$work/output.txt" 2>&1 || {
  tail -30 "$work/output.txt"
  echo "FAILED: a tutorial command failed" >&2
  exit 1
}
grep -F "All 14 example checks passed." "$work/output.txt"
# Facts the tutorial states about the example.
grep -qx "34" "$work/output.txt" || { echo "FAILED: spA no longer has 34 written CDS" >&2; exit 1; }
grep -Fq "example_AvB_q_KEGG_enrichment.tsv" "$work/output.txt" || { echo "FAILED: KEGG output missing" >&2; exit 1; }
grep -Fq "example_AvB_q_IPR_Domain_enrichment.tsv" "$work/output.txt" || { echo "FAILED: IPR output missing" >&2; exit 1; }
echo "Tutorial OK: every command ran in a clean HOME"
