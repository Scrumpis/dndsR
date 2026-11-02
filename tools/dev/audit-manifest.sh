#!/usr/bin/env bash
# Records and checks content hashes for the local audit record. Timestamps there are
# unreliable, because scratch files are touched to survive the purge.
#   audit-manifest.sh          add newly seen versions to .audit/manifest.sha256
#   audit-manifest.sh --check  report whether every archived file has a recorded hash
# The manifest is append-only: transcripts grow, so one path can have several hashes.
set -eu
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
cd "$root/.audit"
manifest=manifest.sha256
touch "$manifest"

current=$(find . -type f ! -name "$(basename "$manifest")" ! -name events.jsonl ! -name git-events.tsv -print0 |
  sort -z | xargs -0 -r sha256sum)

if [ "${1:-}" = --check ]; then
  missing=$(printf '%s\n' "$current" | grep -vxFf "$manifest" | grep -c . || true)
  echo "$(printf '%s\n' "$current" | grep -c .) files, $missing with unrecorded content"
  [ "$missing" -eq 0 ]
else
  printf '%s\n' "$current" | grep -vxFf "$manifest" >>"$manifest" || true
  echo "$(grep -c . "$manifest") hashes recorded in .audit/$manifest"
fi
