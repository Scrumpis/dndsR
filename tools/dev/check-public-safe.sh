#!/bin/sh
# Refuses content unfit for a public repository: tokens, private paths, transcripts,
# containers, and bulky files. Usage: check-public-safe.sh --staged | --all
# Prints only file names and pattern categories, so it never echoes a matched secret.
set -u
MAX_KB=1024
# Files already in the v0.2.0 release that may exceed the size limit.
LARGE_OK='^inst/extdata/interpro_entry\.list\.tsv$'
# Fragments are split so this script does not match itself.
SECRET='(gh[pousr]_[A-Za-z0-9]{20}|github_pat_[A-Za-z0-9_]{20}|sk-ant-[a-z]|AKIA[0-9A-Z]{16}|-----BEGIN [A-Z ]*PRIVATE'' KEY)'
PRIVATE_PATH='(/mnt/(gs''[0-9]+|ffs''[0-9]+|home)/|/home/[a-z0-9_]+/|/Users/[A-Za-z0-9_]+/)'
BAD_NAME='(^|/)(docs/validation/real-genomes/|\.audit|\.git-credentials|\.env|id_rsa|id_ed25519)(/|$)|\.(jsonl|sif|sqfs|pem|key|bam|cram|fastq|fq)(\.gz)?$'

case "${1:---staged}" in
  --staged) files=$(git diff --cached --name-only --diff-filter=ACMR) ; show() { git show ":$1"; } ;;
  --all) files=$(git ls-files) ; show() { cat "$1"; } ;;
  *) echo "usage: $0 --staged|--all" >&2; exit 2 ;;
esac

fail=0
report() { echo "public-safety: $1: $2" >&2; fail=1; }
IFS='
'
for f in $files; do
  printf '%s\n' "$f" | grep -Eq "$BAD_NAME" && report "$f" "file type or name must not be committed"
  [ -f "$f" ] || [ "$1" = --staged ] || continue
  kb=$(($(show "$f" | wc -c) / 1024))
  if [ "$kb" -gt "$MAX_KB" ] && ! printf '%s\n' "$f" | grep -Eq "$LARGE_OK"; then
    report "$f" "${kb} KB exceeds ${MAX_KB} KB limit"
  fi
  show "$f" | grep -EqI "$SECRET" && report "$f" "looks like a credential"
  show "$f" | grep -EqI "$PRIVATE_PATH" && report "$f" "contains a private machine path"
done
exit "$fail"
