#!/usr/bin/env bash
# Setup for each clone: installs the git hooks and creates the audit directory. Safe to repeat.
#   tools/dev/setup.sh [--audit-dir PATH]
# --audit-dir keeps the audit record outside the clone (for example, off a purged
# scratch filesystem) and links it from .audit. Requires a clone of Scrumpis/dndsR.
set -eu
cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
audit_target=
case "${1:-}" in
  --audit-dir) audit_target=${2:?--audit-dir needs a path} ;;
  "") ;;
  *) echo "usage: $0 [--audit-dir PATH]" >&2; exit 2 ;;
esac

git remote get-url origin | grep -Eqi 'github\.com[:/]Scrumpis/dndsR(\.git)?$' ||
  { echo "origin is not Scrumpis/dndsR; refusing to configure" >&2; exit 1; }

# Hooks are copied, not referenced: tools/dev does not exist when main is checked out.
current=$(git config --local --get core.hooksPath || true)
if [ "$current" = tools/dev/githooks ]; then
  git config --local --unset core.hooksPath
elif [ -n "$current" ]; then
  echo "core.hooksPath is already '$current'. Not overriding it; merge the hooks by hand." >&2
  exit 1
fi
hooks=$(git rev-parse --git-common-dir)/hooks
mkdir -p "$hooks"
# Preserve hooks installed by anything else.
for source in tools/dev/githooks/*; do
  target=$hooks/$(basename "$source")
  if [ -e "$target" ] && ! grep -q 'dndsR' "$target"; then
    echo "$target exists and is not a dndsR hook. Merge it by hand; nothing was changed." >&2
    exit 1
  fi
done
cp tools/dev/githooks/* tools/dev/check-public-safe.sh "$hooks/"

if [ -n "$audit_target" ]; then
  if [ -e .audit ] && [ ! -L .audit ]; then
    mkdir -p "$audit_target"
    cp -pRu .audit/. "$audit_target/"
    mv .audit ".audit.moved.$(date +%s)"
    echo "Copied the existing audit record to $audit_target; old copy renamed .audit.moved.*"
  fi
  mkdir -p "$audit_target"
  ln -sfn "$audit_target" .audit
fi
mkdir -p .audit/transcripts
chmod 700 "$(readlink -f .audit)"
echo "hooks installed in $hooks (rerun this script after they change)"
echo ".audit -> $(readlink -f .audit)"
