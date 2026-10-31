# Shared helpers for dndsR git hooks. Sourced, not executed.
ZERO=0000000000000000000000000000000000000000
PROTECTED_BRANCH=refs/heads/main
DEV_BRANCH=refs/heads/dev

block() {
  printf 'dndsR policy: BLOCKED: %s\n' "$*" >&2
  printf 'See docs/dev/RUNBOOK.md. Changes to main or to tags need explicit authorization from the maintainer.\n' >&2
  exit 1
}

# Append "<UTC time> <event> <sha> <branch> <detail>" to the untracked audit log.
log_git_event() {
  top=$(git rev-parse --show-toplevel 2>/dev/null) || return 0
  [ -n "$top" ] || return 0
  mkdir -p "$top/.audit" 2>/dev/null || return 0
  branch=$(git symbolic-ref --short -q HEAD || echo detached)
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$2" "$branch" \
    "$(printf '%s' "$3" | tr '\t\n' '  ')" >>"$top/.audit/git-events.tsv" 2>/dev/null || true
}
