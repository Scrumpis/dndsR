#!/usr/bin/env bash
# PreToolUse guard for Bash. Blocks commands that bypass or break repository policy.
# Pattern matching is best effort; the git hooks in tools/dev/githooks are the second layer.
set -u
command -v jq >/dev/null || { echo "git-guard: jq is required; refusing to run Bash unguarded" >&2; exit 2; }
cmd=$(jq -r '.tool_input.command // empty')
[ -n "$cmd" ] || exit 0

# Collapse global git options so "git -c x=y -C dir push" reads as "git push".
norm=$(printf '%s\n' "$cmd" | sed -E ':a; s/\bgit[[:space:]]+(-c|-C)[[:space:]]+[^[:space:]]+/git/; ta')
deny() { printf 'dndsR policy: %s. (see docs/dev/RUNBOOK.md).\n' "$1" >&2; exit 2; }
has() { printf '%s\n' "$norm" | grep -Eq -- "$1"; }

# Bypassing safeguards
has '\-\-no-verify' && deny 'bypassing git hooks is not allowed'
has '\bgit[[:space:]]+commit([[:space:]]+-[a-zA-Z]+)*[[:space:]]+-[a-zA-Z]*n' && deny 'bypassing git hooks is not allowed'
printf '%s\n' "$cmd" | grep -Eqi 'core\.hookspath' && deny 'changing core.hooksPath is not allowed'
has '\.git/hooks' && deny 'editing .git/hooks is not allowed'

# Pushes belong to the maintainer: commit, then report the exact push command.
has '\bgit[[:space:]]+push\b' && deny 'agents do not push; commit to dev and give the maintainer the push command'

# main and tag reconciliation
has '\bgit[[:space:]]+(merge|rebase|cherry-pick)\b[^;&|]*\b(origin/)?main\b' && deny 'reconciling with main is not allowed'
has '\bgit[[:space:]]+(pull)\b[^;&|]*\bmain\b' && deny 'pulling main into a working branch is not allowed'
has '\bgit[[:space:]]+branch\b[^;&|]*(-[dDfmM]|--delete|--force|--move)[^;&|]*\bmain\b' && deny 'altering local main is not allowed'
has '\bgit[[:space:]]+(checkout|switch)[[:space:]]+(-B|-C)[[:space:]]+main\b' && deny 'altering local main is not allowed'
has '\bgit[[:space:]]+update-ref\b[^;&|]*refs/(heads/main|tags/)' && deny 'altering main or tags is not allowed'
has '\bgit[[:space:]]+tag\b[^;&|]*(-[a-zA-Z]*[fd]\b|--force|--delete)' && deny 'moving or deleting tags is not allowed'
has '\bgit[[:space:]]+filter-(branch|repo)\b' && deny 'rewriting history is not allowed'

# Publishing and repository settings
has '\bgh[[:space:]]+(release|repo[[:space:]]+(edit|delete)|pr[[:space:]]+merge|ruleset)\b' && deny 'publishing or changing GitHub settings is not allowed'
has '\bgh[[:space:]]+api\b.*(-X|--method)[[:space:]]*(POST|PUT|PATCH|DELETE)' && deny 'writing through the GitHub API is not allowed'
has '\b(docker|podman)[[:space:]]+push\b|\b(singularity|apptainer)[[:space:]]+push\b' && deny 'publishing container images is not allowed'

# Protecting the audit record
has '\brm\b[^;&|]*(\.audit|\.git)([/[:space:]]|$)' && deny 'deleting .audit or .git is not allowed'
exit 0
