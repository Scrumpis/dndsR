#!/usr/bin/env bash
# Human-readable view of .audit/events.jsonl.
#   audit-report.sh              one row per session
#   audit-report.sh <session-id> commands, edits, tests, commits, and pushes for a session
set -eu
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
log=$root/.audit/events.jsonl
[ -s "$log" ] || { echo "No audit events at $log"; exit 1; }
command -v jq >/dev/null || { echo "jq is required" >&2; exit 1; }

if [ $# -eq 0 ]; then
  printf 'SESSION\tFIRST\tLAST\tBRANCH\tHEAD\tTOOLS\tFAILED\tTESTS\tCOMMITS\tPUSHES\tENDED\n'
  jq -rs 'group_by(.session_id)[] | sort_by(.ts) as $e | [
    $e[0].session_id, $e[0].ts, $e[-1].ts, $e[-1].branch, ($e[-1].head // "")[:8],
    ([$e[] | select(.event == "PostToolUse" or .event == "PostToolUseFailure")] | length),
    ([$e[] | select(.event == "PostToolUseFailure")] | length),
    ([$e[] | select(.category == "test")] | length),
    ([$e[] | select(.category == "git_commit" and .event == "PostToolUse")] | length),
    ([$e[] | select(.category == "git_push" and .event == "PostToolUse")] | length),
    ([$e[] | select(.event == "SessionEnd") | .reason] | last // "-")
  ] | @tsv' "$log"
  exit 0
fi

jq -r --arg s "$1" 'select(.session_id == $s) |
  [.ts, .event, (.agent_type // "main"), (.head // "")[:8],
   (if .tool then (.tool + " ") else "" end) +
   ((.command // .file // .source // .reason // .agent_id // "") | gsub("\n"; " ") | .[:160])] | @tsv' "$log"
