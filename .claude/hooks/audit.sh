#!/usr/bin/env bash
# Records one JSON line per Claude Code hook event in .audit/events.jsonl and keeps
# .audit/transcripts/ in step with the native transcripts. Logging must never stop
# work, so every path exits 0. Full tool output stays in the native transcript.
set -u
umask 077
input=$(cat)
command -v jq >/dev/null || { echo "audit: jq not found; event not logged" >&2; exit 0; }
root=${CLAUDE_PROJECT_DIR:-$(pwd)}
audit=$root/.audit
mkdir -p "$audit/transcripts" 2>/dev/null || exit 0

ref() { git -C "$root" rev-parse -q --verify "$1" 2>/dev/null; }
branch=$(git -C "$root" symbolic-ref --short -q HEAD 2>/dev/null || echo detached)

jq -c --arg ts "$(date -u +%Y-%m-%dT%H:%M:%SZ)" --arg branch "$branch" \
  --arg head "$(ref HEAD)" --arg origin_dev "$(ref refs/remotes/origin/dev)" '
  (.tool_input.command // "") as $c
  | {ts: $ts, event: .hook_event_name, session_id, agent_id, agent_type,
     branch: $branch, head: $head, origin_dev: $origin_dev,
     cwd, transcript_path, agent_transcript_path, permission_mode, model,
     source, reason, trigger,
     tool: .tool_name, tool_use_id, duration_ms,
     command: .tool_input.command, description: .tool_input.description,
     file: (.tool_input.file_path // .tool_input.notebook_path),
     error: (.error | if type == "string" then .[:2000] else . end),
     summary: (.last_assistant_message | if type == "string" then .[:2000] else . end),
     category: (if $c == "" then null
       elif ($c | test("\\bgit\\b[^;&|]*\\bpush\\b")) then "git_push"
       elif ($c | test("\\bgit\\b[^;&|]*\\bcommit\\b")) then "git_commit"
       elif ($c | test("testthat|devtools::(test|check)|rcmdcheck|R CMD check|test-safeguards|validate-skill")) then "test"
       else "command" end)}
  | with_entries(select(.value != null and .value != ""))' \
  <<<"$input" >>"$audit/events.jsonl" 2>/dev/null

case $(jq -r '.hook_event_name // empty' <<<"$input") in
  SessionStart | SubagentStop | PreCompact | Stop | SessionEnd)
    transcript=$(jq -r '.transcript_path // empty' <<<"$input")
    [ -n "$transcript" ] || exit 0
    dir=$(dirname "$transcript")
    # A subagent event may report .../<session>/subagents/agent-*.jsonl.
    [ "$(basename "$dir")" != subagents ] || dir=$(dirname "$(dirname "$dir")")
    [ -d "$dir" ] || exit 0
    # Compare content, never timestamps: the archive may be recursively touched to
    # survive the scratch purge, and that must not stop updated files being copied.
    (cd "$dir" && find . -type f -print0) | while IFS= read -r -d '' f; do
      target=$audit/transcripts/$f
      cmp -s "$dir/$f" "$target" && continue
      mkdir -p "$(dirname "$target")" && cp -p "$dir/$f" "$target"
    done
    ;;
esac
exit 0
