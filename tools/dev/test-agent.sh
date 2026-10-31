#!/usr/bin/env bash
# End-to-end agent test: a new agent gets a clean clone at HEAD and a plain request, with
# no project hints. Its results are then checked by script, not from its own claims.
#   DNDSR_SIF=/path/image.sif tools/dev/test-agent.sh <claude|codex> [model] [output dir]
# Agent credentials stay in their usual place; only HOME is replaced, so nothing is
# installed into the real home directory. The agent receives a local image to avoid a
# 1.3 GB download. This environment detail is the one hint it is given.
set -u
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
agent=${1:?usage: test-agent.sh claude|codex [model] [output dir]}
model=${2:-}
out=${3:-$(mktemp -d)}
: "${DNDSR_SIF:?Set DNDSR_SIF to the runtime image}"
real_home=$HOME
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/home" "$out"
git clone -q "$root" "$work/dndsR"
cd "$work/dndsR" || exit 1

prompt="I just cloned this repository. Set up dndsR, run its documented example, confirm the results are right, and explain to me what the results mean and what I should be careful about. A copy of the container image is at $DNDSR_SIF, so you do not need to download it."

export HOME=$work/home
case $agent in
  claude)
    CLAUDE_CONFIG_DIR=$real_home/.claude timeout 1500 claude -p ${model:+--model "$model"} \
      --allowedTools "Bash Read Write Edit Skill" --output-format json "$prompt" >"$out/claude.json" 2>"$out/stderr.txt"
    rc=$?
    jq -r '.result // empty' "$out/claude.json" >"$out/answer.md"
    ;;
  codex)
    CODEX_HOME=$real_home/.codex timeout 1500 codex exec --skip-git-repo-check -s danger-full-access "$prompt" >"$out/answer.md" 2>"$out/stderr.txt"
    rc=$?
    ;;
  *) echo "unknown agent: $agent" >&2; exit 2 ;;
esac
echo "agent exit code: $rc"

# Independent checks of what is actually on disk.
pass=0 fail=0
ok() { if eval "$2"; then pass=$((pass + 1)); echo "ok    $1"; else fail=$((fail + 1)); echo "FAIL  $1"; fi; }
find_result() { find "$work" -name example_AvB_dnds.tsv -print -quit; }
result=$(find_result)
ok "dndsr was installed" '[ -x "$HOME/.local/bin/dndsr" ]'
ok "the example produced a dN/dS table" '[ -n "$result" ]'
if [ -n "$result" ]; then
  example=$(dirname "$(dirname "$result")")
  (cd "$example" && "$HOME/.local/bin/dndsr" verify_example) >"$out/verify.txt" 2>&1
  ok "the results pass an independent run of verify_example" 'grep -q "All 14 example checks passed" "$out/verify.txt"'
fi
ok "the answer warns that a ratio above 1 is not proof" 'grep -Eiq "(not|never|no|isn.t|without).{0,60}(proof|prove|evidence|proves)|candidate" "$out/answer.md"'
ok "the answer mentions dS caveats" 'grep -Eiq "\bdS\b" "$out/answer.md"'
ok "the answer says the example is synthetic" 'grep -Eiq "synthetic|designed|simulated" "$out/answer.md"'
ok "the answer says the example does not validate real data" 'grep -Eiq "(your own|real|other).{0,40}(data|genome)" "$out/answer.md"'
echo "passed $pass, failed $fail; answer in $out/answer.md"
[ "$fail" -eq 0 ]
