#!/usr/bin/env bash
# Exercises the git hooks and the Claude Code guard in a disposable clone with a local
# bare remote. Never touches this repository's branches, tags, or real remote.
set -u
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
pass=0 failed=0
ok() { pass=$((pass + 1)); printf 'ok    %s\n' "$1"; }
no() { failed=$((failed + 1)); printf 'FAIL  %s\n' "$1"; }

git clone -q --bare "$root" "$tmp/remote.git"
git clone -q --branch dev "$tmp/remote.git" "$tmp/work"
cd "$tmp/work" || exit 1
git config user.name test; git config user.email test@example.invalid
git branch -q main origin/main
# dev may already contain main, which makes a merge of main a no-op that proves nothing. Keep the
# fixture's main ahead of dev, in this clone and its remote only.
if git merge-base --is-ancestor main dev; then
  git branch -qf main "$(git commit-tree -p main -m "fixture: main ahead of dev" "main^{tree}")"
  git push -q origin main
fi

# Copy the current (possibly uncommitted) safeguards into the clone.
mkdir -p tools .claude
cp -a "$root/tools/dev" tools/
cp -a "$root/skills" "$root/.agents" .; cp -a "$root/.claude/skills" .claude/; cp "$root/.gitignore" .
mkdir -p docs/tutorials; cp -a "$root/docs/tutorials/." docs/tutorials/
git add -A; git commit -qm "test fixture"
git push -q origin dev
cp tools/dev/githooks/* tools/dev/check-public-safe.sh .git/hooks/
git remote set-url origin "$tmp/remote.git"
main0=$(git rev-parse main); tag0=$(git rev-parse v0.2.0); dev0=$(git rev-parse dev)

blocked() { # description, command...
  local desc=$1; shift
  if out=$("$@" 2>&1); then no "$desc (was allowed)"
  elif grep -Eq 'dndsR policy|public-safety' <<<"$out"; then ok "blocked: $desc"
  else no "$desc (failed for another reason: $(tail -1 <<<"$out"))"; fi
  git merge --abort 2>/dev/null; git rebase --abort 2>/dev/null
  git reset -q --hard 2>/dev/null; git clean -qfd 2>/dev/null
  git checkout -q dev 2>/dev/null
}
allowed() {
  local desc=$1; shift
  if out=$("$@" 2>&1); then ok "allowed: $desc"; else no "$desc (was blocked: $(tail -1 <<<"$out"))"; fi
}
sameref() { [ "$(git rev-parse "$1")" = "$2" ] && ok "unchanged: $1" || no "$1 moved"; }

on_main() { git checkout -q main && "$@"; }
commit_file() { printf '%s\n' "$2" >"$1"; git add -f "$1"; git commit -qm "test $1" "${@:3}"; }
secret() { commit_file a.txt "$(printf 'token %s%s' 'ghp_' 'abcdefghijklmnopqrstuvwxyz')"; }
private_path() { commit_file b.txt "$(printf '/mnt/%s%s/user/x' 'gs' '21')"; }
software() { commit_file R/zz_test.R "# test" "$@"; }

blocked "commit on main"              on_main git commit -q --allow-empty -m x
blocked "merge dev into main"         on_main git merge -q --no-edit dev
git branch -q ff-test main
git update-ref refs/heads/ff-test "$(git commit-tree -p main -m ff 'main^{tree}')"
blocked "fast-forward merge into main" on_main git merge -q --ff-only ff-test
blocked "reset main"                  on_main git reset -q --hard dev
blocked "branch -f main"              git branch -f main dev
blocked "branch -D main"              git branch -D main
blocked "update-ref main"             git update-ref refs/heads/main dev
blocked "merge main into dev"         git merge -q --no-edit main
blocked "rebase dev onto main"        git rebase -q main
blocked "rebase rewriting pushed dev" git rebase -q --onto HEAD~3 HEAD~1
blocked "move tag"                    git tag -f v0.2.0 dev
blocked "delete tag"                  git tag -d v0.2.0
blocked "update-ref tag"              git update-ref refs/tags/v0.2.0 dev
blocked "fast-forward push to main"   git push -q origin ff-test:main
blocked "force-push older dev"        git push -q -f origin HEAD~2:dev
git tag v9.9.9-test
blocked "push a tag"                  git push -q origin v9.9.9-test
blocked "delete remote dev"           git push -q origin :dev
blocked "commit a credential"         secret
blocked "commit a private path"       private_path
blocked "commit .sif file"            commit_file big.sif x
blocked "commit software, no skill"   software
blocked "commit a transcript"         commit_file session.jsonl '{}'

allowed "commit on dev"               git commit -q --allow-empty -m "allowed"
allowed "software with Skill-Impact"  software -m "Skill-Impact: none - test"
allowed "fast-forward push to dev"    git push -q origin dev
allowed "fetch"                       git fetch -q origin
allowed "pack-refs"                   git pack-refs --all
allowed "gc"                          git gc -q
blocked "delete tag after packing"    git tag -d v0.2.0
blocked "delete main after packing"   git branch -D main
allowed "checkout main (read only)"   git checkout -q main
blocked "commit on main, second time" git commit -q --allow-empty -m x
git checkout -q dev
git commit -q --allow-empty -m "local"
allowed "rebase unpushed commits"     git rebase -q origin/dev

sameref main "$main0"; sameref v0.2.0 "$tag0"
[ "$(git -C "$tmp/remote.git" rev-parse main)" = "$main0" ] && ok "remote main unchanged" || no "remote main moved"
[ "$(git -C "$tmp/remote.git" rev-parse v0.2.0)" = "$tag0" ] && ok "remote tag unchanged" || no "remote tag moved"

# The maintainer's "Update branch" on GitHub: dev may follow it, but not make its own.
git reset -q --hard origin/dev
git clone -q --branch dev "$tmp/remote.git" "$tmp/github"
git -C "$tmp/github" -c user.name=m -c user.email=m@example.invalid merge -q --no-edit origin/main
git -C "$tmp/github" push -q origin dev
allowed "fetch the maintainer's merge" git fetch -q origin
blocked "merge main locally, once fetched" git merge -q --no-ff --no-edit main
allowed "follow the maintainer's merge" git pull -q --ff-only origin dev
[ "$(git rev-parse dev)" = "$(git rev-parse origin/dev)" ] && ok "dev equals origin/dev" || no "dev did not follow origin/dev"
[ -s .audit/git-events.tsv ] && ok "commit and push events logged" || no "git events not logged"

# The hooks that are really installed in this repository must match their sources.
real=$(git -C "$root" rev-parse --absolute-git-dir)/hooks
for f in "$root"/tools/dev/githooks/* "$root"/tools/dev/check-public-safe.sh; do
  cmp -s "$f" "$real/$(basename "$f")" && ok "installed: $(basename "$f")" || no "stale or missing hook, run tools/dev/setup.sh: $(basename "$f")"
done
[ -z "$(git -C "$root" config --local --get core.hooksPath)" ] && ok "default hooks path in use" || no "core.hooksPath is overridden"

# Claude Code guard
guard=$root/.claude/hooks/git-guard.sh
g() { jq -n --arg c "$1" '{tool_input:{command:$c}}' | "$guard" >/dev/null 2>&1; }
for c in 'git push origin dev' 'git -c a=b push origin dev' 'git commit --no-verify -m x' 'git commit -n -m x' \
  'git -c core.hooksPath=/dev/null commit -m x' 'git config --unset core.hooksPath' 'git merge origin/main' \
  'git rebase main' 'git pull origin main' 'git branch -f main dev' 'git checkout -B main' \
  'git update-ref refs/tags/v0.2.0 HEAD' 'git tag -f v0.2.0' 'git tag -d v0.2.0' 'gh release create v1' \
  'gh pr merge 1' 'docker push scrumpis/dndsr:x' 'cd x && rm -rf .audit' 'rm -rf .git'; do
  g "$c"; [ $? -eq 2 ] && ok "guard denies: $c" || no "guard allowed: $c"
done
for c in 'git status' 'git log main..dev' 'git diff origin/main origin/dev' 'git commit -m "x"' \
  'git commit --no-edit' 'git fetch origin' 'git tag -l' 'Rscript -e "1"' 'rm -f .gitignore.bak'; do
  g "$c"; [ $? -eq 0 ] && ok "guard allows: $c" || no "guard denied: $c"
done

# Read-only hook for audit agents
readonly_hook=$root/.claude/hooks/read-only-bash.sh
r() { jq -n --arg c "$1" --arg a "$2" '{agent_type:$a,tool_input:{command:$c}}' | "$readonly_hook" >/dev/null 2>&1; }
for c in 'git log -1 --oneline' 'grep -n x R/a.R | head -5' 'git diff main dev 2>/dev/null' 'ls -la && wc -l R/*.R'; do
  r "$c" software-auditor; [ $? -eq 0 ] && ok "auditor may run: $c" || no "auditor denied: $c"
done
for c in 'touch x' 'echo a > x' 'cat a >> b' 'Rscript -e 1' 'git commit -m x' 'find . -delete' 'ls; rm x' 'ls $(rm x)'; do
  r "$c" scientific-auditor; [ $? -eq 2 ] && ok "auditor may not run: $c" || no "auditor allowed: $c"
done
r 'touch x' implementation-agent; [ $? -eq 0 ] && ok "other agents unaffected" || no "read-only hook hit another agent"

printf '\n%d passed, %d failed\n' "$pass" "$failed"
[ "$failed" -eq 0 ]
