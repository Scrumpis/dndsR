# Development runbook

## Resume work

```
cd <clone>
tools/dev/setup.sh                 # safe to repeat
git status && git log --oneline -5
tools/dev/audit-report.sh          # sessions so far
claude --continue                  # or: claude --resume <session-id>
```

Then read `docs/dev/STATUS.md`, `docs/dev/DECISIONS.md`, and the latest entry in `docs/dev/SESSION_LOG.md`. A fresh clone has no `.audit` and needs `tools/dev/setup.sh` before committing.

## Commit and push

Agents commit to `dev` and report the push command. The maintainer pushes:

```
git push origin dev
```

Pushing needs GitHub credentials on the machine that pushes.

## Run code in the container

The HPC default R lacks many dependencies, so use the image. Keep images and caches outside the clone. This example installs the current tree into a temporary library:

```
SIF=<path to dndsr image>
LIB=$(mktemp -d)
singularity exec --cleanenv --bind "$PWD,$LIB" "$SIF" R CMD INSTALL --library="$LIB" "$PWD"
```

## Inspect the audit trail

Everything below is in `.audit/`, which is excluded from Git, Docker contexts, and R builds.

| What | Where |
|---|---|
| Hook events: sessions, branch and HEAD, commands, edits, failures, and test, commit, and push categories | `events.jsonl` |
| Commits and push attempts from Git, including those made outside Claude | `git-events.tsv` |
| Copies of native transcripts | `transcripts/<session-id>.jsonl` |
| Subagent transcripts | `transcripts/<session-id>/subagents/` |

```
tools/dev/audit-report.sh                 # one row per session
tools/dev/audit-report.sh <session-id>    # timeline for a session
jq -c 'select(.category == "test")' .audit/events.jsonl
```

Native originals are in `~/.claude/projects/<encoded-project-path>/`, and Claude Code deletes them after `cleanupPeriodDays` (30 days by default). The archive keeps copies: every `SessionStart`, `Stop`, `SubagentStop`, `PreCompact`, and `SessionEnd` copies new or changed files from the project's transcript directory.

Tracked, public-safe records: `STATUS.md`, `DECISIONS.md`, and `SESSION_LOG.md` in `docs/dev/`. They hold summaries, not transcripts.

### Audit limitations

- The log records what the hooks saw, and tool output is only in the transcript. The hooks do not capture work done by a person outside Claude, apart from `git-events.tsv`.
- A command's internal actions are invisible unless it reports them: `Rscript` can change files without an `Edit` event.
- `category` is a text match on the command, so it can mislabel. A Bash failure appears as `PostToolUseFailure`.
- A push made by the maintainer outside Claude appears only as a `push-attempt` in `git-events.tsv`, with no outcome. The next `git fetch` shows the outcome through `origin_dev` in later events.
- Commands denied by a guard or permission rule never run, so they appear only in the transcript, not in `events.jsonl`.
- `SessionEnd` does not run if the process is killed. Earlier `Stop` copies and the next session's sweep cover most of that gap; the event log itself then shows no end.
- Events before the hooks were installed exist only in the transcript. The bootstrap session's early work is one such case.
- `.audit` has no tamper protection. It is a local record, kept private with mode 700. Transcripts may contain anything shown in a session, so treat them as sensitive and do not share them unreviewed.
- If the clone is on a scratch filesystem with a purge policy, the archive can be purged with it. Move it with `tools/dev/setup.sh --audit-dir <durable path>`, which leaves `.audit` as a symlink.

## Safeguards

| Layer | Enforces | Bypass and limits |
|---|---|---|
| GitHub ruleset `protect_main` (active) | On `main`: no deletion or force-push, pull request with 1 approval, required `R-CMD-check` | Repository admins may be able to bypass it; the public API does not show bypass actors. Not applied to `dev` or tags. |
| Git hooks, sourced from `tools/dev/githooks` and copied to `.git/hooks` | Commits and ref updates on `main`; bringing `main` into `dev` by merge or rebase; rewriting pushed `dev`; moving or deleting tags; pushing anything but a fast-forward of `dev`; unsafe staged files; the skill gate | Local to a clone. Rerun `tools/dev/setup.sh` after hooks change, because copies do not update themselves (`test-safeguards.sh` detects stale copies). `--no-verify` skips the commit and push hooks. `-c core.hooksPath=...` skips all of them. Cherry-picks from `main` are not detected. |
| Claude Code guard and `permissions` in `.claude/` | Denies agent pushes and bypass flags, publishing, and reconciliation with `main` | Pattern matching only; commands hidden in scripts or `Rscript -e` can evade it. It reads the whole command text, so a heredoc that merely mentions a blocked command is also denied: write such text with the file tools. Only applies inside Claude Code. |

Nothing protects `dev` or tags on GitHub. Fast-forward merges to local `main` are blocked by the hooks, so read `origin/main` after `git fetch` instead. When the maintainer deliberately updates local `main`, use `git -c core.hooksPath=/dev/null ...`; the guard stops agents from doing so.

Validate with `tools/dev/test-safeguards.sh`, which uses a disposable clone and local remote.

## Permissions

`.claude/settings.json` uses `acceptEdits` with a scoped allowlist, without bypass mode. `Rscript` and container commands are allowed so tests can run, which means they can do more than their names suggest; review is by audit trail. Project `allow` rules are ignored until someone starts `claude` interactively in the clone and accepts the workspace trust prompt; `deny` rules and hooks were confirmed active before that. Personal additions go in the untracked `.claude/settings.local.json`.

## Agents

| Agent | Mode | Writes |
|---|---|---|
| `scientific-auditor` | No edit tools; Bash limited to a read-only list | None; returns a report |
| `software-auditor` | No edit tools; Bash limited to a read-only list | None; returns a report |
| `validation-engineer` | Normal | `tests/`, `docs/dev/`, `tools/dev/` (by instruction only) |
| `implementation-agent` | Normal | Approved changes on `dev` |

Audit reports should be saved by the calling session under `docs/dev/audits/` if they are public-safe.

`.claude/hooks/read-only-bash.sh` enforces the auditors' limit from project settings. In testing, plan mode and a hook in agent frontmatter did not stop a `touch`, so do not rely on them. The list is not a sandbox: an allowed tool could still have an unlisted way to write. The auditors also cannot run R, so numerical checks go through the calling session or `validation-engineer`. The `validation-engineer` path limit is by instruction only.
