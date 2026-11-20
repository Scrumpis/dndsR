# Development runbook

## Resume work

```
cd <clone>
git fetch https://github.com/Scrumpis/dndsR.git '+refs/heads/*:refs/remotes/origin/*'   # or: git fetch origin
git merge --ff-only origin/dev     # see git-policy for unpushed commits
tools/dev/setup.sh                 # after every sync: reinstalls the hooks; safe to repeat
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
export DNDSR_SIF=<path to dndsr image> DNDSR_DEV_LIB=<writable folder outside the clone>
tools/dev/in-container.sh Rscript -e 'packageVersion("dndsR")'   # installs this tree first
tools/dev/run-tests.sh
```

All the checks, side by side, on the CPUs of a Slurm allocation:

```
srun -c 8 --mem=48G --time=10 tools/dev/check.sh          # about 3 min: run before committing
srun -c 8 --mem=24G --time=10 tools/dev/check.sh --fast   # leaves out the tutorial and slowest tests, while developing
srun -c 8 --mem=48G --time=20 tools/dev/check.sh --rcmd-check   # adds R CMD check, as CI runs it: before a pull request
srun -c 8 --mem=48G --time=20 tools/dev/check.sh --regression   # adds the real-genome fixture regression: changes to R/
srun -c 8 --mem=32G --time=15 tools/dev/regression.sh --threads 8   # the fixture regression alone, about 2 min
tools/dev/rcmd-check.sh --as-cran                         # R CMD check alone; extra arguments go to it
tools/dev/run-tests.sh --fast enrichment                  # only the files that match
```

Which checks a change needs, and when whole genomes are rerun: [TESTING.md](TESTING.md). With 8 workers, 24 GB is no longer enough: on 2026-11-20 an orthologr worker was killed for memory (`test-threads.R`, an `unserialize` error), so ask for 48 GB. A compute node's `/tmp` is its own, so give `regression.sh --keep` a folder on shared storage to inspect a run from the login node.

Test files run in parallel, using `SLURM_CPUS_PER_TASK`, or `DNDSR_TEST_CPUS`. A file runs on one worker, so the longest file sets the wall time: keep slow tests in files of their own. `R CMD check` still runs them one at a time.

On a machine with limited memory, or WSL with its default limits, cap the workers and run the checks one at a time. Each R worker holds several hundred MB once Bioconductor and orthologr are loaded, and the default is one worker per CPU:

```
DNDSR_TEST_CPUS=4 tools/dev/check.sh --fast --serial
DNDSR_TEST_CPUS=4 tools/dev/run-tests.sh --fast enrichment   # one filter at a time
```

`--serial` runs the tests, skill, and safeguards (and the tutorial, without `--fast`) one after another instead of side by side.

`--rcmd-check` builds the tarball and runs `R CMD check --no-manual` on it in the image, failing on an ERROR or WARNING as CI does. It catches what only check sees, such as a namespace note or a documentation mismatch. About 95 s alone, longer beside the tests.

On some compute nodes `git` itself fails with "Illegal instruction", which stops the safeguards and the tutorial (both call `git`) before any test runs. Run those two on a login node: `tools/dev/test-safeguards.sh`, `tools/dev/test-tutorial.sh`.

Runtime files can be recreated if purged: pull the image by the digest in `skills/dndsr/SKILL.md`. For the estimator comparison, create a micromamba environment with `paml` and `kakscalculator2` from bioconda, and for regenerating docs, install roxygen2 into `DNDSR_DEV_LIB`.

Individual checks: `tools/dev/run-tests.sh`, `tools/dev/test-tutorial.sh`, `tools/dev/test-safeguards.sh`, `tools/dev/test-agent.sh claude`, and `tools/dev/validation/compare-estimators.R`.

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

Timestamps in the archive are not meaningful, because files may be touched to survive the scratch purge. The archive sync compares content, and `tools/dev/audit-manifest.sh` records SHA-256 hashes; run it with `--check` to look for unrecorded content.

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
- This clone's archive is in the home folder, mode 700, outside the scratch purge: `.audit` is a symlink to `~/.local/share/dndsr-audit`. The earlier copy stays beside the clone as `.audit.moved.*`. A new clone on scratch should do the same: `tools/dev/setup.sh --audit-dir <durable path>`.

## Safeguards

| Layer | Enforces | Bypass and limits |
|---|---|---|
| GitHub ruleset `protect_main` (active) | On `main`: no deletion or force-push, pull request with 0 approvals, required `R-CMD-check`, `validate`, `conda`, and `build-and-test` | Repository admins may be able to bypass it; the public API does not show bypass actors. |
| GitHub ruleset `protect_dev` (active) | On `dev`: no deletion or force-push | The same bypass limit. |
| Git hooks, sourced from `tools/dev/githooks` and copied to `.git/hooks` | Commits and ref updates on `main`; bringing `main` into `dev` by merge or rebase; rewriting pushed `dev`; moving or deleting tags; pushing anything but a fast-forward of `dev`; unsafe staged files; the skill gate | Local to a clone. Agents rerun `tools/dev/setup.sh` after every sync (git-policy), because copies do not update themselves (`test-safeguards.sh` detects stale copies). `--no-verify` skips the commit and push hooks. `-c core.hooksPath=...` skips all of them. Cherry-picks from `main` are not detected. |
| Claude Code guard and `permissions` in `.claude/` | Denies agent pushes and bypass flags, publishing, and reconciliation with `main` | Pattern matching only; commands hidden in scripts or `Rscript -e` can evade it. It reads the whole command text, so a heredoc that merely mentions a blocked command is also denied: write such text with the file tools. Only applies inside Claude Code. |

`protect_tags` is also active; a tag ruleset's rules are not readable without a login. Fast-forward merges to local `main` are blocked by the hooks, so read `origin/main` after `git fetch` instead. When the maintainer deliberately updates local `main`, use `git -c core.hooksPath=/dev/null ...`; the guard stops agents from doing so.

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
