# Git policy

- **Before you edit anything, sync.** `git fetch origin`, then `git pull --rebase` on `dev` (or `git rebase origin/dev` when the pull refuses a divergent branch). A stale tree is how two agents overwrite each other: work lands on `dev` between sessions, and a rebase of your own unpushed commits is expected, not an error. Resolve conflicts; do not substitute a merge commit. Then read `docs/dev/channels/` and `docs/dev/handoffs/` before the first edit, so messages from another side are seen before code is written.
- `main` is protected. `dev` holds all development. Tags and releases are immutable.
- Agents may commit to `dev`, but the maintainer does all pushing. Report the exact push command at the end of the work.
- Unless the maintainer authorizes it for this specific action, do not merge, rebase, or cherry-pick between `main` and `dev`; do not create, move, or delete release tags; do not publish releases or images; do not change GitHub settings.
- Never force-push, rewrite pushed history, or bypass hooks (`--no-verify`, `core.hooksPath`). If a safeguard blocks you, stop and report it rather than working around it.
- Keep commits coherent: software, tests, skill, and docs for one change go together.
- If a commit changes software without changing `skills/`, add the trailer `Skill-Impact: none - <reason>` after actually reviewing the skill.
- End commit messages with the attribution the session specifies.
