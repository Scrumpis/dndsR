# Git policy

- **Before you edit anything, sync.** Fetch `dev` from GitHub with the remote that works on this machine: check `git remote -v`, and if `git fetch origin` fails (for example an SSH key GitHub rejects), fetch over HTTPS, `git fetch https://github.com/Scrumpis/dndsR.git '+refs/heads/*:refs/remotes/origin/*'`. Then bring your local `dev` up to date:
  - Nothing local and unpushed: fast-forward (`git merge --ff-only origin/dev`).
  - Unpushed commits that nothing cites: rebase them onto `origin/dev`, or merge; either is fine.
  - Unpushed commits whose identity matters, because a frozen study, a provenance file, a results root, or a document cites the commit: **merge, do not rebase**. A rebase gives them new hashes and breaks those references.
  - Then run `tools/dev/setup.sh` (safe to repeat), so the installed git hooks match the synced sources. Hook copies do not update themselves, and `test-safeguards.sh` fails on a stale one. Install only hooks that are on `origin/dev`: if `git diff --quiet origin/dev -- tools/dev/githooks tools/dev/check-public-safe.sh` fails, stop and report instead of running it.
  - Then read `docs/dev/channels/` and `docs/dev/handoffs/` before the first edit, so messages from another side are seen before code is written.
- `main` is protected. `dev` holds all development. Tags and releases are immutable.
- Agents may commit to `dev`, but the maintainer does all pushing. Report the exact push command at the end of the work.
- Unless the maintainer authorizes it for this specific action, do not merge, rebase, or cherry-pick between `main` and `dev`; do not create, move, or delete release tags; do not publish releases or images; do not change GitHub settings.
- Never force-push, rewrite pushed history, or bypass hooks (`--no-verify`, `core.hooksPath`). If a safeguard blocks you, stop and report it rather than working around it.
- Keep commits coherent: software, tests, skill, and docs for one change go together.
- If a commit changes software without changing `skills/`, add the trailer `Skill-Impact: none - <reason>` after actually reviewing the skill.
- End commit messages with the attribution the session specifies.
