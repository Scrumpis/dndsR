# Session log

Summaries only. Full transcripts are local: see `RUNBOOK.md`. Newest first.

## 2026-10-31: Phase 0 bootstrap

- **Session:** `61d9f4da-7419-4c58-bb02-9ebdad0c13e6`
- **Start:** fresh clone, clean `dev` at `c8b534f`
- **Done:** git safeguards, Claude Code hooks and permissions, audit trail, four development agents, end-user skill and tutorial, validation tooling, and project docs. No scientific code changed.
- **Verified:** see `skills/dndsr/references/verification.md` and `tools/dev/test-safeguards.sh`.
- **Failed or not done:** this machine has no GitHub credentials, so a push dry-run failed; the maintainer now handles pushes. No analysis was run.
- **Found by testing:** the first hook design left `main` unprotected once it was checked out, missed merges from `main`, and blocked `git pack-refs`. All were redesigned and retested. A real subagent run also showed that plan mode and a frontmatter hook left the auditors able to write through Bash, so a project-level hook now enforces read-only Bash; the probe left one empty file, which was removed. The guard also falsely denied a command that wrote documentation mentioning a blocked flag; those files were written with the file tools and the behavior is documented.
