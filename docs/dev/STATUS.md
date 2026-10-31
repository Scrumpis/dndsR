# Project status

Updated 2026-10-31. Update this at the end of each substantive session.

**Phase:** Phase 0 (infrastructure) complete. Phase 1 (audits) has not started and needs explicit instruction.

## Baseline

| Ref | Commit | Notes |
|---|---|---|
| `v0.2.0` | `b8ea5ff` | Prerelease on GitHub. `DESCRIPTION` version 0.2.0.000 |
| `main` | `bbd8384` | 8 commits past `v0.2.0` |
| `dev` at start | `c8b534f` | 134 commits past `v0.2.0`, version 0.3.0.000 |

`dev` is 128 commits ahead of `main` and 2 behind. The 2 are `d9a2279` and `bbd8384`, which add and then delete `R/dnds_summary.R` on `main`, so `main`'s content equals the merge base (`9094341`). They have no net effect. A read-only `git merge-tree` trial reported no conflict markers. Nothing has been reconciled.

The published image is `scrumpis/dndsr` `v0.2.0`, digest in `skills/dndsr/SKILL.md`. It holds dependencies only.

## Next

1. Maintainer pushes `dev`.
2. On instruction, run the scientific and software audits.
3. Build a small validated example dataset and the first tests. None exist.
4. Prepare the merge and v0.3.0 once authorized.
