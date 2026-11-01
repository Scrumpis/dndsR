# Project status

Update this at the end of each substantive session. Identify things by commit and hash, not by date.

**Phase:** v0.3.0 candidate on `dev`. All CI passes at `f84521f`, and the release awaits the maintainer's setup and approval. Read [RELEASE-READINESS.md](RELEASE-READINESS.md) first. Checkpoints are in [PROGRESS.md](PROGRESS.md).

## Identity

dndsR is a real R package, fit for CRAN one day, and a containerized CLI, for people and AI agents equally, with ease of use first. The path is: audit `dev`, merge to `main`, release v0.3.0, then v1.0.0 for publication.

## Baseline

| Ref | Commit | Notes |
|---|---|---|
| `v0.2.0` | `b8ea5ff` | Prerelease on GitHub, version 0.2.0.000 |
| `main` | `bbd8384` | 8 commits past `v0.2.0`, same content as merge base `9094341` |
| `dev` at start of work | `c8b534f` | 134 commits past `v0.2.0` |

A read-only trial merge of the original `dev` into `main` had no conflicts. `main` has not been touched.

## Next

1. The maintainer completes `GITHUB-SETUP.md`, then follows `RELEASE.md`. Agents do not push; this machine has no GitHub credentials.
2. After release, work through `ROADMAP-1.0.md`. Its first decision is item 5, per-gene uncertainty.
