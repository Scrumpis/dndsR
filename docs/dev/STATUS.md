# Project status

Update this at the end of each substantive session. Identify things by commit and hash, not by date.

**Phase:** v0.3.0 candidate on `dev`. Software is at `eed1096`. All CI workflows passed at `f2da1f1`. The maintainer accepted GO at 9.5% under the `strata` null as a documented limitation. Not merged, tagged, or published. Read [RELEASE-READINESS.md](RELEASE-READINESS.md) first. Checkpoints are in [PROGRESS.md](PROGRESS.md).

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

1. The maintainer completes `GITHUB-SETUP.md`, then opens the pull request and follows `RELEASE.md`. Agents do not push; this machine has no GitHub credentials.
2. After v0.3.0, continue `ROADMAP-1.0.md`, starting with gene-length and dS confounding, independent estimator benchmarking, and empirical calibration of enrichment statistics. Threading comes after that: [performance.md](../validation/performance.md).
