# Project status

Update this at the end of each substantive session. Identify things by commit and hash, not by date.

**Phase:** v0.3.0 is public: tag and `main` at `d5f4b21`. `dev` holds v1.0.0 validation work, unpushed until the maintainer pushes it. `covariate_enrichment` awaits approval. Checkpoints are in [PROGRESS.md](PROGRESS.md).

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

Open maintainer decisions are in [DECISIONS.md](DECISIONS.md).

1. The maintainer pushes `dev`, and CI runs on these commits. It has not run on them yet.
2. The maintainer decides on `covariate_enrichment` and on a tagged release. `docs/dev/release-notes-next.md` is a draft.
3. A publication analysis follows `docs/publication-protocol.md`, from that tag. Do not put unpublished findings in this repository.
4. For v1.0.0: tests chosen in advance against independent evidence, isoform effects on pairing, removing gene-level `dnds_state_contrast`, and threading in `performance.md`.
