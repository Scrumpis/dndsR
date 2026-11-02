# Project status

Update this at the end of each substantive session. Identify things by commit and hash, not by date.

**Phase:** v0.3.0 is public: tag `d5f4b21`, image `sha256:feb5e4b1`, DOI 10.5281/zenodo.23879469. `main` is `5ef6fa7`. `dev` holds the v0.4.0 candidate, software at `5983405`, which needs a merge and a tag. Checkpoints are in [PROGRESS.md](PROGRESS.md).

## Identity

dndsR is a real R package, fit for CRAN one day, and a containerized CLI, for people and AI agents equally, with ease of use first.

## Refs

| Ref | Commit | Notes |
|---|---|---|
| `v0.2.0` | `b8ea5ff` | Prerelease |
| `v0.3.0` | `d5f4b21` | Release |
| `main` | `5ef6fa7` | `v0.3.0`, and the DOI badge and correctness notice from PR #3 |
| `dev` | See `git log` | `main` was merged into it at `3f86d42`, by the maintainer |

## Next

1. All workflows passed at `5983405`: see [RELEASE-READINESS.md](RELEASE-READINESS.md). The maintainer pushes the records commit after it.
2. The maintainer opens the pull request to `main`, merges, and tags `v0.4.0`. Then `publish.yaml` asks for approval. Follow `docs/dev/RELEASE.md`; notes are in `docs/dev/release-notes-0.4.0.md`.
3. After that, this software is frozen for the Chenopodium analysis. Its work and results live outside this repository, and nothing from them comes here.
4. For v1.0.0: `docs/dev/ROADMAP-1.0.md`. Do not start it for the current manuscript.
