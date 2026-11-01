# Release readiness: v0.3.0 candidate

**Recommendation: NO-GO until CI passes on the commits after `b59cf23`.** The code passes every check available locally. What is missing is CI on the fixed workflows, and maintainer settings. Unpublished real-genome results stay outside Git.

## Blocking the release

1. **Push `dev`, and get CI green,** including `R-CMD-check` on macOS, Windows, and Linux, which has never reached package checks, and the Docker and Apptainer container steps.
2. **Rulesets and Docker Hub setup:** `GITHUB-SETUP.md`.
3. **Final approval** to merge, tag, and publish.

## CI at `b59cf23`, the only push so far

| Workflow | Result | Note |
|---|---|---|
| `skill` | Pass | |
| `install` (clean conda, example) | Pass | |
| `container`: real `docker build` | **Pass** | |
| `container`: example in the image | Fail | Very likely a folder-permission bug in the workflow (logs need a login). Fixed, not yet confirmed |
| `R-CMD-check`, all platforms | Fail before any package check | `upgrade: true` rejected. Fixed, not yet confirmed |

## Local evidence at the final code

| Check | Result |
|---|---|
| `R CMD check --no-manual`, published image, `lme4` absent | Status OK |
| Package tests | 56 tests, 264 expectations, 0 failed, 0 skipped |
| Tutorial, clean `HOME`, Singularity | Pass |
| R tutorial, every block | Pass |
| Candidate image from the Dockerfile's steps, via Apptainer, `container-checks.sh` | All pass, on R 4.6.1 |
| Clean conda, CLI and R | Pass at `fedd59e` locally, and on GitHub at `b59cf23` |
| Rootless Apptainer tutorial | Pass at an earlier commit |
| Fresh Claude agent at `c1eda04` | Correct; results recheck 14 of 14 |
| Git safeguards | 88 of 88 |
| Public real-genome pair | `docs/validation/public-genome.md` |
| Regression coverage | `docs/validation/regression-coverage.md` |

## Not done or not tested

- macOS and Windows: `R CMD check` will run in CI, without DIAMOND, so dN/dS tests skip. Runtime use is untested.
- R from CRAN and Bioconductor failed on the tested HPC module, because of `rtracklayer` and OpenSSL. Conda is the recommended route.
- **Codex:** partial, and not rerun at the final commit. In its sandbox, Singularity could not start, so it ran the image under Bubblewrap.
- Docker Hub push: cannot be tested without publishing.
- Enrichment on a public annotated genome: constructed cases only. See `ROADMAP-1.0.md`.
- Open items: `docs/validation/known-issues.md`.
