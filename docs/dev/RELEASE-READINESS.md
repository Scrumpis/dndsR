# Release readiness: v0.3.0 candidate

**Recommendation: GO on the evidence, subject to the maintainer's setup and approval.** No code-side blocker remains. Unpublished real-genome results stay outside Git.

## Blocking the release

Only the maintainer can do these:

1. **Rulesets and Docker Hub setup:** `GITHUB-SETUP.md`. As read from the public API at `f84521f`, `protect_main` still requires 1 approval, and no `dev` or tag ruleset is visible.
2. **Authorization** for the one "Update branch" merge of `main` into `dev`.
3. **Final approval** to merge, tag, and publish.

## CI at `f84521f`

| Workflow | Result |
|---|---|
| `R-CMD-check`: macOS, Windows, Ubuntu release, oldrel, devel, and the `R-CMD-check` rollup job | **Pass.** DIAMOND is absent there, so dN/dS tests skip |
| `install`: clean conda, CLI, and example | **Pass** |
| `skill` | **Pass** |
| `container`: real `docker build` from a clean checkout | **Pass** |
| `container`: full container checks in the Docker image | **Pass** |
| `container`: R tutorial in the Docker image | **Pass** |
| `container`: package tests in the Docker image | **Pass** |
| `container`: Apptainer 1.5.4 from conda, checks on the image | **Pass** |
| `container`: documented installer with Apptainer | **Pass** |

Step and job results are from the public API. Job logs need a login and were not read.

`R/`, `inst/`, `skills/`, `tools/Dockerfile`, `tools/dndsr-install`, `DESCRIPTION`, and `NAMESPACE` have no differences between `c1eda04` and `f84521f`, so the agent test and the local checks below still describe this code. `R/`, `inst/`, `tools/Dockerfile`, `DESCRIPTION`, and `NAMESPACE` also have no differences from `59faf30`, the commit used for the public-genome comparison.

## Local evidence at the final code

| Check | Result |
|---|---|
| `R CMD check --no-manual`, published image, `lme4` absent | Status OK |
| Package tests | 56 tests, 264 expectations, 0 failed, 0 skipped |
| Tutorial, clean `HOME`, Singularity | Pass |
| R tutorial, every block | Pass |
| Candidate image from the Dockerfile's steps, via Apptainer, `container-checks.sh` | All pass, on R 4.6.1 |
| Package tests in that image | 56 tests, 264 expectations, 0 failed, 0 skipped |
| Clean conda, CLI and R | Pass at `fedd59e` locally, and on GitHub at `b59cf23` |
| Rootless Apptainer tutorial | Pass at an earlier commit |
| Fresh Claude agent at `c1eda04` | Correct; results recheck 14 of 14 |
| Git safeguards | 88 of 88 |
| Public real-genome pair | `docs/validation/public-genome.md` |
| Regression coverage | `docs/validation/regression-coverage.md` |

## Not done or not tested

- macOS and Windows: `R CMD check` passes in CI, without DIAMOND, so dN/dS tests skip. Full analyses have not been run there.
- R from CRAN and Bioconductor failed on the tested HPC module, because of `rtracklayer` and OpenSSL. Conda is the recommended route.
- **Codex:** partial, and not rerun at the final commit. In its sandbox, Singularity could not start, so it ran the image under Bubblewrap.
- Docker Hub push: cannot be tested without publishing.
- Enrichment on a public annotated genome: constructed cases only. See `ROADMAP-1.0.md`.
- Open items: `docs/validation/known-issues.md`.
