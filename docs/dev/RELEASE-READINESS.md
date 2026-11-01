# Release readiness: v0.3.0 candidate

**Recommendation: NO-GO until the Apptainer steps pass in CI, then ready for the maintainer's decisions.** R checks and the real Docker image pass on GitHub. Unpublished real-genome results stay outside Git.

## Blocking the release

1. **Push `dev`, and confirm the `container` workflow,** whose Apptainer install step failed at `fca3306`, so its Apptainer tests and the new test step have not run on GitHub.
2. **Rulesets and Docker Hub setup:** `GITHUB-SETUP.md`.
3. **Final approval** to merge, tag, and publish.

## CI at `fca3306`

| Workflow | Result |
|---|---|
| `R-CMD-check`: macOS, Windows, Ubuntu release, oldrel, devel, and the `R-CMD-check` rollup job | **Pass.** DIAMOND is absent there, so dN/dS tests skip |
| `install`: clean conda, CLI, and example | **Pass** |
| `skill` | **Pass** |
| `container`: real `docker build` from a clean checkout | **Pass** |
| `container`: full container checks in the Docker image | **Pass** |
| `container`: R tutorial in the Docker image | **Pass** |
| `container`: installing Apptainer | Fail: third-party action could not install on the runner. Replaced with conda Apptainer. Not yet confirmed |
| `container`: Apptainer checks, documented installer | Not run on GitHub. Tested locally |

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
