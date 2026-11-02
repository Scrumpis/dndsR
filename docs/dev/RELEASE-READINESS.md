# Release readiness: v0.3.0 candidate

**Local candidate ready for the maintainer's review, at the `dev` commit that adds this text. The software is at `1a363fb`.** It has not been through CI: the CI table below is from `f84521f`, before the enrichment correction. The maintainer's earlier GO at `f84521f` predates these changes, and final approval is still needed to merge, tag, and publish. Unpublished real-genome results stay outside Git.

**Corrected since `f84521f`, with maintainer approval.** Multiple-testing correction in `term_enrichment` and `ipr_enrichment`; `#` in GFF attributes; the `go_enrichment(include_definition = FALSE)` crash. See [enrichment calibration](../validation/enrichment.md).

**Limitations to read before approving.** Enrichment is corrected for multiple testing, not for gene length, dS, or annotation completeness. Default `ipr_enrichment` corrects each InterPro type separately, so across tables 9–17% of null runs have a significant term, documented and left as the default. No enrichment result has been compared with independent biological evidence.

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

## Local evidence at `1a363fb`, the corrected code

Later commits change only documentation and validation scripts: `R/`, `inst/`, `DESCRIPTION`, and `NAMESPACE` are unchanged.

| Check | Result |
|---|---|
| `tools/dev/check.sh` at final HEAD: safeguards, skill, tests, tutorial | Pass. 78 tests, 333 expectations, 0 failed, 0 skipped |
| `R CMD check --no-manual`, published image | Status OK. 328 passed, 0 failed, 5 expected warnings from tests that exercise no network and a missing file |
| Candidate image from the Dockerfile's steps, via Apptainer, `container-checks.sh` | All pass |
| Package tests in that image | 78 tests, 333 expectations, 0 failed, 0 skipped |
| R tutorial in that image | 4 blocks ran |
| Full CLI workflow on the public pair, against the previous candidate | `AtAl_dnds.tsv`, summary, and ideogram SVGs byte-identical. Annotated table differs on the 73 recovered rows only |
| Null models and independent recalculation | [enrichment.md](../validation/enrichment.md) |
| Fresh Claude agent at `cf4e61a`, `tools/dev/test-agent.sh claude` | 7 of 7 checks; its results pass an independent `verify_example`, 14 of 14 |
| Thread determinism on the real table | Annotated table and all `term_enrichment`, `ipr_enrichment`, and `go_enrichment` tables byte-identical at 1 and 8 threads |

## Local evidence at earlier code

These were not rerun after the correction, unless also listed above.

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
- Biological comparison of enrichment terms, and enrichment at other divergences: open. See `ROADMAP-1.0.md`.
- Open items: `docs/validation/known-issues.md`.
