# Release readiness: v0.3.0 candidate

**Local candidate, with one statistical result for the maintainer to accept or reject before it is called ready.** Software is at `eed1096`. It has not been through CI: the CI table below is from `f84521f`, before the enrichment changes. The maintainer's earlier GO at `f84521f` predates these changes, and final approval is still needed to merge, tag, and publish. Unpublished real-genome results stay outside Git.

**Corrected since `f84521f`, with maintainer approval.** Multiple-testing correction in `term_enrichment` and `ipr_enrichment`; `#` in GFF attributes; the `go_enrichment(include_definition = FALSE)` crash. **Defaults changed at `eed1096`, as the maintainer directed:** `ipr_enrichment` corrects across all entry types, and `go_enrichment` applies BH across all requested ontologies. See [enrichment calibration](../validation/enrichment.md).

**The open result.** Default `go_enrichment`, with no true association, has a false term in 1.0% and 0.8% of runs under `permute` and `rotate`, and in **9.5% (7.0–12.8) under `strata`**, above nominal. GO:0005576, extracellular region, is in 28 of those 38 runs; without it the rate is 3.0% (1.7–5.2). The cause is confounding with length and dS, not the correction. No correction was added, as directed. Recommended: accept as a documented limitation, with modeling on the v1.0.0 roadmap.

**Other limitations.** Enrichment is corrected for multiple testing, not for gene length, dS, or annotation completeness. False discovery with real associations has not been measured. `term_enrichment` corrects each term type separately. No enrichment result has been compared with independent biological evidence.

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

## Local evidence at `eed1096`, the final defaults

Later commits change only records. `R/`, `inst/`, `DESCRIPTION`, and `NAMESPACE` are unchanged.

| Check | Result |
|---|---|
| `tools/dev/check.sh`: safeguards, skill, tests, tutorial | Pass. 85 tests (the earlier 78 and 7 for these defaults), 382 expectations, 0 failed, 0 skipped |
| `R CMD check --no-manual`, published image | Status OK, no notes or warnings from the check. Tests: 377 passed, 0 failed, 5 warnings, below |
| Candidate image from the Dockerfile's steps, via Apptainer, `container-checks.sh` | All pass |
| Package tests in that image | 85 tests, 382 expectations, 0 failed, 0 skipped |
| R tutorial in that image | 4 blocks ran |
| Full CLI workflow on the public pair, against `1a363fb` | `AtAl_dnds.tsv`, annotated table, summary, ideogram SVGs, and CDS FASTA byte-identical |
| `--adjust-scope per_type`, on real data | Tables equal `1a363fb`'s on every shared column, with parent-child too. Null runs on shared seeds reproduce exactly |
| Null models and independent recalculation | [enrichment.md](../validation/enrichment.md) |
| Thread determinism on the real table | Annotated table and all `ipr_enrichment` and `go_enrichment` tables byte-identical at 1 and 8 threads |
| Fresh Claude agent at `eed1096`, `tools/dev/test-agent.sh claude` | 7 of 7 checks; its results pass an independent `verify_example`, 14 of 14 |

The five test warnings, each confirmed in `tests/testthat.Rout`. All come from conditions these tests create, and none affects a result:

| Test | Warning | Cause |
|---|---|---|
| `test-enrichment-analytic.R:53` | `release_notes.txt`: `Couldn't connect to server` | `download.file`. The test sets a refusing proxy on purpose, to show that `ipr_enrichment` works offline |
| `test-enrichment-analytic.R:53` | `entry.list`: `Couldn't connect to server` | The same, for the second download |
| `test-example-pipeline.R:158` | `path[1]="no_q.fasta": No such file or directory` | `normalizePath(qs)` in `.resolve_inputs`, in `calculate_dnds`, before the "not found" error the test expects |
| `test-example-pipeline.R:158` | `path[1]="no_s.fasta": No such file or directory` | The same, `normalizePath(ss)` |
| `test-statistics.R:125` | `pos_threshold (1) is < neutral_upper (1.1)...` | `dnds_state_contrast`'s own defaults overlap, so positive overrides neutral above 1, and it warns about this at its defaults. Behavior is documented; not changed |

## Local evidence at `1a363fb`

Superseded above where the same check appears. Evidence: [enrichment-1a363fb.md](../validation/enrichment-1a363fb.md).

| Check | Result |
|---|---|
| `tools/dev/check.sh` | Pass. 78 tests, 333 expectations, 0 failed, 0 skipped |
| `R CMD check --no-manual`, published image | Status OK. 328 passed, 0 failed, the same 5 warnings |
| Candidate image, package tests, R tutorial | All pass |
| Fresh Claude agent at `cf4e61a` | 7 of 7 checks; 14 of 14 |

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
