# Release readiness

## Next release candidate: `dev` at `b65bf8b` and later (not yet a candidate)

**Assessed 2026-11-14 (US Eastern). Not ready for release; blockers below.** 58 commits since v0.4.0. Version `0.4.0.9000`; no release number chosen.

| Check | Run | Result |
|---|---|---|
| `R-CMD-check`: macOS, Windows, Ubuntu release, oldrel, devel | CI at `b65bf8b` | Pass. Windows failed at `3a62ad2` and `60339ce` on a test expectation (path separators); fixed in the test, `be1a7fa` |
| `install`, `skill`, `container` | CI at `b65bf8b` | Pass. `container` failed once at `3a62ad2`, in the documented-installer step, and passed at the next two commits with the same code; logs need a login, so the cause is unconfirmed |
| `R CMD check --as-cran` on the tarball, local development image (HyPhy, KaKs_Calculator) | Local at `3a62ad2` | 0 errors, 0 warnings, 2 notes (new submission; no pandoc) |
| `tools/dev/check.sh` | Local at `b65bf8b` | Tests 195, 0 failed; imports, skill ok. Safeguards: 1 failure, this clone's installed `pre-push` hook predates PR #8's change: the maintainer reruns `tools/dev/setup.sh` |
| v0.4.0 against `dev`, real genomes, data read-only, the same v0.4.0 image for both (`tools/dev/validation/`, results outside Git) | Public *A. thaliana*–*A. lyrata* pair and two unpublished pairs, at `b65bf8b` | Every CDS FASTA, CDS QC table, dN/dS table, and annotated table **byte-identical** to v0.4.0. Differences: the new `_dnds_reliability.tsv`; provenance records `dnds_method` and the development version; logs differ in times and process IDs. On the public pair, `compare-versions.R` with gffread 0.12.9: every CDS compared matches gffread (0 differ, 83 not comparable by name, as before); the same 22,149 pairs with the same partners; all 22,094 comparable dN/dS values identical. Details: [public-genome.md](../validation/public-genome.md#v040-against-dev) |
| Phylogenetic calibration, stages 2A and 2B | Simulation, `3603efc` | Reported: `docs/validation/phylo-calibration-stage2b.md`; guidance approved by the maintainer |

**What changed for users since v0.4.0** (NEWS): the experimental phylogenetic commands; `--dnds-method` alternatives through KaKs_Calculator; the reliability table; the development McDonald-Kreitman test; header rows in comparison files now work. Pairwise results with default settings are unchanged, as shown above.

### Blocking the release

1. **McDonald-Kreitman: keep experimental, or hold back.** Its method and output columns changed on 2026-11-14 (PR #8: allele-aware classes, BED needs ALT in column 4, `snps_midfold` and `div_steps_midfold` removed). It is tested on designed examples only, not on real population data. Real-data testing and the maintainer's approval come before it is released as anything but experimental. Open integration point: SNP-file errors appear only after dN/dS (channel message 7).
2. **Phylogenetic module: still experimental.** The calibration covers one five-taxon tree in simulation. No published image includes HyPhy, and the release that carries the module is an open decision. `p_adj` corrects `p_value`, not the recommended `p_best_fits`; whether to add a corrected column is open (`DECISIONS.md`).
3. **Author metadata.** Confirm the second author's preferred name and optional ORCID before publication (`CITATION.cff`, `.zenodo.json`, `DESCRIPTION`).
4. **Maintainer steps.** Rerun `tools/dev/setup.sh`; choose the version number; then `RELEASE.md` (merge, tag, image, digests in the docs). Each needs the maintainer's approval.

### Not blocking, but open

- `check.sh` does not run `R CMD check`; CI does, and the import check covers the failure that reached CI on 2026-11-14.
- The `pre-push` hook now allows `feat/*` branches (PR #8); no safeguard test covers that path yet.


## v0.4.0 candidate

**Software `5983405`. All four workflows passed there.** Later commits change records only.

| Check | Run | Result |
|---|---|---|
| `R-CMD-check`: macOS, Windows, Ubuntu release, oldrel, devel | CI at `5983405` | Pass |
| `install`, `skill` | CI at `5983405` | Pass |
| `container`: `docker build --no-cache` and container checks | CI at `5983405` | Pass |
| `R CMD check` on the tarball, v0.3.0 image | Local at `5983405` | `Status: OK`, no warnings or notes |
| `tools/dev/check.sh`: 102 tests, 468 expectations, safeguards, skill, tutorial | Local at `5983405`, v0.3.0 image | Pass, 0 skipped |
| Tutorial with Apptainer | Local at `4e2161d`, same software but for notice wording | Pass |
| `tools/dev/test-agent.sh claude`, fresh agent | Local at `4e2161d` | 7 of 7 |
| Public-pair pipeline, `tools/dev/validation` | v0.3.0 image | Thread counts identical; differences from the older stack explained: `docs/validation/stack-comparison.md` |
| Independent reproduction of `covariate_enrichment` | 1,200 tables, 2 thresholds | 4,800 of 4,800 comparisons agree: `docs/validation/enrichment-confounding.md` |
| Real genomes, 9 comparisons, private | At `5983405` | All steps ran; CDS identical to gffread; thread counts and reruns byte-identical. Results are outside this repository |
| `tools/dev/test-public-tutorial.sb`, public genomes, clean HOME | Local at `5983405` | Pass, on the third try: UniProt's stream endpoint returned HTTP 500 for about four hours before that, which the first two tries stopped on |

**Scientific correctness:** no known material defect is open. Limits accepted and documented: `covariate_enrichment` is calibrated in simulation at `dNdS > 1` on 4 plant pairs, and is not validated for detecting real selection; `q0.95` exceeded 5% in one setting; `go_enrichment` and `term_enrichment` keep their defaults and warn; Comeron stays the default with its known bias; reciprocal best hits can mispair paralogs.

**Not tested:** Docker on a user's machine, beyond CI; macOS beyond CI.

**Maintainer steps:** `RELEASE.md`. Before merging, this commit's CI must pass.

## v0.3.0 candidate

**Candidate passed CI at `f2da1f1`, with software at `eed1096`.** Later commits change records and documentation only, and need their own green CI before the pull request. Final approval is still needed to merge, tag, and publish. Unpublished real-genome results stay outside Git.

**Corrected since `f84521f`, with maintainer approval.** Multiple-testing correction in `term_enrichment` and `ipr_enrichment`; `#` in GFF attributes; the `go_enrichment(include_definition = FALSE)` crash. **Defaults changed at `eed1096`, as the maintainer directed:** `ipr_enrichment` corrects across all entry types, and `go_enrichment` applies BH across all requested ontologies. See [enrichment calibration](../validation/enrichment.md).

**Accepted limitation.** Default `go_enrichment`, with no true association, has a false term in 1.0% and 0.8% of runs under `permute` and `rotate`, and in **9.5% (7.0–12.8) under `strata`**, above nominal. GO:0005576, extracellular region, is in 28 of those 38 runs; without it the rate is 3.0% (1.7–5.2). The cause is confounding with length and dS, not the correction. The maintainer accepted this as a documented limitation for v0.3.0: no terms removed, no thresholds changed, no further correction. BH does not guarantee nominal error control under all conditions, and significant terms remain exploratory.

**Other limitations.** Enrichment is corrected for multiple testing, not for gene length, dS, or annotation completeness. False discovery with real associations has not been measured. `term_enrichment` corrects each term type separately. No enrichment result has been compared with independent biological evidence.

## Blocking the release

Only the maintainer can do these. As read from the public API at `f2da1f1`:

1. **`protect_main`** still requires 1 approval, which an author cannot give, and only the `R-CMD-check` check. No `protect_tags` or `protect_dev` ruleset exists. See `GITHUB-SETUP.md`.
2. **Authorization** for the one "Update branch" merge of `main` into `dev`, which strict status checks require: `main` is 2 commits ahead.
3. **Docker Hub setup.** No `dockerhub-release` environment exists. Whether `DOCKERHUB_PUBLISH_ENABLED` is set cannot be read without a login. `scrumpis/dndsr:v0.3.0` does not exist yet, as the workflow requires.
4. **Final approval** to merge, tag, and publish.

## CI at `f2da1f1`

Run on push to `dev`; `main` unchanged at `bbd8384`.

| Workflow | Result |
|---|---|
| `R-CMD-check`: macOS, Windows, Ubuntu release, oldrel, devel, and the `R-CMD-check` rollup job | **Pass.** DIAMOND is absent there, so dN/dS tests skip |
| `install`: clean conda, CLI, and example | **Pass** |
| `skill` | **Pass** |
| `container`: real `docker build --no-cache` from a clean checkout | **Pass** |
| `container`: full container checks in the Docker image | **Pass** |
| `container`: R tutorial in the Docker image | **Pass** |
| `container`: package tests in the Docker image, which include the scientific regression tests with DIAMOND | **Pass** |
| `container`: Apptainer from conda, checks on the image | **Pass** |
| `container`: documented installer with Apptainer | **Pass** |

Step and job results are from the public API. Job logs need a login and were not read, so test counts in CI are not confirmed. There is no separate scientific-regression workflow. The null models and real-genome comparisons are not in CI: they ran locally, as recorded below. `container` failed on its build step at `5d53ced`, about 90 seconds in, and its logs were not read, so the cause is unknown. That build passed here with the same Dockerfile: `tools/Dockerfile` has no differences between those commits.

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
