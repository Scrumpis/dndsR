# Release readiness

## v0.5.0 candidate: `dev` at `010ab1e` and its records commit

**Assessed 2026-11-20. Ready to merge once the two items under "Blocking" are done; the maintainer approved the scope, the push of `dev`, the clone references, and the agent test.** Not merged, tagged, or published. Release notes draft: [release-notes-0.5.0.md](release-notes-0.5.0.md).

### Checklist

| Step | State |
|---|---|
| Version 0.5.0 in `DESCRIPTION`, skill, `CITATION.cff`; NEWS | Done (`0d00046`) |
| `p_best_fits_adj`, independently tested; docs reconciled with the calibration | Done (`0d00046`) |
| Clone references `--branch v0.5.0` (RELEASE.md step 5) | Done (`1ea7faf`); they fail until the tag exists, as approved |
| Interim image documented as v0.4.0, without HyPhy and KaKs_Calculator | Done (`010ab1e`); the digest is replaced after publication (publish step 3) |
| Fixture kept, attributed, outside dndsR's license | Done (`1ea7faf`, `tools/dev/fixtures/atal/NOTICE.md`) |
| Second author's metadata | Kept provisionally, no ORCID (maintainer) |
| Local gate: `check.sh --rcmd-check --regression` | Pass before `0d00046` (236 tests, 0 failed, 4 skipped for HyPhy and KaKs_Calculator); `check.sh` with the tutorial at `010ab1e`'s tree: pass |
| Safeguards | 97 of 97 after the version bump |
| HyPhy tests and the recommended settings, with HyPhy | Pass (3 of 3; independent BH check) |
| New container checks (KaKs_Calculator through `--dnds-method YN`; `phylo_busted` with the recommended settings) | Pass against this tree on Slurm, with the tools from Bioconda |
| Agent test (`test-agent.sh claude`) | 7 of 7 at `010ab1e` |
| Whole-genome regression | Not needed: no pairwise code changed since the rerun of the gene-ID fix (`f529235`), and the fixture is byte-identical |
| **CI on the exact candidate**, including the container workflow with `DNDSR_FAIL_ON_SKIP=1` | **Pending the push.** Read its annotations after the push: `Package tests` notice with `skipped=0`, no `Skipped test` warning |
| MK `gaps.md` correction | **Pending the MK maintainer** (channel message 11) |
| Merge `dev` to `main`, tag, image, digests, GitHub release, Zenodo | Waiting for approval (`RELEASE.md`, "Publish, once authorized") |

### Blocking

1. CI green on the pushed candidate, with the container workflow showing that no test skipped in the freshly built image.
2. The MK maintainer's correction of the MK line in `skills/dndsr/references/gaps.md`; then rerun CI on that commit.

### Evidence

| Check | Run | Result |
|---|---|---|
| `check.sh --rcmd-check --regression`, Slurm, 8 CPUs, 48 GB, v0.4.0 image | This tree before commit | Tests 236, 0 failed, 0 errors, 4 skipped (HyPhy, KaKs_Calculator absent from the image); R CMD check ok; tutorial ok; skill ok; fixture regression byte-identical. Safeguards: 3 allowed pushes blocked, because the test clones the last commit (version 0.4.0.9000) and copies the working tree's skill (0.5.0). Rerun after committing the bump (`0d00046`, Slurm): 97 passed, 0 failed |
| BUSTED tests with HyPhy 2.5.103 (`test-phylo-busted.R`, the 3 tests the gate skips) | Slurm, this tree | 3 tests, 26 expectations, 0 failed, none skipped |
| Recommended settings on the example, with HyPhy (`--multiple-hits Double --extra-seeds 2,3`) | Slurm, this tree | 8 tested, 3 starts; `p_adj` and `p_best_fits_adj` equal an independent BH (6.7e-16) |
| CI at `3943ffc` (pushed) | GitHub | All five workflows pass, `regression` included; the container workflow's fixture step passes in the freshly built image |
| `R-CMD-check`: macOS, Windows, Ubuntu release, oldrel, devel | CI at `b65bf8b` | Pass. Windows failed at `3a62ad2` and `60339ce` on a test expectation (path separators); fixed in the test, `be1a7fa` |
| `install`, `skill`, `container` | CI at `b65bf8b` | Pass. `container` failed once at `3a62ad2`, in the documented-installer step, and passed at the next two commits with the same code; logs need a login, so the cause is unconfirmed |
| `R CMD check --as-cran` on the tarball, local development image (HyPhy, KaKs_Calculator) | Local at `3a62ad2` | 0 errors, 0 warnings, 2 notes (new submission; no pandoc) |
| `tools/dev/check.sh` | Local at `b65bf8b`; again after `fda5dbe` with `--rcmd-check` | At `b65bf8b`: tests 195, 0 failed; safeguards 1 failure (stale installed `pre-push` copy). After `fda5dbe`: tests 196, 0 failed; R CMD check `Status: OK`; tutorial OK; after `setup.sh`, safeguards 97 passed, 0 failed |
| v0.4.0 against `dev`, real genomes, data read-only, the same v0.4.0 image for both (`tools/dev/validation/`, results outside Git) | Public *A. thaliana*–*A. lyrata* pair and two unpublished pairs, at `b65bf8b` | Every CDS FASTA, CDS QC table, dN/dS table, and annotated table **byte-identical** to v0.4.0. Differences: the new `_dnds_reliability.tsv`; provenance records `dnds_method` and the development version; logs differ in times and process IDs. On the public pair, `compare-versions.R` with gffread 0.12.9: every CDS compared matches gffread (0 differ, 83 not comparable by name, as before); the same 22,149 pairs with the same partners; all 22,094 comparable dN/dS values identical. Details: [public-genome.md](../validation/public-genome.md#v040-against-dev) |
| v0.4.0 against `dev`, real genomes again, after PRs #9 and #10 (data read-only, same v0.4.0 image, the v0.4.0 outputs of the row above reused as the baseline; its provenance confirmed: run log and provenance say `957cfe1` and dndsR 0.4.0, inputs match their recorded SHA-256) | Public pair and two unpublished pairs, at `9a5c758` | Every CDS FASTA, CDS QC table, dN/dS table, and annotated table **byte-identical** to v0.4.0 on all three pairs (7 and 13 files). Differences: the reliability tables v0.4.0 does not write; provenance (version, method); logs (times, process IDs). Checksums of both sets are kept beside the results, outside Git. `compare-versions.R` against gffread was not rerun: the CDS are byte-identical to those it checked at `b65bf8b` |
| v0.4.0 against the gene-ID fix (tree snapshot `2b7385b`, identical in `R/`, `DESCRIPTION`, `NAMESPACE`, `man/`, and `tests/` to the commit that adds this row), whole genomes, same image and v0.4.0 baseline | Public pair and two unpublished pairs | Every CDS FASTA, dN/dS table, and annotated table byte-identical; the unpublished pairs' QC tables too. Public pair: `Alyrata_CDS_qc.tsv` differs only in the `transcript` column of 84 chloroplast genes whose RefSeq records lack `transcript_id`, now named by ID instead of a Bioconductor fallback (`B9N47_pgp001.tx53` to `B9N47_pgp001`). No `append_annotations` run reported missing IDs |
| Real-genome fixture regression (`tools/dev/regression.sh`, new, [TESTING.md](TESTING.md)) | v0.4.0 image, `9a5c758` plus this session's change, and v0.4.0 itself | Baseline made with the current tree; v0.4.0 reproduces every output it writes (508 pairs, 0 values changed). Fixture pairs equal the full-genome run's for all 503 shared pairs: `tools/dev/fixtures/README.md` |
| Phylogenetic calibration, stages 2A and 2B | Simulation, `3603efc` | Reported: `docs/validation/phylo-calibration-stage2b.md`; guidance approved by the maintainer |

### Not blocking

- MK: message 9's documentation of `--unique-pairs` with MK rows, and message 7's SNP-file pre-check (MK maintainer).
- The `mamba-org/setup-micromamba@v2` action runs on a deprecated Node.js version (CI warning only).

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
