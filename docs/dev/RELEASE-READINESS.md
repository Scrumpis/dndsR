# Release readiness: v0.3.0 candidate

**Recommendation: NO-GO today, with a short path to GO.** The software evidence is strong, but one required gate could not run: nothing has been pushed, because this machine has no GitHub credentials. So CI, a real `docker build`, and macOS and Windows checks have never happened. The code is not the blocker; the remaining work is running what is already written.

## Evidence at the final code

Commit `fedd59e` for code. Later commits change only documentation.

| Check | Result |
|---|---|
| `R CMD check --no-manual`, in the published image, `lme4` absent | Status OK: no errors, warnings, or notes |
| Package tests | 53 tests, 256 expectations, 0 failed, 0 skipped |
| `dndsr verify_example` | 14 of 14 |
| Tutorial commands, clean `HOME`, Singularity | Pass (`fedd59e`) |
| Tutorial commands, clean `HOME`, rootless Apptainer 1.5.4 | Pass at an earlier commit (`7ae5a6a` or before); installer and tutorial unchanged in behavior since |
| Candidate image, built from the Dockerfile's `RUN` steps through Apptainer | Built and passes 14 of 14; dN and dS match the reference exactly, on R 4.6.1 |
| Clean conda install | Pass, 14 of 14, CLI and R |
| Fresh Claude agent, clean clone, plain request | 7 of 7 script-checked criteria |
| Fresh sandboxed Codex agent | Partial; see below |
| Git safeguards | 88 of 88 |
| Real genomes | `docs/validation/real-genomes.md` |
| Estimators | `docs/validation/estimators.md` |

## Blocking the release

| Gate | Status | Who |
|---|---|---|
| Push `dev`, then CI green on all workflows | **Not run. No GitHub credentials on this machine.** Expect some first-run failures | Maintainer provides credentials or pushes; then fix findings |
| `docker build` of `tools/Dockerfile` | **Not run.** Its steps ran through Apptainer, but Docker has not parsed the file or applied `.dockerignore`. `container.yaml` does this | CI |
| Regenerate affected Chenopodium analyses | **Needs a decision.** Project results use the old code | Maintainer |
| Publishing decisions below | Waiting | Maintainer |

## Decisions for the maintainer

1. **Authorize the v0.1.0 and v0.2.0 warnings.** Drafts: `release-notice-draft.md`. Recommendation: post them as soon as possible, since those releases are public and the effects are large.
2. **Review your own results.** Figures and enrichment built on `Chenopod_dnds` used affected code.
3. **Default estimator.** Comeron stays for v0.3.0, as directed. It read low in simulation, and agreed well with other methods on real genes. Recommendation: add KaKs_Calculator to the image and decide before v1.0.0.
4. **Gene-level `dnds_state_contrast`.** Confirmed by running that it gives no usable result. It stays experimental. Recommendation: `state-contrast-design.md`, option 1 or 3, before v1.0.0.
5. **`main` ruleset.** Whether the new `R-CMD-check` job satisfies it is unconfirmed until a pull request runs, and the required approval cannot come from a pull request's author.
6. **Is the aggregate real-genome summary fit to publish?** `docs/validation/real-genomes/` has per-comparison counts from *C. formosanum* and *C. album* data. No sequences or gene IDs are included. Please confirm before pushing if those genomes are unpublished.
7. **Agent pushes.** You authorized pushing to `dev`. The project guard still denies it, because agents cannot push here anyway. To allow it, remove the push rules from `.claude/settings.json` and `.claude/hooks/git-guard.sh`.
8. **`.audit` location.** Its scratch location depends on regular touching.

## Bugs in earlier versions

| Bug | Versions | Evidence |
|---|---|---|
| Wrong exon order on the minus strand | v0.1.0, v0.2.0, `main`, `dev` to `c8b534f` | All such genes in three plant genomes |
| Chimeric multi-isoform CDS | Same | All 10,695 such genes in *A. thaliana* |
| Wrong-genome annotations at one thread | Same | `docs/validation/real-genomes/annotation-leak.md` |
| Empty annotations for RefSeq-style IDs | Same | *A. thaliana* |
| Failed runs reporting success, ignored options, misread values | Same | Tests |
| Statistical fixes | Same | `NEWS.md` |

Consequences and rerun steps: `docs/migrating-to-0.3.md`.

## Regression caught during this work

A pattern fix made during this release lost every GO term after the first in `%2C`-separated lists, the standard GFF3 encoding. Unit tests did not catch it, but comparison with real annotations did, and it was fixed with a test. It shows that real-genome comparison should remain part of release checks.

## Not done or not tested

- Docker, macOS, Windows, GitHub Actions: not available here.
- R from CRAN and Bioconductor failed on the tested HPC module, because of `rtracklayer` and OpenSSL. Conda is the recommended route.
- **Codex:** the sandbox blocked Singularity (no `/dev/fuse`), so the agent unpacked the image and ran it under Bubblewrap. It reports passing all 14 checks, but our disk checks could not confirm it, so this is not a clean pass. Its run found two real defects, both fixed.
- Real-genome tests cover three plant genomes. Multi-isoform dN/dS was not compared end to end, because only *A. thaliana* has isoforms here, and there was no suitable partner.
- `dnds_ideogram` content on real genomes was only confirmed to render.
- `docs/validation/known-issues.md` lists open items.

## Path to GO

1. Push `dev` (see the report in `PROGRESS.md`), and fix what CI finds.
2. Decide items 1, 2, 5, and 6.
3. Run the steps in `RELEASE.md`.
