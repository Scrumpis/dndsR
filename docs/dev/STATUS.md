# Project status

Update this at the end of each substantive session. Identify things by commit and hash, not by date.

**Phase:** v0.4.0 is public: tag on `main` at `957cfe1` (merge of PR #4), image `scrumpis/dndsr@sha256:f1b42d9d338b406ec23daca1bb146ed0329db29d21a24a6871158996aa0572da`, version DOI 10.5281/zenodo.23903253 (concept DOI 10.5281/zenodo.23879469, used in `README.md` and `CITATION.cff`, resolves to the latest). Its software trees (`R`, `inst`, `NAMESPACE`, `DESCRIPTION`) are identical to `5983405`. `dev` names the v0.4.0 image. Checkpoints are in [PROGRESS.md](PROGRESS.md).

## Identity

dndsR is a real R package, fit for CRAN one day, and a containerized CLI, for people and AI agents equally, with ease of use first.

## Refs

| Ref | Commit | Notes |
|---|---|---|
| `v0.2.0` | `b8ea5ff` | Prerelease |
| `v0.3.0` | `d5f4b21` | Release, image `sha256:feb5e4b1` |
| `v0.4.0` | `957cfe1` | Release, image `sha256:f1b42d9d` |
| `main` | `957cfe1` | `v0.4.0` |
| `dev` | See `git log` | `main` was merged into it at `3f86d42`, by the maintainer |

## Resume here: v0.5.0

Updated 2026-11-20. The maintainer set the next release as **v0.5.0** (`DECISIONS.md`): pairwise pipeline, alternative estimators, reliability flags, and, marked experimental, `phylo_prepare`/`phylo_busted` and McDonald-Kreitman. The version is bumped on `dev`; nothing is merged, tagged, or published.

1. Done: `p_best_fits_adj` (BH over `p_best_fits`, same family as `p_adj`, which is unchanged); docs, tutorial, and skill recommend `--multiple-hits Double --extra-seeds 2,3` and `p_best_fits_adj`; release labels; release notes draft (`release-notes-0.5.0.md`); readiness (`RELEASE-READINESS.md`).
2. Blocking, for the maintainer: the second author's preferred name and ORCID; accepting the fixture's NCBI terms (`tools/dev/fixtures/README.md`). For the MK maintainer: the stale MK line in `skills/dndsr/references/gaps.md` (channel message 10).
3. Then, with approval: `RELEASE.md` step 5 (clone commands to `--branch v0.5.0`), `tools/dev/test-agent.sh claude`, merge, tag, image, digests.
4. Future, not in v0.5.0: lineage-level selection for a comparative project (`phylo_branch_omega`, `phylo_absrel`, `phylo_relax`, `phylo_summarize`). Review and recommended order: [phylo-lineage-plan.md](phylo-lineage-plan.md). First decision: the engine for branch ω (PAML, a new dependency, or HyPhy).
5. Heavier commands go through Slurm, never the login node (maintainer, 2026-11-20).

## Next

1. Done: the maintainer pushed `dev`, which swaps the v0.4.0 digest into `README.md`, `INSTALL.md`, the tutorial, and the skill (RELEASE.md step 3), and merges it to `main` by pull request when convenient. Until then, `main`'s documents name the v0.3.0 image, which has the same dependencies.
2. The v0.4.0 GitHub release has no assets. RELEASE.md step 4 expects the workflow's `release-manifest.json` and skill tarball attached; the maintainer can attach them from the `publish` run, 41565720097.
3. This software is frozen for the Chenopodium analysis, which needs no rerun: `verify_release_equivalence.sh` reports the analysis's software identical to `v0.4.0`. Its work and results live outside this repository, and nothing from them comes here.
4. For v1.0.0: `docs/dev/ROADMAP-1.0.md`. Do not start it for the current manuscript.
5. Phylogenetic module (`phylo_prepare`, `phylo_busted`, HyPhy BUSTED): P0 to P2 are on `dev`, experimental and unreleased, at development version 0.4.0.9000. P3 pilot and real-genome benchmark done: [docs/validation/phylo-calibration-pilot.md](../validation/phylo-calibration-pilot.md). The maintainer approved a staged calibration study: plan in [phylo-calibration-stage2.md](phylo-calibration-stage2.md). Stage 2A (21,270 fits, 0 failed) is done and reported in [phylo-calibration-stage2a.md](../validation/phylo-calibration-stage2a.md), with a preliminary recommendation. The maintainer approved stage 2B's essential blocks: [phylo-calibration-stage2b.md](phylo-calibration-stage2b.md). Its smoke test passed (386 core-hours, recorded as spent; 1,806 is the prospective ceiling for the rest of 2B). A pre-launch budget and execution audit moved filtered refits to an 8-core array, added duplicate-submission protection, and fixed two Windows-only test failures; see the design's change log. Stage 2B is complete: 115 of 115 units, 0 failed, 1,475.5 core-hours (B2 with 3 starts, since D4 failed in one cell). Report, with the stage 2 recommendations for the maintainer's review: [phylo-calibration-stage2b.md](../validation/phylo-calibration-stage2b.md). No default changed, and masking is not in the package. No default changed. Design and decisions: [phylo-design.md](phylo-design.md), [phylo-calibration-design.md](phylo-calibration-design.md). No published image has HyPhy; development images were built and tested locally only. The v1.0.0 roadmap is unchanged.
6. `calculate_dnds --dnds-method` offers NG, LWL, LPB, MLWL, YN, and MYN besides the default Comeron, through KaKs_Calculator 2.0.1, development code only. Validation: [docs/validation/estimators.md](../validation/estimators.md), "Other methods". No published image has KaKs_Calculator.
7. McDonald-Kreitman test (development code, owned by the MK maintainer, a collaborator; PR #5-#8 merged). Since PR #8 (`b4522b2`), sites are classed by allele for VCF and BED (BED needs the ALT allele in column 4), and `snps_midfold`/`div_steps_midfold` are gone; this is the MK owner's choice in their area. Half-specified MK columns fail fast, with a regression test. Discussion: `docs/dev/channels/mcdonald-kreitman.md`.
8. Local gate: `tools/dev/check.sh --rcmd-check` adds `R CMD check`, and is for pull requests that touch `R/`. The safeguard tests now cover `feat/*` pushes. CI: every workflow passed at `b65bf8b`, Windows included. Next-release readiness, with the real-genome rerun (v0.4.0 against `dev`, byte-identical pairwise outputs) and the blockers: [RELEASE-READINESS.md](RELEASE-READINESS.md). The maintainer approved the stage 2 phylogenetic guidance (`DECISIONS.md`), now in `docs/phylo-methods.md`; no default changed.
