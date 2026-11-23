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

## Resume here: study B completion, alignment pilot, multiple-testing study (2026-11-22)

The maintainer's seven rulings on the study B core results are in `DECISIONS.md` (2026-11-22). Reported, with proposals awaiting the maintainer: the alignment pilot ([report](../validation/alignment-pilot.md): no aligner nominated, cost criterion; option of documenting codon-aware user alignments) and the multiple-testing study ([report](../validation/lineage-multiple-testing.md): per-lineage `p_bh_lineage` and gene `p_adj` control false discoveries on correct alignments, every aBSREL output fails after realigning indels; proposed optional pooled-BH column, not implemented). Study B is complete (remaining arms reported 2026-11-23 in [lineage-robustness-core.md](../validation/lineage-robustness-core.md#remaining-arms-16-taxa-and-busted-2026-11-23): 16-taxon units confirm the 8-taxon conclusions; BUSTED is broken by a foreground corrupted segment, masking restores size at 12 points of power, and it is nearly unaffected by realigned indels). Jobs (all finished):

| Work | Design | Root | Job |
|---|---|---|---|
| Study B, remaining arms: 16-taxon units, every BUSTED arm, BE3h for the unrun indel arms (about 950 core-hour ceiling) | `lineage-robustness-study.md`, amendment of 2026-11-22 | `study-B3` (lib `49fd39f`), units `phylo-study-B-units-rest.tsv` | array 21578592 |
| Alignment pilot: MAFFT (A0) against PRANK codon (A1) and MACSE (A2), true alignment for diagnosis; BE3h genes 1-150 (151-300 held out) and BE3hS (seed 321) | `alignment-pilot.md` (40 core-hour ceiling) | `aln-pilot` (lib `6962b4a`) | array 21580170 |
| Multiple-testing mixture study: BMIX and BMIXh, 1,000 genes each, plus aBSREL on BMIXh's true alignments | `lineage-multiple-testing.md` (50 core-hour ceiling) | `study-Q` (lib `48aa75d`) | arrays 21579083, 21581106 |

Summaries: `phylo-study-B-summary.R` and `phylo-study-B-exclude.R` (B), `aln-pilot-summary.R`, `lineage-mt-summary.R`. Done this session: scoped `Alignment check` warning and provenance entry (`820abf3`); docs on input errors, masking and the aBSREL p-value families; two corrections to the study B core report: `p_bh_lineage` and aBSREL's gene `p_adj` do not control false discoveries under realignment error (`48aa75d`), and the mixture's per-branch truth (`81a3007`). Not approved: full study L, the RELAX confirmation.

## Resume here: lineage-level analyses (`feat/lineage-selection`)

Updated 2026-11-21. The four lineage commands go into v0.5.0, keeping that number (maintainer, 2026-11-21, `DECISIONS.md`). They are, on `feat/lineage-selection` (not pushed; the maintainer pushes and opens a pull request into `dev`): `phylo_branch_omega` (codeml), `phylo_absrel`, `phylo_relax`, `phylo_summarize`, CLI and R, PAML 4.10.10 in the Dockerfile and conda environment.

- Engine equivalence: [validation/phylo-lineage.md](../validation/phylo-lineage.md). Candidate image built from `0c7c813` with PAML and HyPhy inside: container checks pass, 248 tests, 0 failed, 0 skipped.
- Calibration in simulation, complete (27 of 27 units, 234.6 core-hours): [validation/phylo-lineage-calibration.md](../validation/phylo-lineage-calibration.md). aBSREL's `p_raw` and `p_bh_lineage` calibrated, its gene `p_value` on T5 only. Follow-up ([design](phylo-lineage-followup.md), 110.6 core-hours, same report): RELAX on 16 taxa anticonservative and reproduced (8.3-9.3%), not explained by dndsR, the simulator, or starts; `p_fg_omega_gt1` at omega 1 5.5% and 6.8%, not established; `p_free_vs_m0` with missing taxa not reproduced (4%), not established. These carry a `calibration` column, warning, and provenance entry. R1 (direct HyPhy on dndsR's files within 1e-3) failed on 14 of 300 genes: HyPhy's RELAX optimum is not reproducible to 1e-3.
- Validation studies approved 2026-11-21 (designs fixed before running in `docs/dev/`): N (enrichment error criterion) done, [report](../validation/enrichment-error-calibration.md); R stage 1 done, [report](../validation/relax-calibration-stage1.md): no-go for the 930 core-hour design, a ~60 core-hour confirmation recommended and not approved; L pilot failed gate 3 (LS0 t7 −10.5% against 10%), full L held for the maintainer; B core arms done, [report](../validation/lineage-robustness-core.md): a corrupted segment breaks aBSREL and the branch models, is detected by `dissimilar_windows`, and masking restores calibration; realigned indels inflate aBSREL's `p_raw` and gene `p_value` (BE3h, homologous insertions: 7.1% and 21%, calibrated on the true alignment, so alignment error), `p_bh_lineage` calibrated; paralogs, discordance, missing taxa, recombination robust; every failure reproduced by HyPhy or codeml run directly. Four recommendations (warn, masking, document, none) await the maintainer; 16-taxon units and BUSTED arms held. A simulator defect in `add_indels()` (non-homologous internal insertions) was found, recorded, and given an option without changing the default (`49fd39f`). P and E are drafts. No interface or default changed.
- Still open: push and the pull request into `dev` (maintainer); CI on the pushed branch; a candidate image rebuilt from the final commit (the one checked, from `0c7c813`, predates the calibration flags); whether RELAX should get a recalibrated null (for example a parametric bootstrap), a method change for the maintainer. `rst`, `rst1`, `rub` deleted with approval (empty, unreferenced).

## Resume here: v0.5.0 (on hold)

Updated 2026-11-20. The maintainer approved the v0.5.0 scope and its preparation; the checklist is in [RELEASE-READINESS.md](RELEASE-READINESS.md). Nothing is merged, tagged, or published.

1. After the maintainer pushes `dev`: read the container workflow's annotations on the pushed commit (public API, `check-runs/<job>/annotations`). It must show `Package tests` with `skipped=0` and no `Skipped test` warning, because `DNDSR_FAIL_ON_SKIP=1` now makes any skip fail in the freshly built image, and the new container checks run KaKs_Calculator (`--dnds-method YN`) and `phylo_busted` with the recommended settings. All workflows must pass.
2. The MK maintainer corrects the MK line in `skills/dndsr/references/gaps.md` (channel message 11); then CI again on that commit.
3. Then, with approval: merge to `main`, tag `v0.5.0`, image, digests, release, Zenodo (`RELEASE.md`).
4. Separate from v0.5.0: lineage-level work in the maintainer's order, PAML branch models first ([phylo-lineage-plan.md](phylo-lineage-plan.md)); aBSREL's baseline per-branch ω is a candidate cross-check, not a substitute, until an equivalence study.
5. Heavier commands go through Slurm, never the login node. Agents do not push (the guard refuses even when asked); give the maintainer the command.

## Next

1. Done: the maintainer pushed `dev`, which swaps the v0.4.0 digest into `README.md`, `INSTALL.md`, the tutorial, and the skill (RELEASE.md step 3), and merges it to `main` by pull request when convenient. Until then, `main`'s documents name the v0.3.0 image, which has the same dependencies.
2. The v0.4.0 GitHub release has no assets. RELEASE.md step 4 expects the workflow's `release-manifest.json` and skill tarball attached; the maintainer can attach them from the `publish` run, 41565720097.
3. This software is frozen for the Chenopodium analysis, which needs no rerun: `verify_release_equivalence.sh` reports the analysis's software identical to `v0.4.0`. Its work and results live outside this repository, and nothing from them comes here.
4. For v1.0.0: `docs/dev/ROADMAP-1.0.md`. Do not start it for the current manuscript.
5. Phylogenetic module (`phylo_prepare`, `phylo_busted`, HyPhy BUSTED): P0 to P2 are on `dev`, experimental and unreleased, at development version 0.4.0.9000. P3 pilot and real-genome benchmark done: [docs/validation/phylo-calibration-pilot.md](../validation/phylo-calibration-pilot.md). The maintainer approved a staged calibration study: plan in [phylo-calibration-stage2.md](phylo-calibration-stage2.md). Stage 2A (21,270 fits, 0 failed) is done and reported in [phylo-calibration-stage2a.md](../validation/phylo-calibration-stage2a.md), with a preliminary recommendation. The maintainer approved stage 2B's essential blocks: [phylo-calibration-stage2b.md](phylo-calibration-stage2b.md). Its smoke test passed (386 core-hours, recorded as spent; 1,806 is the prospective ceiling for the rest of 2B). A pre-launch budget and execution audit moved filtered refits to an 8-core array, added duplicate-submission protection, and fixed two Windows-only test failures; see the design's change log. Stage 2B is complete: 115 of 115 units, 0 failed, 1,475.5 core-hours (B2 with 3 starts, since D4 failed in one cell). Report, with the stage 2 recommendations for the maintainer's review: [phylo-calibration-stage2b.md](../validation/phylo-calibration-stage2b.md). No default changed, and masking is not in the package. No default changed. Design and decisions: [phylo-design.md](phylo-design.md), [phylo-calibration-design.md](phylo-calibration-design.md). No published image has HyPhy; development images were built and tested locally only. The v1.0.0 roadmap is unchanged.
6. `calculate_dnds --dnds-method` offers NG, LWL, LPB, MLWL, YN, and MYN besides the default Comeron, through KaKs_Calculator 2.0.1, development code only. Validation: [docs/validation/estimators.md](../validation/estimators.md), "Other methods". No published image has KaKs_Calculator.
7. McDonald-Kreitman test (development code, owned by the MK maintainer, a collaborator; PR #5-#8 merged). Since PR #8 (`b4522b2`), sites are classed by allele for VCF and BED (BED needs the ALT allele in column 4), and `snps_midfold`/`div_steps_midfold` are gone; this is the MK owner's choice in their area. Half-specified MK columns fail fast, with a regression test. Discussion: `docs/dev/channels/mcdonald-kreitman.md`.
8. Local gate: `tools/dev/check.sh --rcmd-check` adds `R CMD check`, and is for pull requests that touch `R/`. The safeguard tests now cover `feat/*` pushes. CI: every workflow passed at `b65bf8b`, Windows included. Next-release readiness, with the real-genome rerun (v0.4.0 against `dev`, byte-identical pairwise outputs) and the blockers: [RELEASE-READINESS.md](RELEASE-READINESS.md). The maintainer approved the stage 2 phylogenetic guidance (`DECISIONS.md`), now in `docs/phylo-methods.md`; no default changed.
