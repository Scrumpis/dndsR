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

## Next

1. The maintainer pushes `dev`, which swaps the v0.4.0 digest into `README.md`, `INSTALL.md`, the tutorial, and the skill (RELEASE.md step 3), and merges it to `main` by pull request when convenient. Until then, `main`'s documents name the v0.3.0 image, which has the same dependencies.
2. The v0.4.0 GitHub release has no assets. RELEASE.md step 4 expects the workflow's `release-manifest.json` and skill tarball attached; the maintainer can attach them from the `publish` run, 41565720097.
3. This software is frozen for the Chenopodium analysis, which needs no rerun: `verify_release_equivalence.sh` reports the analysis's software identical to `v0.4.0`. Its work and results live outside this repository, and nothing from them comes here.
4. For v1.0.0: `docs/dev/ROADMAP-1.0.md`. Do not start it for the current manuscript.
5. Phylogenetic module (`phylo_prepare`, `phylo_busted`, HyPhy BUSTED): P0 to P2 are on `dev`, experimental and unreleased, at development version 0.4.0.9000. P3 pilot and real-genome benchmark done: [docs/validation/phylo-calibration-pilot.md](../validation/phylo-calibration-pilot.md). The maintainer approved a staged calibration study: plan in [phylo-calibration-stage2.md](phylo-calibration-stage2.md). Stage 2A (21,270 fits, 0 failed) is done and reported in [phylo-calibration-stage2a.md](../validation/phylo-calibration-stage2a.md), with a preliminary recommendation. The maintainer approved stage 2B's essential blocks: [phylo-calibration-stage2b.md](phylo-calibration-stage2b.md). Its smoke test passed (386 core-hours, recorded as spent; 1,806 is the prospective ceiling for the rest of 2B). A pre-launch budget and execution audit moved filtered refits to an 8-core array, added duplicate-submission protection, and fixed two Windows-only test failures; see the design's change log. Main blocks finished (101 of 101 units, 1,194.5 core-hours); D4 failed in one of 12 cells (B6 M4, C1), so B2 runs with 3 starts as prespecified, approved and launched. Job IDs: the design's change log and `<root>/jobs.txt` (root `dndsr-runtime/p3s2b`). The 2A/2B report follows B2. No default changed. Design and decisions: [phylo-design.md](phylo-design.md), [phylo-calibration-design.md](phylo-calibration-design.md). No published image has HyPhy; development images were built and tested locally only. The v1.0.0 roadmap is unchanged.
6. `calculate_dnds --dnds-method` offers NG, LWL, LPB, MLWL, YN, and MYN besides the default Comeron, through KaKs_Calculator 2.0.1, development code only. Validation: [docs/validation/estimators.md](../validation/estimators.md), "Other methods". No published image has KaKs_Calculator.
7. McDonald-Kreitman test (development code, owned by the MK maintainer, a collaborator; PR #5, #6, #7 merged). `dev` counts 0-fold/4-fold for VCF and BED, and fails fast on half-specified MK columns (`87e5c02`). The method is the MK owner's to develop: an allele-aware rewrite (`9880f07`) is on their branch, not in `dev`, and comes through pull-request review. Discussion and open requests: `docs/dev/channels/mcdonald-kreitman.md` (message 4: a regression test for the fail-fast is requested from the MK side).
8. CI: every workflow passed at `c78abfa`; at `6b3245f` (PR #7) R-CMD-check, install, and skill pass. `check.sh` now also runs `tools/dev/check-imports.R` (undeclared `pkg::` use). Fixed after PR #7: `calculate_dnds` ignored a header row in comparison files since v0.2.0 (NEWS).
