# Roadmap to v1.0.0

Goal: a scientifically defensible, reproducible, and accessible dndsR, suitable for a preprint. This lists what evidence shows is missing, not features. Items are dropped if they prove unnecessary.

## Critical scientific validation

| # | Item | Why | Status |
|---|---|---|---|
| 1 | Characterize the estimator | Interpretation depends on it | **Done:** `docs/validation/estimators.md`. No error found; compresses ratios, property of the method family |
| 2 | Quantitative criteria | Evidence should be judged by stated rules | **Done:** `docs/validation/criteria.md`; revisit as evidence grows |
| 3 | Public real-genome workflow | Publishable evidence without private data | **Done for CDS through dN/dS:** `docs/validation/public-genome.md` |
| 4 | Isoform sensitivity of dN/dS | Longest isoform is a convention | **Done for one pair:** range above 0.1 in 22% of multi-isoform pairs; 0.3% cross 1 (`docs/validation/public-genome.md`) |
| 5 | Uncertainty for per-gene dN/dS | Half of neutral genes exceed 1, so a bare ratio cannot support claims about individual genes | **Open, and the most important scientific decision.** Options below |
| 6 | Enrichment on a public annotated genome | Public evidence so far uses constructed cases | Open. Needs public GO or InterPro annotations for *A. thaliana* |
| 7 | Null calibration of enrichment | Selecting on `min_pos`, plus gene dependence, may inflate significance | Open. Permuting positive labels on the public pair gives empirical false-positive rates |
| 8 | Orthology accuracy | Best-hit pairing is a heuristic | Open. Compare RBH pairs with a public ortholog resource for *A. thaliana* and *A. lyrata* |
| 9 | Independently verify biological results | Consistency between tools is not biological truth | Open. Compare the public pair's dS and gene-level rank against published values, citing only what is confirmed |
| 10 | Gene-level `dnds_state_contrast` | Disabled | Remove it from the v1.0.0 interface unless a design passes null simulation. Decision for the maintainer: `state-contrast-design.md` |
| 11 | Saturation warning | `Comeron95` truncates undefined terms silently | Open. Small: flag high dS and truncated values in QC |

### Options for item 5

The smallest defensible remedy comes first:

1. **Strengthen language only.** Done in this release: above 1 is a candidate list.
2. **Report a per-gene bootstrap interval,** resampling codons, and base positive calls on the interval's lower bound as an optional stricter filter. It changes no default, and is simple and method-independent. **Recommended.**
3. **Add a pairwise likelihood-ratio test** with PAML codeml. More formal but adds a heavy dependency, and pairwise tests have low power and ω across sites is unmodeled.

Keeping the default as is but changing what enrichment means is a scientific decision for the maintainer.

## Critical reliability

| # | Item | Status |
|---|---|---|
| 12 | CI green on all platforms | **Done** at `f84521f`, including Apptainer. Keep green |
| 13 | Docker-built image passes container checks | **Done** at `f84521f`, including package tests in the image |
| 14 | Determinism | **Done for dN/dS** across thread counts. Add downstream steps with #7 |
| 15 | Documented examples run | **Done:** tutorial, R tutorial, and container checks. `your-own-data.md` uses placeholder files |
| 16 | Failures are explicit | **Done** for known paths, with tests |
| 17 | Skill tested at the release commit | Repeat at each release |

## Publication readiness

| # | Item | Status |
|---|---|---|
| 18 | Methods description | **Draft:** `docs/methods.md` |
| 19 | Reproducibility materials | Scripts exist. Package the public pair's commands as one runnable recipe |
| 20 | Benchmark of time and memory | Open. Slurm accounting for the public pair is a start |
| 21 | Zenodo DOI, CITATION | Prepared. Awaits approval |
| 22 | Keep unreleased findings out of public artifacts | In force, with a check |

## Optional, only if users need it

- Docker on macOS and Windows testing.
- A plugin or marketplace package for the skill.
- CRAN submission: `orthologr` is GitHub-only and blocks it for now.
- A `dndsr update` command.
- Other estimators in the image, such as KaKs_Calculator's. Useful for sensitivity analyses, not required.

## Not planned

New analysis types or new plots, before the above is complete.
