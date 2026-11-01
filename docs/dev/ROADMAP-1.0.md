# Roadmap to v1.0.0

v0.3.0 establishes a corrected, working foundation. v1.0.0 establishes a scientifically validated, reproducible, publication-ready tool. There are no intermediate versions for their own sake.

Scientific validity comes before new functionality. Existing behavior that is correct, tested, and documented is left alone. Items are dropped if the evidence shows they are unnecessary. The publication must state what dndsR can reliably estimate, and what biological conclusions cannot be established from those estimates alone.

## Priority 1: functional enrichment validation

Is the enrichment workflow statistically defensible on real annotated data, and what are the smallest corrections needed?

| Item | Approach | Status |
|---|---|---|
| Real annotations | Public GO or InterPro annotations for *A. thaliana*, through the whole workflow | **Done:** [enrichment.md](../validation/enrichment.md) |
| Null calibration | Empirical false-positive rates. Not only unrestricted label permutation, which assumes genes are exchangeable | **Done for one pair. Filter-then-BH fails; the correction needs a decision** |
| Candidate selection | How the threshold on a point estimate, `min_pos`, and `max_dnds` affect results | **Measured,** except `max_dnds` |
| Confounders | Annotation coverage, gene length, expression or conservation proxies, and dS. Selection on dN/dS may favor certain genes regardless of function | **Measured** for coverage, length, and dS. Expression and conservation are open |
| Dependence | Tandem clusters and shared terms, which Fisher's test assumes away | **Measured:** collapse and genomic rotation |

## Priority 2: orthology accuracy

| Item | Approach | Status |
|---|---|---|
| Accuracy of reciprocal best hits | Benchmark against independent, established orthology resources for the public pair | Open |
| Failure modes | Paralogs, gene duplication, isoform effects | Partly: isoform effect on dN/dS is measured, not its effect on pairing |
| Downstream effect | Whether errors materially change dN/dS conclusions | Open |

## Priority 3: independent estimator validation

Existing comparisons with PAML and KaKs_Calculator stay, and are expanded only where necessary.

| Item | Approach | Status |
|---|---|---|
| Bias and variance | Across biologically meaningful simulation regimes | Partly: `docs/validation/estimators.md` |
| Saturation and failure conditions | Where estimates stop being reliable, using dS and gene length | Open |
| Confidence intervals | First investigate whether intervals are useful and calibrated. **No bootstrap feature unless that is established** | Open, investigation only |
| Estimation against testing | State clearly that dN/dS is a descriptive estimate, not a test of selection | Started: `docs/interpretation.md` |

A neutral gene's point estimate falling above 1 about half the time is expected behavior, not an estimator defect. The issue to manage is users reading dN/dS above 1 as statistically significant evidence of positive selection.

## Priority 4: end-to-end public demonstration

| Item | Status |
|---|---|
| One fully reproducible public-genome workflow, including annotation and enrichment | Done through dN/dS: `docs/validation/public-genome.md`. Annotation and enrichment are open |
| Biological plausibility against independent published evidence, citing only what is confirmed | Open |
| The whole analysis as a simple, executable tutorial | Open |

## Priority 5: interface cleanup and QC

| Item | Status |
|---|---|
| Remove gene-level `dnds_state_contrast` from the supported interface, unless a statistically valid design is demonstrated | Open. See `state-contrast-design.md` |
| dS saturation, unreliable-estimate, and other QC warnings where evidence supports them | Open. `Comeron95` truncates undefined terms silently |

No unnecessary statistical machinery or dependencies.

## Keep true throughout

| Item | Status |
|---|---|
| CI passes on all platforms, including Apptainer and package tests in the image | Done at `f84521f`. Keep it so |
| Determinism | Done for dN/dS across thread counts. Add downstream steps with Priority 1 |
| Documented examples run, and failures are explicit | Done for known paths |
| Skill tested at the release commit | Repeat at each release |
| Human-native and agent-native: simple install, clear workflows, reproducible examples | Ongoing: `AGENT_NATIVE.md` |
| Unreleased findings stay out of public artifacts | In force, with a check |

## Publication materials

| Item | Status |
|---|---|
| Methods description | Draft: `docs/methods.md` |
| Results evidence | `docs/validation/`, growing with the priorities above |
| Benchmark of time and memory | Open. Slurm accounting for the public pair is a start |
| Zenodo DOI, CITATION | Prepared. Needs approval |

## Deferred, only if users need them

- Docker on macOS and Windows testing.
- A plugin or marketplace package for the skill.
- CRAN submission: `orthologr` is GitHub-only and blocks it for now.
- A `dndsr update` command.
- Other estimators in the image.
- New analysis types or new plots.
