# Roadmap to v1.0.0

v0.3.0 establishes a corrected, working foundation. v1.0.0 establishes a scientifically validated, reproducible, publication-ready tool. There are no intermediate versions for their own sake.

Scientific validity comes before new functionality. Existing behavior that is correct, tested, and documented is left alone. Items are dropped if the evidence shows they are unnecessary. The publication must state what dndsR can reliably estimate, and what biological conclusions cannot be established from those estimates alone.

**Maintainer's order after v0.3.0:** gene-length and dS confounding, independent estimator benchmarking, and empirical calibration of enrichment statistics come first.

## Priority 1: functional enrichment validation

Is the enrichment workflow statistically defensible on real annotated data, and what are the smallest corrections needed?

| Item | Approach | Status |
|---|---|---|
| Real annotations | Public GO or InterPro annotations for *A. thaliana*, through the whole workflow | **Done:** [enrichment.md](../validation/enrichment.md) |
| Null calibration | Empirical false-positive rates. Not only unrestricted label permutation, which assumes genes are exchangeable | **Done, at 3 divergences and a held-out pair, with labels from simulated evolution:** [enrichment-confounding.md](../validation/enrichment-confounding.md) |
| Candidate selection | How the threshold on a point estimate, `min_pos`, and `max_dnds` affect results | **Measured,** except `max_dnds` |
| Confounders | Annotation coverage, gene length, expression or conservation proxies, and dS. Selection on dN/dS may favor certain genes regardless of function | **Length and dS modeled:** opt-in `covariate_enrichment`. Defaults are unchanged, pending the maintainer's decision. The excess in defaults is larger than the 9.5% recorded earlier. Annotation coverage, expression, and conservation are open |
| Dependence | Tandem clusters and shared terms, which Fisher's test assumes away | **Measured:** collapse and genomic rotation |

## Priority 2: orthology accuracy

| Item | Approach | Status |
|---|---|---|
| Accuracy of reciprocal best hits | Benchmark against independent, established orthology resources for the public pair | **Done, 2 pairs against OMA:** [orthology.md](../validation/orthology.md) |
| Failure modes | Paralogs, gene duplication, isoform effects | Measured for dS and for a triplicated genome. Isoform effect on pairing is open |
| Downstream effect | Whether errors materially change dN/dS conclusions | **Measured** on term lists: most terms stay, a few do not |

## Priority 3: independent estimator validation

Existing comparisons with PAML and KaKs_Calculator stay, and are expanded only where necessary.

| Item | Approach | Status |
|---|---|---|
| Bias and variance | Across biologically meaningful simulation regimes | **Done:** 540 settings, 6 estimators, with site variation: `docs/validation/estimators.md` |
| Saturation and failure conditions | Where estimates stop being reliable, using dS and gene length | **Done** |
| Confidence intervals | First investigate whether intervals are useful and calibrated. **No bootstrap feature unless that is established** | **Investigated. Not established, so no feature** |
| Estimation against testing | State clearly that dN/dS is a descriptive estimate, not a test of selection | Started: `docs/interpretation.md` |

A neutral gene's point estimate falling above 1 about half the time is expected behavior, not an estimator defect. The issue to manage is users reading dN/dS above 1 as statistically significant evidence of positive selection.

## Priority 4: end-to-end public demonstration

| Item | Status |
|---|---|
| One fully reproducible public-genome workflow, including annotation and enrichment | **Done:** [public-genomes.md](../tutorials/public-genomes.md), run by `tools/dev/test-public-tutorial.sb` |
| Biological plausibility against independent published evidence, citing only what is confirmed | **Started, one family, chosen post hoc:** supporting, not validation. Needs a list chosen in advance |
| The whole analysis as a simple, executable tutorial | **Done** |

## Priority 5: interface cleanup and QC

| Item | Status |
|---|---|
| Remove gene-level `dnds_state_contrast` from the supported interface, unless a statistically valid design is demonstrated | Open. See `state-contrast-design.md` |
| dS saturation, unreliable-estimate, and other QC warnings where evidence supports them | **Done as notices** from `calculate_dnds`. `Comeron95` still truncates undefined terms silently |

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
| Methods description | Draft: `docs/methods.md`. Protocol for a publication: `docs/publication-protocol.md` |
| Results evidence | `docs/validation/`, growing with the priorities above |
| Benchmark of time and memory | Measured for the public pair, in [the tutorial](../tutorials/public-genomes.md) and [performance](../validation/performance.md) |
| Zenodo DOI, CITATION | Prepared. Needs approval |

## Deferred, only if users need them

- Docker on macOS and Windows testing.
- A plugin or marketplace package for the skill.
- CRAN submission: `orthologr` is GitHub-only and blocks it for now.
- A `dndsr update` command.
- Other estimators in the image.
- New analysis types or new plots.

## Scaling

Every step should use the CPUs it is given, as far as the work allows. Measurements and order: [performance.md](../validation/performance.md). Changes must leave every output byte-identical, checked by `test-threads.R` and on the public pair. Start after v0.3.0 is tagged.
