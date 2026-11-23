# Enrichment: correction or restricted domain, and an independent validation plan (proposal)

Status: **proposal for the maintainer's decision (2026-11-23); nothing here is implemented or run, except the worst-cell follow-up already running (`enrichment-error-study.md`).** Requested by the maintainer: the 43-70% complete-null rates of `go_enrichment`, `term_enrichment` and `ipr_enrichment` are genuine problems; revised documentation does not resolve them; dndsR needs either an effective validated correction or a clearly restricted operating domain, with an independent validation plan for each affected method.

## Why the running follow-up targets `covariate_enrichment`

The failures of the three Fisher-type commands are already established with narrow intervals (complete-null false-term rates of 33-42% with exact labels and 43-70% with estimator labels at *A. thaliana*-*A. lyrata*; mixed-null false-discovery proportion 20-34%; studies [N](../validation/enrichment-error-calibration.md) and [confounding](../validation/enrichment-confounding.md)). More null runs of those commands cannot change their classification. The one open *calibration* question was about the proposed correction itself: `covariate_enrichment`'s worst cell sat exactly at the limit (5.0%, 2.2-11.2, 100 runs). The follow-up tests whether the correction holds there. It does not, and was not meant to, resolve the three commands' failures; this plan addresses them.

## Where each method stands

| Method | Hypothesis it tests correctly | For selection related to terms | Gaps |
|---|---|---|---|
| `term_enrichment` (Fisher, any term family) | call rate differs between term and universe (study N: FWER ≤1.3%, FDP ≤2.5% with exchangeable labels) | **fails** wherever the call depends on length or dS: 33-70% complete-null false-term rate at AtAl, 6% at AtCr, 3% at AtBr | no corrected counterpart exists for custom term families |
| `ipr_enrichment` (Fisher, InterPro) | as above | **fails** as above (InterPro 42% exact labels at AtAl) | `per_type` adds family-wise inflation (9-17% of runs) |
| `go_enrichment` (topGO, default `weight01`) | as above, with topGO's dependent tests (no proven BH error rate) | **fails**: 21% (exact labels), 43% (estimator labels) at AtAl | topGO's elimination algorithms cannot take a covariate model |
| `covariate_enrichment` (Poisson-binomial with length and dS, GO and InterPro) | selection related to terms given length and dS | **passes at `dNdS > 1`** on 4 plant pairs in 2 families (worst cell 5.0%, follow-up running), FDP ≤1.7% | **fails at `q0.95`** (11.0%, 6.3-18.6); polyploids, other clades, annotation completeness, gene dependence untested |

## Options for the three Fisher-type commands

1. **Correction inside the commands.** Add an optional `adjust = "length_ds"` to `term_enrichment` and `ipr_enrichment` that uses `covariate_enrichment`'s validated model (`pb_floor`) for any term family, so custom terms also get a corrected test; for `go_enrichment`, the adjustment can only apply to the classic (all-ancestor) test, not to `weight01`. Interface change; needs approval, tests and the validation below.
2. **Restricted operating domain.** Keep the commands unadjusted, but state their domain as *association of the dN/dS call with a term*, not selection, and add a diagnostic computed on every run: a likelihood-ratio test of whether the call depends on log length and log dS (the same spline model). Runs where it does are outside the domain for selection claims and say so in the output. The diagnostic itself must be validated: among cells it classes as in-domain, Fisher's complete-null rate must be at most 5%.
3. **Both:** option 1 as the route for selection questions, and option 2's diagnostic and wording for the unadjusted tests.

**Recommendation (for decision):** option 3. It gives a validated corrected test for every term family (including custom terms, which today have none), keeps the unadjusted tests honest about their hypothesis, and changes no default until the validation below passes. `covariate_enrichment` stays restricted to `dNdS > 1` unless a quantile threshold passes on independent data; `q0.95` is documented as not calibrated.

## Independent validation plan (criteria fixed before any run; data never used in development)

**Data.** At least three new public genome pairs outside Brassicaceae and Solanaceae, spanning divergence (median dS about 0.05-0.5) and including one recent polyploid: candidates *Glycine max*-*Phaseolus vulgaris* (paleopolyploid), *Oryza sativa*-*Brachypodium distachyon* (monocots), *Drosophila melanogaster*-*D. simulans* (animal, curated GO), each with GO and InterPro from the same pipeline; accessions confirmed before approval. Each pair's genes, lengths, dS and terms come from the documented dndsR pipeline.

**Scenarios per pair.** Complete nulls: estimator-noise labels (`sim-genome.R`: `long`, `bylen`), exact labels with graded dependence on length and dS (study N's N0 and N1 at half, full and double strength), reassignments (`permute`, `rotate`, `strata5`). Mixed nulls: spikes of 2%, 5% and 20% of genes in 5 terms. Thresholds `dNdS > 1` and `q0.95`.

**Methods.** `go_enrichment` (default and classic), `term_enrichment`, `ipr_enrichment` (global and per type), `covariate_enrichment`, and, if option 1 is approved, the adjusted modes; the option 2 diagnostic on every run.

**Criteria.** Per method, pair, scenario and threshold: complete-null share of runs with any false term at most 5% with the upper Wilson bound at most 7.5% (300 runs per cell: interval about ±2.5 points); mixed-null mean false-discovery proportion at most 5%; power reported, not judged. A method is *validated for selection inference* on a pair and threshold only if every null cell passes; the diagnostic is *validated* if Fisher's rate in its in-domain cells is at most 5% and it classes AtAl-like dependence as out of domain.

**Cost.** Dominated by the estimator on simulated sequences: about 3 pairs × 10 null cells × 300 runs plus spikes; to be measured on a pilot of 20 runs per cell before the full runs; roughly 300-600 core-hours. Requires approval.

## What happens to the release meanwhile

Until this is decided and validated, the three Fisher-type commands are classed as **calibration insufficient for selection inference**, usable only as descriptive call-term association with the existing notice; `covariate_enrichment` is **stable within limits** (`dNdS > 1`, the tested clades, pending the worst-cell follow-up).
