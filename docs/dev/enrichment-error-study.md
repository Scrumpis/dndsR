# Study N: which error the enrichment tests make, and why (design)

Status: **design fixed before any run, 2026-11-21.** Requested by the maintainer after the [validation audit](validation-audit.md) reported Fisher-type enrichment as failing calibration, with "a false term in up to 43-70% of runs".

## What the earlier figure measured

From the code that produced it (`enrichment-methods-summary.py`; [enrichment-confounding.md](../validation/enrichment-confounding.md)), the 43-70% is the **share of runs with at least one significant term (BH 0.05)** in scenarios where omega has no link to any term. In such a run every discovery is false, so this share is both the family-wise error rate and the FDR (FDR = E[V/R; R > 0] = P(V > 0) when every hypothesis is null). Spiked scenarios report a separate mean false discovery proportion, and power on spiked InterPro terms.

So the earlier figure is an FDR, but of **which null**. Labels came from estimator noise, which depends on gene length and dS. Label and term are then genuinely dependent for any term whose genes are short or slowly evolving, so the hypothesis Fisher's test examines (the label is independent of the term) is **false** for those terms, while the scientific one (selection is unrelated to the term) is true. Existing evidence already separates the two in part: under exchangeable reassignments (`permute`, `rotate`), Fisher's test with BH gave a term in 0.7-2.0% of runs at AtAl and `go_enrichment`'s defaults in 1.0% (`enrichment_methods_*.tsv`, `enrichment_go_null.tsv`), while under length-linked labels (`bylen`) it gave 43-70%.

What is missing, and what this study measures with exact truth (no estimator, no shuffle): the complete-null FWER, the mixed-null FDR, and power, separately, under (a) exchangeable labels and (b) labels that depend on length and dS only; and whether each discovery is a true association of the label with the term (the test's own hypothesis) or not.

## Hypotheses

- **H1.** With exchangeable labels, Fisher's test with BH (`term_enrichment`, `ipr_enrichment`, `go_enrichment(algorithm = "classic")`), `go_enrichment`'s default `weight01`, and `covariate_enrichment` control the complete-null FWER and the mixed-null FDR at 5%, despite nested GO terms and overlapping InterPro entries (dependence) and discrete p-values.
- **H2.** With labels that depend only on length and dS, Fisher's test reports terms whose label rate truly differs (relative risk above 1): it is right about its own hypothesis and wrong about the scientific one. `covariate_enrichment` controls both error rates for the scientific hypothesis.
- Alternatives: H1 fails (a problem in correction or dependence, a software or statistical defect), or H2 fails for `covariate_enrichment` (its adjustment is inadequate when the true dependence is not its spline form).

## Design

Real genes of *A. thaliana* against *A. lyrata* (the table of [enrichment.md](../validation/enrichment.md), `56321306d1b52a23`; 22,091 genes, GO with ancestors, InterPro), with lengths, dS, and terms unchanged. Labels are drawn directly (`enrichment-direct-labels.R`):

| Scenario | P(positive) | Terms with an effect | Runs |
|---|---|---|---|
| N0 | constant, the observed rate (206 / 22,091) | none | 300 |
| N1 | logistic in log length and log dS (quadratics and their product, fitted to the observed labels; not the spline basis `covariate_enrichment` uses), scaled to the same expected count | none | 300 |
| N0s_0.02, N0s_0.05 | N0, then 2% or 5% of the genes of 5 InterPro terms (20-300 genes) made positive | 5 | 200 each |
| N1s_0.02, N1s_0.05 | N1, the same spikes | 5 | 200 each |

Tests: the candidates of `enrichment-methods.R` (`fisher`, which equals `term_enrichment` and `go_enrichment(algorithm = "classic")`; `pb_floor`, which equals `covariate_enrichment`; and the others, reported for completeness), and the real `go_enrichment()` defaults on the first 100 runs of N0 and N1 and the first 50 of each mixed scenario. BH within a family (GO: all ontologies; InterPro), as dndsR's defaults do. Seeds 960,000 + 1,000 x scenario index + run, new to this study.

## Criteria (fixed before running)

For a method and family, in a scenario:

- **Complete null (N0, N1):** FWER (= FDR) **controlled** if the point estimate is at most 5% and the 95% Wilson upper bound at most 7.5% (300 runs: an interval of about plus or minus 2.5 points at 5%); **exceeded** if the lower bound is above 5%; otherwise inconclusive.
- **Mixed null:** mean FDP controlled by the same rule, with a normal interval over runs.
- **Power:** share of spiked terms found, reported with no threshold; compared between methods on the same runs.
- **H2, the test's own hypothesis:** among Fisher's false discoveries in N1, the median true relative risk, and the share above 1.1. H2 holds for Fisher if the median exceeds 1.1.

Conclusions follow from these: H1 holding for Fisher and H2 holding means its excess is confounding by length and dS (it answers a different question), not a failure of multiple-testing correction; H1 failing would be a defect to fix.

## Also audited (code, no runs)

Universe (genes with at least one term of the family, per GO ontology as in topGO); how `min_total` and `min_pos` enter the family (before and after correction); per-type correction (`ipr_enrichment` per entry type: each type is its own family, so the chance of a term somewhere grows with the number of types, while the FDR in each family is controlled); `qvalue` and IHW options on discrete p-values; the exclusion of genes with dS = 0 or dN/dS at least 10 (gene-selection effects that act through length and dS).

## Cost

About 1,400 runs of `enrichment-methods.R`, measured in the smoke test below; projected under 50 core-hours.

## Follow-up: the worst cell (2026-11-23, fixed before running)

The audit's gap 7: `covariate_enrichment`'s worst null cell in the confirmation, InterPro under the `strata5` reassignment at *A. thaliana*-*A. lyrata*, was 5.0% (2.2-11.2) on 100 runs, meeting the rule with no margin. The maintainer asked that the remaining enrichment calibration be finished. **Design:** 1,000 new runs of that cell (`strata5:S`, seeds 970,000-970,999, never used), through the frozen study N code (`760ca03`) with `enrichment-methods.R`, so every candidate is computed on the same labels; `pb_floor` is `covariate_enrichment`. Runner `enrichment-worstcell.sb`. **Measure:** share of runs with at least one significant InterPro term (BH 0.05; in a reassignment null every term is false, so this is both the family-wise error rate and the FDR), Wilson 95% interval; GO family and Fisher's test reported beside it. **Criterion:** holds if the point estimate is at most 5% and the upper bound at most 7.5%; exceeds if the lower bound is above 5%; otherwise inconclusive. **Cost:** about 78 core-seconds per run, about 22 core-hours.
