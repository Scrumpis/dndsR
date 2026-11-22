# Enrichment: which error the tests make (study N)

**Question.** The validation audit reported Fisher-type enrichment (`term_enrichment`, `ipr_enrichment`, `go_enrichment`) as failing calibration, with a false term in 43-70% of simulated genomes ([enrichment-confounding.md](enrichment-confounding.md)). The maintainer asked which error criterion that was, and whether it is a failure of FDR control or something else. Design, fixed before running: [docs/dev/enrichment-error-study.md](../dev/enrichment-error-study.md).

**Answer.** It is not a failure of multiple-testing correction. With labels exchangeable across genes, every test controls the complete-null family-wise error and the mixed-null FDR at 5%, despite nested GO terms, overlapping InterPro entries, and discrete p-values. The excess appears only when the chance of a gene passing dN/dS > 1 depends on its length and dS. Then the hypothesis the Fisher-type tests examine (whether the label is independent of the term) is genuinely false for terms of short or slowly evolving genes, and the tests correctly report those terms: their false discoveries have a median true relative risk of 5 to 10. The scientific hypothesis (selection unrelated to the term, given length and dS) is true, so for that question they are confounded. `covariate_enrichment`, which adjusts for length and dS, controls both error rates on the same labels, with power at least equal.

## What the earlier figure measured

`enrichment-methods-summary.py` counts a run as false if any term is significant at BH 0.05 in a scenario where omega has no link to annotation. Every discovery is then false, so that share is the family-wise error rate and also the FDR (FDR = P(V > 0) when all hypotheses are null). Its labels came from the real estimator applied to simulated sequences, so they inherit the estimator's dependence on length and dS. That figure measures confounding, not BH.

## Design

The *A. thaliana*-*A. lyrata* table (`56321306d1b52a23`) keeps its genes, lengths, dS, GO terms (with ancestors) and InterPro entries. Labels are drawn directly with known probabilities (`enrichment-direct-labels.R`), so there is no estimator and the truth is exact:

- **N0**: P(positive) constant (206 / 22,091, as observed). Every term is null in every sense.
- **N1**: logistic in log length and log dS, with quadratic and product terms. It is fitted to the observed labels and scaled to the same expected count, using a basis unlike `covariate_enrichment`'s splines. No term has an effect.
- **N0s_F, N1s_F**: the same, then a share F (2% or 5%) of the genes of 5 InterPro terms with 20-300 genes are made positive.

Runs: 300 per complete null, 200 per mixed scenario, and the real `go_enrichment()` defaults on 100 and 50 of them. Seeds 960,000 and up, new. The tests are the candidates of `enrichment-methods.R`, unchanged from the earlier study: `fisher` equals `term_enrichment`, `ipr_enrichment` and `go_enrichment(algorithm = "classic")`, and `pb_floor` equals `covariate_enrichment`. There are 1,400 runs, 18.7 core-hours, image `dndsr_v0.4.0_sha256-f1b42d9d`, code at `760ca03`. Summary: `enrichment-direct-summary.R`; table: [data/enrichment_direct_cells.tsv](data/enrichment_direct_cells.tsv).

## Results

Complete nulls: share of runs with at least one significant term (FWER, equal to FDR here), 95% Wilson interval. Criterion fixed in advance: controlled if at most 5% with upper bound at most 7.5%; exceeded if the lower bound is above 5%.

| Test | N0 GO | N0 InterPro | N1 GO | N1 InterPro |
|---|---:|---:|---:|---:|
| `fisher` (`term_enrichment`, `ipr_enrichment`, classic GO) | 1.3% (0.5-3.4) | 1.3% (0.5-3.4) | **33.0% (27.9-38.5)** | **41.7% (36.2-47.3)** |
| `go_enrichment()` defaults (`weight01`) | 0 of 100 (0-3.7) | not applicable | **21% (14.2-30.0)** | not applicable |
| `covariate_enrichment` (`pb_floor`) | 0.3% (0.1-1.9) | 1.0% (0.3-2.9) | **1.0% (0.3-2.9)** | **0.3% (0.1-1.9)** |
| `cmh_10x10` (stratified exact) | 1.3% | 0.7% | 2.0% (0.9-4.3) | 2.0% (0.9-4.3) |

- **H1, holds.** Under exchangeable labels every test is controlled, with large margins. Correction, term dependence and discreteness cause no excess. Every other candidate in the table file is controlled under N0 too.
- **H2, holds.** Under N1, Fisher's false discoveries are terms whose label rate truly differs. The median true relative risk (mean P(positive) in the term over the rest of its universe) is 5.0 for GO and 9.8 for InterPro; 96.5% and 99.6% are above 1.1. The `go_enrichment` defaults' false discoveries likewise have a median relative risk of 3.3, all above 1.1. `covariate_enrichment` is controlled.

Mixed nulls: mean false discovery proportion (95% interval) and share of spiked InterPro terms found.

| Scenario | `fisher` FDP, GO | `fisher` FDP, InterPro | `covariate_enrichment` FDP, GO | `covariate_enrichment` FDP, InterPro | Power, `fisher` | Power, `covariate_enrichment` |
|---|---:|---:|---:|---:|---:|---:|
| N0s_0.02 | 1.5% | 2.5% | 1.5% | 2.5% | 0.3% | 0% |
| N0s_0.05 | 1.2% | 1.4% | 0.6% | 1.8% | 10.3% | 8.8% |
| N1s_0.02 | **28.0%** | **34.0%** | 1.0% | 0.5% | 0.8% | 0.6% |
| N1s_0.05 | **27.4%** | **20.2%** | 1.7% | 0.6% | 9.5% | 17.8% |

Under exchangeable labels, FDP is controlled for both. Under covariate-driven labels, Fisher's FDP is 20-34% and `covariate_enrichment`'s at most 1.7%. Power is low for every test at these spike sizes, because 2-5% of a term's genes is a small signal against 0.9% positives. It is reported for comparison only; the earlier study measured power at larger effects. `go_enrichment`'s defaults gave 4% FDP (interval to 9.5%) on 50 N0s_0.05 runs, which is inconclusive at that run count.

## Code audit (no runs)

- **Universe.** For `go_enrichment`, genes with at least one GO term of the ontology for that side, as topGO does; for `term_enrichment` and `ipr_enrichment`, genes with at least one term of the family; for `covariate_enrichment`, the same, also requiring dS > 0 and a known length. Positives are restricted to the universe. This matches the simulation's universes.
- **Family.** BH runs once across all requested GO ontologies (`adjust_scope = "global"`, the default) or across all InterPro entries. `ipr_enrichment(adjust_scope = "per_type")` makes each entry type its own family, which controls FDR within each type. The chance of some false term in the whole output then grows with the number of types (9.6% per run in [enrichment.md](enrichment.md)'s `final_per_type` rows). Results report the scope.
- **Filters.** `min_total` is applied before correction and does not depend on the outcome. `term_enrichment`'s `min_pos` is applied after correction, so p-values stay valid.
- **Gene selection.** Genes with dS = 0 or dN/dS ≥ `max_dnds` are excluded. That exclusion depends on length and dS, not on terms, so it acts only through the confounding measured here.
- **Optional `qvalue` and IHW** assume continuous null p-values. Hypergeometric p-values are discrete and stochastically larger under the null, so they tend to be conservative. They were not studied, and BH remains the default.

## Conclusions for the release

- **Calibration of `term_enrichment`, `ipr_enrichment`, `go_enrichment`: pass, for the hypothesis they test.** That hypothesis is whether the share of genes above the dN/dS threshold differs between a term and its universe. They do not test whether selection differs, and they should not be described as doing so.
- **For the scientific question, they fail, by confounding, at the size measured here and in [enrichment-confounding.md](enrichment-confounding.md).** The correction exists and is now validated twice, with estimator-noise labels and with exact labels: `covariate_enrichment`.
- **Recommendation, for the maintainer's decision:** keep the three commands, and document their result as an association of the dN/dS call with the term, warning (as every run already does) that length and dS can produce it. Present `covariate_enrichment` as the test for selection related to terms. Any change to defaults, or to which command the tutorials lead with, needs the maintainer's approval.

## Limits

One pair (*A. thaliana*-*A. lyrata*). Labels are independent across genes (no tandem families). The covariate effect is one smooth function. The earlier study's hard case, few positive genes in a divergent pair (*B. rapa*), was not repeated here. Small spikes only. `qvalue`, IHW and quantile thresholds were not tested.
