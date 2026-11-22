# RELAX calibration, stage 1: diagnosis and a genome-scale correction (study R)

**Question.** RELAX's `p_value` (chi-square(1) on the LRT) is anticonservative on 16 taxa ([phylo-lineage-calibration.md](phylo-lineage-calibration.md)). The maintainer approved a diagnostic stage before any expensive study: does a bootstrap correct it on independent held-out genes, where does the excess come from, can a correction support genome-wide BH, and does it keep power. Design, questions RD1-RD7 and the go/no-go rule, all fixed before running: [docs/dev/relax-bootstrap-study.md](../dev/relax-bootstrap-study.md), "Stage 1".

**Answer.** The excess comes from HyPhy itself, not from dndsR, and it shrinks as genes get longer. It is shared across the genes of a dataset, not specific to some of them. One shared scale factor per dataset, estimated from 1,000 null bootstrap fits on 50 of its genes, restores size on new held-out genes without losing power beyond the size-adjusted reference, and keeps BH's false discovery proportion within q. Both candidate p-values meet every criterion fixed in advance. **Recommendation: go for a roughly 60 core-hour confirmation on a further held-out configuration; no-go for the 930 core-hour per-gene Monte Carlo design**, which cannot support genome-wide BH by arithmetic.

## Runs

Simulated with the lineage simulator (cell C0: no relaxation, k = 1; C1 k = 0.5; C2 k = 2). RELAX through `phylo_relax()` with HyPhy 2.5.103. Code at `760ca03`, image `dndsr_v0.4.0_sha256-f1b42d9d`, units `tools/dev/validation/phylo-study-R-units.tsv` (seeds 401-405, 411-412, new). Calibration sets by `phylo-relax-calib.sb`: 50 tested genes per dataset, sampled with a seed from the unit's name; 20 replicates each, simulated by HyPhy from that gene's fitted RELAX null and refitted (1,000 null fits per dataset). The bootstrap replicates never use the observed LRTs they correct. Summary: `phylo-study-R-summary.R`; tables in `data/relax_stage1_*.tsv`.

| Unit | Truth | Taxa | Genes | Codons |
|---|---|---:|---:|---:|
| R-N16 | k = 1 | 16 | 600 | 300 |
| R-N8 | k = 1 | 8 | 300 | 300 |
| R-N5 | k = 1 | 5 | 300 | 300 |
| R-K05 | k = 0.5 | 16 | 200 | 300 |
| R-K2 | k = 2 | 16 | 200 | 300 |
| RD-L1000, RD-L3000 | k = 1 | 16 | 300, 200 | 1,000, 3,000 |

## Diagnosis

**Wrapper against upstream HyPhy (settled earlier, rechecked here).** HyPhy 2.5.103 run directly on dndsR's inputs for R-N16 reproduced dndsR's LRTs in all 599 genes (largest difference 1e-14). dndsR introduces nothing.

**RD1, version.** HyPhy 2.5.104 gives identical LRTs to 2.5.103. HyPhy 2.5.62 gives a similar excess (8.7%, interval 6.6-11.2, against 8.2%), but individual LRTs agree within 0.01 in only 34% of genes (largest difference 5.1). Criterion: the share moves by less than 2 points (met), and LRTs agree in at least 90% of genes (not met). So the version changes individual optima but not the inflation, which is present in all three versions. dndsR pins 2.5.103, so its results are reproducible within that version.

**RD2, finite sample.** Mean null LRT, chi-square(1) expects 1:

| Codons | Genes | Mean LRT (se) | Share above 3.84 (95% CI) |
|---:|---:|---:|---:|
| 300 | 599 | 1.29 (0.07) | 8.2% (6.1-10.7) |
| 1,000 | 300 | 1.03 (0.08) | 4.3% (2.3-7.3) |
| 3,000 | 200 | 1.07 (0.10) | 6.5% (3.5-10.9) |

The mean falls below 1.1 at 3,000 codons, which meets the criterion for a finite-sample effect. The share at 3,000 codons (6.5%) is not separable from either 5% or 8% with 200 genes, so the evidence comes from the mean, not the share. The inflation is a small-sample property of the asymptotic test at typical gene lengths, not a model or software defect.

**Optimizer and simulator misspecification** were settled before this stage (best of several starts still 9.3%; data simulated by HyPhy from RELAX's own null reproduce it).

**RD3, shared or per-gene.** Variance of per-gene mean bootstrap LRTs against replicate noise alone (10,000 permutations):

| Dataset | Mean bootstrap LRT, c | Share above 3.84 | p, heterogeneity |
|---|---:|---:|---:|
| R-N16 | 1.25 | 8.5% | 0.84 |
| R-N8 | 1.34 | 8.7% | 0.084 |
| R-N5 | 1.20 | 6.5% | 0.050 |
| R-K05 | 1.27 | 7.4% | 0.28 |
| R-K2 | 1.31 | 8.6% | 0.18 |

None reaches p < 0.01, so a shared factor per dataset is supported. R-N5 and R-N8 are close enough to the threshold that heterogeneity among genes of a small tree cannot be excluded at a larger size.

## Genome-scale procedures

Two p-values, both continuous enough for BH, computed for every gene from the dataset's calibration set:

- `p_scaled`: P(chi-square(1) > LRT / c).
- `p_pool`: (1 + number of the 1,000 bootstrap LRTs at least the observed) / 1,001; beyond the largest, the scaled chi-square tail times 1,000 / 1,001.

**RD4, held-out size at 0.05.** Criterion: on R-N16 the 95% Clopper-Pearson interval inside Bradley's 2.5-7.5%; on R-N8 and R-N5, lower bound not above 5%.

| Dataset | HyPhy `p_value` | `p_scaled` | `p_pool` |
|---|---:|---:|---:|
| R-N16 | 8.2% (6.1-10.7), anticonservative | **5.2% (3.5-7.3), within** | **5.2% (3.5-7.3), within** |
| R-N8 | 7.3% (4.7-10.9) | 2.7% (1.2-5.2), not anticonservative | 5.3% (3.1-8.5), not anticonservative |
| R-N5 | 7.0% (4.4-10.5) | 4.7% (2.6-7.7), not anticonservative | 4.7% (2.6-7.7), not anticonservative |

Both pass. `p_scaled` is conservative on R-N8, where the bootstrap tail is lighter than a scaled chi-square.

**RD5, held-out power at 0.05.** Criterion: at least the size-adjusted power of HyPhy's chi-square (its power at the critical value that gives 5% in R-N16) minus 5 points. Each power cell used its own calibration set, which in a real genome would include selected genes.

| Cell | HyPhy `p_value` (unadjusted) | Size-adjusted reference | `p_scaled` | `p_pool` | Right direction |
|---|---:|---:|---:|---:|---:|
| R-K05 (k = 0.5) | 64.5% | 51.5% | 51.0% | 52.5% | 98% |
| R-K2 (k = 2) | 31.0% | 21.0% | 21.0% | 21.0% | 98% |

Both pass. The loss against HyPhy's unadjusted figure is the size it had wrongly gained, not power.

**RD6, genome-wide BH.** 2,000 resampled 660-gene mixtures: all 600 R-N16 nulls plus 30 genes each from R-K05 and R-K2, all scored against R-N16's calibration set. BH at q = 0.05; mean (95% range) false discovery proportion and power:

| p-value | FDP | Power |
|---|---:|---:|
| HyPhy `p_value` | 5.2% (0-14.3) | 9.2% (0-20.0) |
| `p_scaled` | 0% (0-0) | 3.2% (0-8.3) |
| `p_pool` | 0% (0-0) | 3.2% (0-8.3) |

Both keep FDP within q (criterion met). HyPhy's own p-value already exceeds q on average. Power under genome-wide BH is low for every p-value at these effect sizes and 300 codons: RELAX finds strong relaxation in single genes, but a genome-wide scan at 9% alternatives will report few genes. That limitation belongs to the test at this gene length, not to the correction.

**Monte Carlo per gene** (the main design): smallest possible p is 1 / (B + 1); BH over 10,000 genes needs about 5e-6, about 200,000 replicates per gene. Not usable for genome-wide BH, as fixed in the design.

**RD7, cost.** Plain RELAX on R-N16: 49 core-seconds per gene, about 136 core-hours for 10,000 genes. A calibration set is 1,000 fits regardless of genome size: 37-44 core-hours here. That is about 30% on top of a 10,000-gene run, against about 4,000 core-hours per 1,000 genes for per-gene Monte Carlo at B = 100 (130-160 core-seconds per bootstrap fit).

## Go / no-go, under the rule fixed in advance

- **Main study as designed (per-gene Monte Carlo, ~930 core-hours): no-go**, by RD6's arithmetic. Holding it for a dataset's few discoveries as a per-gene check remains possible.
- **Genome-scale procedure: go.** `p_scaled` and `p_pool` both pass RD4 on R-N16 without being anticonservative on R-N8 or R-N5, meet RD5 in both power cells, keep RD6's FDP within q, and RD3 shows no heterogeneity at p < 0.01.
- **Next allocation, for approval: about 60 core-hours**, a confirmation on a further held-out configuration (another tree and omega distribution), evaluating both p-values from the same fits, with criteria fixed before it runs. Any column added to `phylo_relax()` output would sit beside HyPhy's `p_value`, not replace it, and needs the maintainer's approval as an interface change.

## Limits

- One simulator, two omega configurations, three trees, 300 codons in the size and power cells. Real alignments with indels, misassembly or recombination were not part of this stage (study B covers input quality).
- c is estimated from 50 genes (standard error about 0.05, so about 4% in c). It is dataset-specific and must be estimated on each dataset, at the cost above.
- The calibration genes were simulated from their own fitted nulls; if a dataset's genes are mostly under relaxation, its fitted nulls are misspecified in a way this stage did not test, beyond the two power cells.
- RD3's p-values for the 5- and 8-taxon trees (0.05, 0.08) are close enough to the threshold that a larger study could find heterogeneity there.
- HyPhy 2.5.62 moves individual LRTs; any correction is valid only for the HyPhy version that produced its calibration set.
