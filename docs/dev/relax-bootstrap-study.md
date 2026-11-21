# Study R: calibrating RELAX by parametric bootstrap (design)

Status: **pilot done; main study proposed for the maintainer's approval; not run.** RELAX's `p_value` (chi-square(1) on the LRT) is anticonservative on 16 taxa: 8.3-9.3% at 0.05 in every check of the follow-up ([report](../validation/phylo-lineage-calibration.md)); the LRT behaves like about 1.27 times chi-square(1) (audit, E2).

## Procedure under test

For each gene: RELAX as dndsR runs it, with `--save-fit`; the saved alternative fit loaded, `relax.K` constrained to 1, and refitted (the RELAX null); `B` alignments simulated from that fitted null by HyPhy (`SimulateDataSet`), so the bootstrap uses exactly RELAX's model and the gene's estimates; RELAX run on each with the same tree and settings. Scripts: `tools/dev/validation/relax-null-sim.bf`, `phylo-relax-boot.sb`, `phylo-relax-boot-summary.R`. HyPhy's own statistics are kept unchanged beside any new p-value.

Three p-values, all from the same replicates (so they cost nothing extra to compare):

1. **Monte Carlo p** (Besag-Clifford sequential): draw replicates until 10 have an LRT at least the observed one, or 99 are drawn; p = 10 / b at stopping, else (k + 1) / 100. Valid at any number of replicates; smallest possible value 0.01.
2. **Per-gene scaled LRT**: LRT divided by the mean of that gene's first 20 bootstrap LRTs, against chi-square(1). Continuous p-values (needed for BH over thousands of genes, which the Monte Carlo p cannot provide), but valid only if the scaled chi-square shape holds in the tail.
3. **Dataset-level scaled LRT**: LRT divided by the mean bootstrap LRT over a random 50 genes of the dataset (20 replicates each, 1,000 fits). Cheapest per gene; assumes one factor fits all genes of a dataset.

## Pilot (development genes; done 2026-11-21)

16 null genes of the first study's C0 T16 cell (already used, so development only), 20 replicates each, Slurm job of 16 cores:

- The refitted null equals RELAX's own null fit to within 0.026 log-likelihood units (median 5e-5); 0 of 320 replicate fits failed.
- **The bootstrap reproduces the inflation**: bootstrap LRTs averaged 1.25 (SE 0.10), quantiles 0.65 / 1.62 / 3.13 / 4.65 at 50/75/90/95% against chi-square(1)'s 0.45 / 1.32 / 2.71 / 3.84, and 7.2% exceeded 3.84. The simulated null genes themselves gave 1.27, and 4.86 at the 95th percentile. So the excess is a property of RELAX's model at this sample size, which HyPhy's own simulation shows; it is not an artifact of the study's simulator.
- Cost: about 44 s per gene for the observed fit, null refit, and simulation, and about 60 core-seconds per replicate fit.

## Questions and criteria (fixed before the main runs)

- **R5, size.** At 0.05, each p-value's size in held-out null cells. Criterion: within Bradley's tolerance (95% interval inside 2.5-7.5%) on T16 at 600 genes, and not anticonservative on T8b and T5. HyPhy's chi-square p is reported beside it, recomputed for the same genes.
- **R6, power.** Detection of k = 0.5 and k = 2 on T16. Criterion: each calibrated p-value's power at least the *size-adjusted* power of the chi-square test (its power at the critical value that gives 5% in the null cell) minus 5 points. Raw chi-square power is reported but is not the benchmark, since part of it comes from the excess false positives.
- **R7, direction.** Share of detections with k on the right side of 1, as in the first study.
- **R8, cost.** Core-hours per gene for each p-value, against the plain RELAX run.

## Cells (held-out, new seeds 401-405)

| Cell | Truth | Tree | Genes |
|---|---|---|---|
| R-N16 | k = 1 (C0) | T16 | 600 |
| R-N8 | k = 1 | T8b | 300 |
| R-N5 | k = 1 | T5 | 300 |
| R-K05 | k = 0.5 (C1) | T16 | 200 |
| R-K2 | k = 2 (C2) | T16 | 200 |

## Precision

600 null genes: an interval of about 3.4-7.1% at a true 5%. 200 genes per power cell: plus or minus 7 points near 50% power.

## Cost

Replicates needed: under the null, about 33 per gene on average (10 + 10 ln 9.9, for p uniform); for genes that are significant, the full 99. With 60 core-seconds per replicate on 16 taxa (about 30 on 8 and 20 on 5): R-N16 about 330 core-hours, R-N8 and R-N5 about 130 together, R-K05 and R-K2 about 470 (roughly half the genes significant): **about 930 core-hours**. A leaner version (stop at 5 exceedances; about 20 replicates per null gene) would cost about 600 and gives coarser p-values without affecting their validity; the proposal is the 930-core-hour version, because power and the per-gene scale estimate (20 replicates) both need it.

## What a positive result would and would not mean

Calibration in this simulation, for the RELAX configuration dndsR runs (Minimal, 3 rates). A bootstrap p-value would be a new inference procedure: if adopted, it would be an additional column, not a replacement for HyPhy's p-value, and only with the maintainer's approval.
