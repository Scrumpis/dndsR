# Study R: calibrating RELAX by parametric bootstrap (design)

Status: **pilot done; stage 1 done ([results](../validation/relax-calibration-stage1.md)): no-go for the ~930 core-hour main design; both genome-scale p-values pass, and a ~60 core-hour confirmation is recommended, not yet approved.** RELAX's `p_value` (chi-square(1) on the LRT) is anticonservative on 16 taxa: 8.3-9.3% at 0.05 in every check of the follow-up ([report](../validation/phylo-lineage-calibration.md)); the LRT behaves like about 1.27 times chi-square(1) (audit, E2).

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

## Stage 1: diagnosis and a genome-scale procedure (approved 2026-11-21; fixed before running)

The maintainer approved R in stages: first establish that a bootstrap correction works on independent held-out genes, diagnose the discrepancy against upstream HyPhy, and show that a correction can support genome-scale multiple-testing at feasible cost, then recommend go or no-go for the expensive phase.

### What existing evidence already settles

- **Wrapper against upstream HyPhy, identical inputs:** settled. HyPhy run directly on the simulator's FASTA files and a hand-written tree gave 9.0 / 9.0 / 8.7% at three seeds, against dndsR's 9.0% on the same genes; on dndsR's own inputs, LRTs agreed ([report](../validation/phylo-lineage-calibration.md), follow-up part R). The excess is not introduced by dndsR.
- **Optimizer:** largely settled. The best of several seeds per gene gave 9.3% (a better optimum does not remove it); the null refitted from RELAX's own alternative fit agrees with RELAX's null to 0.026 log-likelihood units (pilot); 14 of 300 genes moved by more than 1e-3 between identical runs, too few to explain an excess of about 4 points.
- **Simulator misspecification:** settled. Data simulated by HyPhy from RELAX's own fitted null reproduce the inflation (pilot: mean LRT 1.25, 7.2% above 3.84), and per-branch site classes, which RELAX's random-effects model assumes, gave 8.5% (C0i).

### Questions

- **RD1, version.** Is the inflation specific to HyPhy 2.5.103? R-N16's dndsR inputs rerun directly with HyPhy 2.5.104 (latest) and 2.5.62 (older), same settings. Criterion: the share above 3.84 differs from 2.5.103's by less than 2 points, and LRTs agree to within 0.01 in at least 90% of genes; otherwise the version matters and is reported.
- **RD2, finite sample.** Does the excess shrink with more data per gene? C0 on T16 at 1,000 codons (300 genes, seed 411) and 3,000 codons (200 genes, seed 412), against 300 codons. Read on the mean LRT (standard error about 0.10 at 300 genes): finite-sample if the mean falls toward 1 (to below 1.1 at 3,000 codons), and the share above 3.84 with it; otherwise a persistent property of the model on these trees.
- **RD3, shared or per-gene.** Is the inflation the same for every gene of a dataset? From each dataset's calibration set (50 genes, 20 replicates each): the variance of per-gene mean bootstrap LRTs against the variance expected from replicate noise alone (permutation of replicates across genes, 10,000 permutations). Heterogeneous if p < 0.01 and the excess variance exceeds the replicate variance; then a shared factor cannot be exact, and per-gene features (alignment length, tree, a phantom omega class) are examined.
- **RD4, held-out size of the genome-scale procedures.** On the held-out datasets R-N16 (600 genes, seed 401), R-N8 (300, 402), R-N5 (300, 403): each dataset's calibration set is 50 of its tested genes, sampled by a seed from the unit's name, with 20 replicates each (1,000 null fits). Two p-values, computed for every gene, both continuous enough for BH:
  1. *Shared scaled chi-square*: p = P(chi-square(1) > LRT / c), c = the mean of the 1,000 bootstrap LRTs.
  2. *Pooled bootstrap null*: p = (1 + number of the 1,000 bootstrap LRTs at least the observed) / 1,001, resolution about 0.001, and beyond the largest bootstrap LRT the scaled chi-square's tail times 1,000 / 1,001.
  Criteria as R5: on R-N16, size within Bradley's tolerance (95% Clopper-Pearson interval inside 2.5-7.5%); on R-N8 and R-N5, not anticonservative (lower bound not above 5%). HyPhy's own p-value is reported beside them, unchanged.
- **RD5, held-out power.** R-K05 (k = 0.5, 200 genes, seed 404) and R-K2 (k = 2, 200, seed 405), each with its own calibration set (in a real genome the calibration genes include selected ones). Criterion as R6: each procedure's power at least the size-adjusted power of HyPhy's chi-square (its power at the critical value giving 5% in R-N16) minus 5 points. Direction (R7) reported.
- **RD6, multiple testing at genome scale.** The Monte Carlo p-value of the main design has a smallest value of 1 / (B + 1): BH over m genes needs p below 0.05 / m for its first rejection, which for 10,000 genes is 5e-6 and would need about 200,000 replicates per gene. So the Monte Carlo p-value is not usable for genome-wide BH, by arithmetic, and is kept only as a per-gene check. RD6 measures, in 660-gene mixtures formed from all 600 R-N16 genes and 30 each of R-K05 and R-K2 (9% alternatives; all with R-N16's calibration set, as a genome would have one), the false discovery proportion and power of BH at q = 0.05 on each continuous p-value, with 2,000 resampled mixtures (the alternatives resampled) to give their spread.
- **RD7, cost.** Core-hours of the calibration set against the plain RELAX run, projected to 10,000 genes.

### Go / no-go (fixed now)

- **No-go for the main study as designed** (per-gene Monte Carlo on every gene) whatever stage 1 shows, because of RD6's arithmetic: it cannot provide genome-wide BH. A smaller per-gene confirmation for a dataset's discoveries remains possible.
- **Recommend a genome-scale procedure for approval** if one of RD4's p-values passes on R-N16 and is not anticonservative on R-N8 and R-N5, meets RD5 on both power cells, keeps RD6's FDP within q (point estimate at most 5%), and RD3 shows no heterogeneity that makes it fail in a subgroup. It would be an additional column beside HyPhy's p-value, not a replacement, and only with the maintainer's approval; then a confirmation on a further held-out configuration (another tree and omega distribution, about 60 core-hours) is the next allocation, not 930.
- **If none passes**, RELAX stays flagged, its use is limited to descriptive k estimates, and the report says that calibration could not be established.

### Cost

Plain RELAX on the held-out cells: about 15 core-hours; RD2: about 36; calibration sets (5 datasets x 1,000 fits): about 64; RD1: about 8. **About 125 core-hours.** Units: `phylo-study-R-units.tsv` (seeds 401-405, 411-412) and a 10-gene smoke copy (`phylo-study-R-pilot-units.tsv`, seeds 491-497); calibration sets by `phylo-relax-calib.sb`.
