# Phylogenetic module: calibration stage 2A

Simulated data only. Plan, fixed before any stage 2 job: [phylo-calibration-stage2.md](../dev/phylo-calibration-stage2.md). This is stage 2A, the bounded intermediate batch. It reports and applies the plan's rules for stage 2B. Recommendations are for the maintainer to decide on. **No default has changed**, and nothing here is added to the package.

Tables: `docs/validation/data/phylo_stage2a_*.tsv` (`units`, `paired`, `starts_summary`, `filtering`, `cost`). Unit names read `block__scenario__config__filter__seeds`.

## Run

- 108 units, 21,270 BUSTED fits, all finished. **0 genes failed and 0 were excluded.** 635 core-hours.
- Commit `5020a91`, installed into a read-only library; development image `sha256:e24b394df53597922308224fd2507956fb6f70bd146e470540b1c91179f1d17f` (HyPhy 2.5.103, not published), checked by every job. Every fit ran on AMD EPYC 9654.
- Inputs: `params.json` md5 `6a678fadefcbcf6b4d6d6377505997a9`; empirical segment lengths md5 `30766db96c312d825e0c49e331db23bc`; simulation seeds in `tools/dev/validation/phylo-stage2-scenarios.tsv`; HyPhy seeds 1, 2, 3.
- Configurations: C0 BUSTED (dndsR's default), C1 `--multiple-hits Double`, C2 `--error-sink Yes` (BUSTED-E), C3 BUSTED-E with `Double`, C4 `Double+Triple`. All with `--srv Yes`. Configuration and filter blocks use one start (seed 1), as dndsR reports by default.

## Conclusions

1. **Every configuration is calibrated under correctly specified nulls.** Across 4 null cells of 300 genes, the rate at α = 0.05 is 0-1.3%, with every 95% upper bound at or below 3.4%. C2 and C3 made no call at all in 1,200 null genes each.
2. **Each protection costs power, and BUSTED-E costs a lot.** At ω+ = 10 on 5% of sites, power is 39% (C0), 34% (C1), 11% (C2), 8% (C3). At ω+ = 10 on 10%: 67%, 61%, 43%, 39%. Double hits cost little except at ω+ = 10 on 20% (91% to 70%).
3. **Multinucleotide changes (M1) need double hits.** False positives: C0 97%, C2 56%, C1 0.7%, C3 0%, C4 0%.
4. **BUSTED-E only partly protects against realistic alignment corruption; masking flagged windows nearly removes it.** One empirical-length segment in At: C0 91-92%, C1 89-91%, C2 and C3 24-28%. After masking the flagged windows: 0-0.7% in every configuration tested (C0, C3). The same masking cut the design's random misaligned window (M3) from 87% to 0.7%.
5. **Masking cost no detectable power on clean alternatives.** It changed 0 calls in clean A(10, 0.1) and A(10, 0.2), under C0 and C3. With selection plus corruption, it restored power to the clean level (C0: 68% masked, 67% clean; C3: 38% and 39%), removing 36 (C0) and 48 (C3) calls, most of them added by the corruption.
6. **Corruption in a background taxon (Br) caused no false positives** in any configuration, as expected for a foreground test.
7. **Realigned indels (M4) inflate C0 and C1 modestly, and masking does not help.** C0 8.0% (4.2-13.6%), C1 8.7%; C2 and C3 0%. The window rule flagged the foreground in only 1.3% of M4 genes. Dropping gappy columns as well changed 6 calls (net 2 fewer, McNemar p 0.69).
8. **One optimizer start is not enough under selection, especially for BUSTED-E.** With the best alternative and best null fits over 3 seeds, C3's power at A(10, 0.1) rose from 39% to 53% (21 calls changed, 19 gained). For C0 it rose from 67% to 69%. Going from 2 to 3 starts still changed up to 2.5% of calls. Extra starts changed no null call (0 of 600), and moved M1 under C0 from 97% to 99%.
9. **Triple hits did not differ from double hits in net rate.** C4 and C1 had the same rates in every cell, but 13-14% of calls swapped in each power cell (9 and 7, 8 and 9; McNemar p 0.8 and 1). That swap rate is the size of the start-to-start noise in finding 8, not a directional effect.
10. **Cost per fit (core-seconds, null L50 / power):** C0 42 / 97, C2 48 / 120, C1 78 / 169, C3 111 / 235, C4 131 / 318.

## Results

### Calibration and power, one start, α = 0.05

Rate, with Clopper-Pearson 95% interval. Genes: 300 per null cell, 120 per alternative, 150 per misspecified or corrupted cell.

| Scenario | C0 | C1 | C2 | C3 |
|---|---|---|---|---|
| N0 L50 | 0.7% (0.1-2.4) | 0.3% (0-1.8) | 0% (0-1.2) | 0% (0-1.2) |
| N0 L25 | 1.3% (0.4-3.4) | 1.3% (0.4-3.4) | 0% | 0% |
| N0 L50, divergence ×2 | 0.3% | 0% | 0% | 0% |
| N1 L50 | 0.3% | 0% | 0% | 0% |
| A(3, 0.2) | 12% (6-19) | 10% | 4% | 4% |
| A(5, 0.1) | 24% (17-33) | 23% | 8% | 8% |
| A(10, 0.05) | 39% (30-49) | 34% | 11% | 8% |
| A(10, 0.1) | 67% (58-75) | 61% | 43% | 39% |
| A(10, 0.2) | 91% (84-95) | 70% | 65% | 59% |
| M1, multinucleotide, δ = 0.2 | 97% | 0.7% | 56% | 0% |
| M3, random misaligned window | 87% | 84% | 9% | 8% |
| M4, realigned indels | 8.0% | 8.7% | 0% | 0% |
| R1, internal segment in At | 92% | 89% | 27% | 24% |
| R2, terminal segment in At | 91% | 91% | 27% | 28% |
| R3, frameshift segment in At | 92% | 89% | 24% | 24% |
| R4, internal segment in Br | 0% | 0% | 0% | 0% |
| A(10, 0.1) with R1 | 98% | 98% | 75% | 77% |

C4: N0 L50 0%, M1 0%, A(10, 0.1) 59%, A(10, 0.2) 71%.

Paired comparisons on the same genes (`phylo_stage2a_paired.tsv`): every power loss of C2 or C3 against C0 has McNemar p below 0.004, and every gain against corruption p below 1e-11. C3 and C2 differ only in M1 (84 calls) and within noise elsewhere (p ≥ 0.07).

### Filtering (flagged windows masked: F1; also gappy columns dropped: F2)

| Scenario | Flagged foreground | C0 none | C0 F1 | C3 none | C3 F1 |
|---|---:|---:|---:|---:|---:|
| N0 L50 | 0% | 0.7% | 0.7% | 0% | 0% |
| A(10, 0.1) | 0% | 67% | 67% | 39% | 39% |
| A(10, 0.2) | 5% | 91% | 91% | 59% | 59% |
| M3 | 99% | 87% | 0.7% | 8% | 0.7% |
| M4 | 1.3% | 8.0% | 8.0% | 0% | 0% |
| R1 | 100% | 92% | 0.7% | 24% | 0% |
| R2 | 100% | 91% | 0.7% | 28% | 0% |
| R3 | 100% | 92% | 0% | 24% | 0% |
| R4 | 0% (Br 100%) | 0% | 0.7% | 0% | 0% |
| A(10, 0.1) with R1 | 100% | 98% | 68% | 77% | 38% |

F2 equals F1 everywhere except M4, where no window was masked but 22 gappy columns per gene were dropped on average: C0 8.0% to 6.7%. On clean alternatives, the flag marked no foreground sequence except at A(10, 0.2), 5% of genes, where strong selection itself looks dissimilar; any sequence was flagged in 0.8-1.7% of genes, and 6.7% at A(10, 0.2). Masking those windows changed no call.

### Optimizer starts

Best alternative and best null fits over seeds 1, 2, 3 (never the largest LRT).

| Scenario, configuration | 1 start | 2 | 3 | Calls changed, 1 vs 3 | 2 vs 3 | Alternative fit improved > 0.5 by 3 starts |
|---|---:|---:|---:|---:|---:|---:|
| N0 L50, C0 | 0.7% | 0.7% | 0.7% | 0 / 300 | 0 | 7% |
| N0 L50, C3 | 0% | 0% | 0% | 0 / 300 | 0 | 10% |
| M1, C0 | 96.7% | 97.3% | 98.7% | 3 / 150 | 2 | 32% |
| M1, C3 | 0% | 0% | 0% | 0 / 150 | 0 | 2% |
| A(10, 0.05), C0 | 39% | 38% | 38% | 9 / 120 | 1 | 32% |
| A(10, 0.05), C3 | 8% | 10% | 9% | 7 / 120 | 3 | 41% |
| A(10, 0.1), C0 | 67% | 69% | 69% | 5 / 120 | 0 | 26% |
| A(10, 0.1), C3 | 39% | 53% | 53% | 21 / 120 | 2 | 57% |

HyPhy's reported LRT differed from twice its final likelihood difference in 50 of 1,380 start-1 fits, 47 of them by less than 0.1. In the other three it reported LRT 0 after a refit while the likelihoods gave 0.39, 0.57, and 2.10. In the last (A(10, 0.05), C3), the best of 3 starts gave LRT 4.79, p 0.046.

### Cost

Core-seconds per fit, from wall time and threads, one processor class (`phylo_stage2a_cost.tsv`). Relative to C0, by scenario kind: C1 1.5-2.1×, C2 1.1-1.6×, C3 2.2-2.9×, C4 1.8-3.3×. Three starts cost three fits.

## The plan's rules for stage 2B

1. **Redundant configurations: none.** Every configuration differs from each cheaper one by more than 1% of calls in some cell, with paired p below 0.001.
2. **Unsafe configurations: none.** No clean null cell has a 95% lower bound above 0.05 (the largest is 0.4%).
3. **Starts: k = 3.** Going from 2 to 3 starts changed more than 1% of calls in 3 of 8 cells (up to 2.5%). Whether 3 is enough is not shown.
4. **Triple hits: met, but only by swaps.** C4 changed 13-14% of calls against C1 in both power cells, so the rule allows scaling it. The changes are balanced and of the size of start noise (finding 9), so they do not show a C4 effect. C4 did not lower M1 meaningfully (1 of 150 to 0).
5. **Filtering: F1 met.** It lowered false positives in every corrupted foreground cell and in M3. F2 meets the rule only through its F1 step; its own step changed no corrupted cell, and changed M4 within noise.

## What stage 2A does not settle, and why 2B needs a decision

- **The rules as written exceed the budget.** Keeping C0-C3 with 3 starts across the full study's scope costs very roughly 12,000 core-hours, against about 3,100 approved (4,000 with reruns), of which 2A used 635, leaving about 2,450. Scaling therefore needs choices that the rules did not fix: which configurations to carry, and where one start is enough.
- **Rule 1 and rule 4's 1% threshold is below the optimizer's own noise.** In power cells, one start and three starts of the same model differ in 4-17.5% of calls. Paired thresholds at 1% will flag noise as a difference.
- **Untested combinations:** C1 with masking, the likely candidate for routine use, was not run. Filtering ran only with C0 and C3.
- **M4 remains open:** realigned indels inflate C0 and C1 (8%), and masking does not catch them.
- **The 10-taxon tree's parameters are not yet defined.**

### Proposed stage 2B, for the maintainer's decision

About 2,400 core-hours, within what remains, estimated from 2A's costs:

| Block | Cells | Genes | Configurations, starts | Core-hours |
|---|---|---:|---|---:|
| Null grid | N0, N1 × L25, L50, L75 × divergence ×0.5, 1, 2 | 500 | C0, C1, C3; 1 start (2A: 0 of 600 null calls changed with 3) | 580 |
| Starts on nulls | N0 L50 | 1,000 | C1, C3; seeds 2, 3 | 110 |
| Power grid | ω+ 2, 3, 5, 10 × p+ 0.02, 0.05, 0.1, 0.2, L50 | 150 | C1, C3 with 3 starts; C0 with 1 | 870 |
| Misspecification | M1 δ 0.1, 0.3; M2; M5 | 300 | C1, C3 | 90 |
| Filtering | R1-R3, M3, M4, and A(10, 0.1) with R1 | 300 | C1, with and without F1 | 110 |
| Triple hits | M1 δ 0.3; A(10, 0.1), A(10, 0.2) | 300 | C4 with 3 starts | 190 |
| Genome-wide BH | 10 sets of 500, corruption at the P3A per-taxon rates | 500 | C1 with F1, 3 starts | 460 |
| **Total** | | | | **about 2,400** |

It departs from the plan in four places, each needing approval: it drops C2 from scaling (it does not protect against M1, 56%, and in every other cell is within noise of C3, which does), against rule 1 as written; it uses one start in null cells, against rule 3; it halves null cells to 500 genes (a 0.5% rate then has an interval of 0.1-1.7%) and the genome-wide sets to 10; and it leaves the 10-taxon tree for later. Nothing in 2B is launched before the maintainer decides.

## Preliminary recommendation (to be confirmed by 2B)

For the maintainer, not a default change. Masking is a validation script here, not a dndsR feature. Adding it would be a separate, approved change.

- **Default (C0) stays honest only on clean data.** It is calibrated under correct nulls and has the most power, but gives 91-97% false positives under multinucleotide changes or one corrupted foreground segment.
- **Exploratory scans:** `--multiple-hits Double` (C1), with the foreground's flagged windows masked, 3 starts combined by likelihood (`--extra-seeds 2,3`, reading `lrt_best_fits`), BH within the scan. In 2A, C1 handled M1 and lost little power except at very strong selection. C1 with masking was not tested directly; C0 with masking was.
- **Conservative genome-wide scans:** as above, and report a gene as robust only if it is also significant under BUSTED-E with double hits (C3). That removes the M4 inflation and most residual corruption, at roughly half the power.
- **Not recommended:** BUSTED-E without double hits (C2), which fails under M1.

## Limits

- The simulator is close to the model BUSTED fits, with five taxa, one tree, κ = 2, and one foreground clade.
- Corrupted cells carry exactly one segment per gene. Real genes may carry none or several.
- Cell sizes detect gross differences only: ±5% at 150 genes near 90%, about ±9% at 120 genes near 50%.
- Configuration and filter results use one start. C3's power in particular is underestimated by one start (finding 8).
- Which share of real flagged segments are errors, rather than true divergence, is not known. Simulation knows the truth, real data do not.
