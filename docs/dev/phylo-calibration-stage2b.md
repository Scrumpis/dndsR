# Phylogenetic module: calibration, stage 2B (proposed design)

**Status: approved by the maintainer on 2026-11-14 for the essential blocks B1-B7 (about 1,566 core-hours, up to about 1,806 if the B1 audit requires three starts in B2), with the clarifications below. No optional block is approved.** Changes after approval go in the change log at the end, with their reason.

This revises the stage 2B scope in [phylo-calibration-stage2.md](phylo-calibration-stage2.md), following the maintainer's priorities after [stage 2A](../validation/phylo-calibration-stage2a.md). It replaces that page's "Rules for stage 2B", whose 1% paired threshold 2A showed to be below optimizer noise. Everything else in the stage 2 plan stands: the simulator, parameters, corruption profile, filters F1 and F2, measures, and reproducibility.

The goal is a reliable, reproducible, defensible selection-testing workflow for dndsR, not the largest number of fits. No package default changes, and masking is not added to the package.

## Maintainer's clarifications at approval

- Masking: all-sequence masking (F1) is primary, as in 2A. Foreground-only masking (F1fg) is compared wherever it gives a different alignment. False-positive rates, power, and the fraction of positions removed are reported for both, and no scope is recommended until the evidence supports one.
- Masking's power cost is judged by non-inferiority (D6, D7), with paired effect sizes and intervals; weak and clustered selection get particular attention.
- Stage 2B's replicates are simulated independently of 2A's, with their own documented seeds. Reusing fits of identical alignments is acceptable only for paired methodological comparisons.
- Model averaging is evaluated with three starts, including whether better BUSTED-E optimization changes the Akaike weights, power, and false-positive rates.
- B6 also asks whether the gap-flank mask removes genuine selection along with alignment artifacts.
- B7 reports the false discovery proportion of one simulated mixture, with its uncertainty, not an expected genome-wide FDR.
- Optimization diagnostics, per-gene results, software versions, alignment checksums, and simulation metadata are kept. No dndsR default or feature changes. B2 may launch automatically by the B1 audit rule, after the audit and decision are recorded here.

## Budget

| | Core-hours |
|---|---:|
| Approved for the study | about 3,100 |
| Used by stage 2A | 635 |
| Remaining | about 2,465 |
| **This design: essential blocks** | **about 1,570** |
| Contingent: null grid with 3 starts, if the audit fails | +240 |
| **Committed at most** | **about 1,810** |
| Left uncommitted, for follow-up and reruns | at least 650 |
| Optional blocks, each needing its own go-ahead | 42-354 each |

Estimates use 2A's measured core-seconds per fit, by configuration and scenario kind, on the same processor class, scaled by gene length where it differs. 2A ran every fit with 0 failures, so no failure margin is added beyond the uncommitted reserve. Expect ±20%.

## What 2A leaves open

1. Double hits with masking (C1 + F1), the likely routine workflow, was never run.
2. Realigned indels (M4) inflate C0 and C1 (8%), and the window rule does not flag them.
3. Whether masking removes genuine selection when selection and corruption coexist, at weak as well as strong selection, and when selection is concentrated in one region.
4. Whether one optimizer start is enough for null cells, including difficult ones.
5. Whether Double+Triple differs from Double once optimization is stabilized.
6. How agreement between BUSTED and BUSTED-E compares with the model averaging of the BUSTED-E paper.
7. False discovery after BH in a genome-like mixture.

**Model averaging, preview from 2A.** Selberg et al. (2024, bioRxiv 10.1101/2024.11.13.620707, "Minus the Error"; PMC11601313) fit BUSTED and BUSTED-E to the same data, compute Akaike weights from AICc, and combine p-values: p(MA) = w(BUSTED) p(BUSTED) + w(BUSTED-E) p(BUSTED-E). It needs no extra fits: each alternative fit's AICc is in the HyPhy JSON that 2A kept. Applied to 2A's one-start fits (C1 with C3), MA tracked BUSTED closely: clean nulls 0-1.3%; A(10, 0.1) 62% (C1 61%, C3 39%); R1-R3 85-88% (C1 89-91%, C3 24-28%); M3 83%; M4 8.7%. Agreement (significant under both) tracked C3. In this simulation, then, MA recovers power but not protection from corruption. It is a one-start preview only: C3's fits are the least converged (2A, finding 8), which lowers its likelihood and so its weight.

## Procedures compared

| | Model | Alignment | Role |
|---|---|---|---|
| C1 | `--multiple-hits Double` | as prepared | routine candidate |
| C1+F1 | C1 | flagged windows masked | routine candidate |
| C1+F1fg | C1 | flagged windows masked in foreground sequences only | variant of F1 |
| C1+F3 | C1 | F1, and codons within 3 of an internal gap masked | indel candidate |
| C3 | `--error-sink Yes --multiple-hits Double` | as prepared | conservative candidate |
| C3+F1 | C3 | masked | conservative candidate |
| C1∧C3 | significant under both | | conservative rule |
| MA(C1, C3) | Akaike-weighted p, as Selberg et al. | | conservative rule |
| C4 | `--multiple-hits Double+Triple` | as prepared | targeted |
| C0 | BUSTED, dndsR's default | as prepared | baseline, power cells only |

- **Dropped:** C2 (BUSTED-E without double hits). It does not protect against multinucleotide changes (56% in 2A) and was within noise of C3 elsewhere. F2 (gappy columns dropped) is replaced by F3, which targets the indel problem directly.
- **Masking uses only what a real analysis has.** F1 applies dndsR's own dissimilar-window rule to the alignment. F3 uses only gap positions in the alignment. Neither reads the simulator's records of which codons were corrupted or selected. The true alignment of indel cells is used only as a diagnostic control, never as a filter.
- **F1 or F1fg.** F1 masks flagged windows in every sequence, as in 2A. Masking background sequences had no effect on false positives in 2A (R4: 0% unmasked), but it removes data. F1fg masks only the foreground's flagged windows. Both are run where they differ, and the recommendation names one.

## Optimizer, likelihoods, and records

- **Three starts (HyPhy seeds 1, 2, 3)** in every cell with selection, every misspecified or corrupted null, and the null audit. Combined by the best alternative and the best null likelihoods, never the smallest p-value or largest LRT. This is how dndsR's `lrt_best_fits` works.
- **One start for the clean null grid only if the audit passes** (D4 below).
- **Noise floor:** a second, independent three-start replicate (seeds 4, 5, 6) of C1 in two power cells measures how many calls differ between two stabilized runs of the same model. Differences between procedures are judged against it, never against a fixed 1%.
- **Likelihood consistency:** for every fit, the reported LRT is compared with twice the likelihood difference. Counts above 0.01 and 0.5 are reported by procedure, as are refits after negative LRTs. Inference uses the likelihood-based LRT.
- **Preserved:** every start's HyPhy JSON and log, every alignment (simulated, corrupted, filtered), the corruption and filter logs, the unit specifications, and the summary tables.

## Matrix

3 starts unless stated. Genes per cell in brackets. "Reuse" means 2A fits of byte-identical alignments with the same options and seed (see Reproducibility).

### Essential

| Block | Scenarios [genes] | Procedures | Core-h |
|---|---|---|---:|
| B1. Null audit for one start | N0 L25, N0 L50 div ×2, N1 L50, M1 δ 0.2, R1 with F1 [300 each] (M4 in B6) | C1, C3 | 171 |
| B2. Null grid | N0 × L25, L50, L75 × divergence ×0.5, 1, 2 [300]; 2 cells are B1's | C1, C3; 1 start if D4 passes | 104 |
| B3. Power, stratified | Weak: A(2, 0.2), A(3, 0.05), A(3, 0.1), A(3, 0.2). Moderate and strong: A(5, 0.1), A(5, 0.2), A(10, 0.05), A(10, 0.1), A(10, 0.2) [150] | C1, C3; F1 refits where any window is flagged | 374 |
| B3b. Baseline | A(3, 0.2), A(10, 0.05), A(10, 0.1), A(10, 0.2) [150] | C0 | 23 |
| B4. Triple hits, stabilized | A(10, 0.1), A(10, 0.2) [150]; M1 δ 0.3 [150]; noise floor: C1 seeds 4-6 on both power cells | C4, C1 | 131 |
| B5. Selection with corruption | R1 in At on A(3, 0.2), A(5, 0.1), A(10, 0.1), A(10, 0.2) [150] | C1, C1+F1; C3, C3+F1 on A(10, 0.1) only | 211 |
| B5b. Clustered selection | New cells: the selected sites as one contiguous foreground block, A(10, 0.05) (21 codons) and A(5, 0.1) (42 codons) [150] | C1, C1+F1, C3 | 122 |
| B6. Realigned indels | M4 [300]: realigned, the true alignment (diagnostic), F1, F3. New cell A(10, 0.1) with indels, realigned [150] | M4: C1 all four, C3 realigned; A(10, 0.1) with indels: C1, C1+F3 | 158 |
| B7. Genome-like mixture | 1,000 genes: 900 null and 100 selected (25 each of A(3, 0.2), A(5, 0.1), A(10, 0.05), A(10, 0.1)), each sequence corrupted independently at the empirical per-taxon rates (At 3.2%, Al 2.0%, Cr 5.6%, Es 8.3%, Br 15.4%), segment lengths and types from the P3A profile | C1, C3, each with and without F1, and F1fg (refits only where masking changes the alignment) | 272 |
| | | **Essential total** | **about 1,570** |
| Contingent | B2 with 3 starts, if D4 fails | | +240 |

C1∧C3 and MA(C1, C3) are computed for every cell that has both C1 and C3, with and without F1, at no extra cost.

### Optional, each needing its own go-ahead

| Block | Core-h | Why it might be worth it |
|---|---:|---|
| O1. A(3, 0.2) with indels, C1 and C1+F3 | 42 | F3's power cost at weak selection |
| O2. C3 and C3+F1 on A(3, 0.2) with R1 | 64 | The conservative candidate's weak-selection power under corruption |
| O3. N1 over the length and divergence grid, 1 start | 132 | Calibration with stronger purifying selection beyond L50 |
| O4. The rest of the power grid (ω+ 2, 3, 5, 10 at p+ 0.02; ω+ 2 at 0.05, 0.1; ω+ 5 at 0.05) | 354 | Complete power surface; low value, since power there is near the null rate |

## Decision criteria, fixed before any 2B job

- **D1. Calibration class** of a procedure in a null-like cell, at α = 0.05 with an exact 95% interval: *calibrated* if the upper bound is at most 5%; *unsafe* if the lower bound exceeds 5%; *inconclusive* otherwise. At 300 genes, *calibrated* needs at most 7 calls.
- **D2. Robustness:** D1 applied to M1 (δ 0.2, 0.3), M3, M4, R1-R3, and R1 with F1.
- **D3. A real difference between two procedures** on the same genes needs both an exact McNemar p below 0.05 and more discordant calls than the 95th percentile of the replicate discordance (B4's noise floor, for the matching scenario type). Power is reported by stratum (weak: ω+ ≤ 3; moderate: ω+ = 5; strong: ω+ = 10) with paired differences and their 95% intervals. A loss that meets D3 is *material*.
- **D4. One start in the clean null grid** is accepted only if, in every audited cell (B1, and M4 in B6), going from 1 to 3 starts changes at most 2 of 300 calls and leaves the D1 class unchanged, for both C1 and C3. Otherwise B2 runs with 3 starts. The outcome is recorded in the change log before B2 is launched.
- **D5. Triple hits:** C4 is a candidate only if, against C1 at 3 starts, it meets D3 in some cell or improves a D1 or D2 class. Otherwise, Double is enough in the settings tested.
- **D6. Masking, by non-inferiority.** For each selection cell, the paired difference in power (masked minus unmasked on the same genes; for corrupted cells, masked corrupted minus clean) with a Newcombe (1998, method 10) 95% interval. Margin: 5 percentage points absolute. *Non-inferior* if the lower bound is above -5 points; *inferior* if the upper bound is below -5 points; otherwise *inconclusive (insufficient precision)*, which is not evidence of non-inferiority. Reported by stratum (weak, moderate, strong, clustered). F1 (or F1fg) is recommended for evaluation as a package feature only if it moves corrupted nulls to *calibrated* and is non-inferior in clean, corrupted-with-selection, and clustered cells. If any of those is inferior, masking is recommended only as a sensitivity analysis, with both results reported; if some are inconclusive, that is stated as a limit.
- **D7. F3:** a candidate only if it makes M4 *calibrated* under C1 and is non-inferior (as D6) in A(10, 0.1) with indels, the test of whether it removes genuine selection along with alignment artifacts. The true-alignment control says whether M4's false positives come from alignment error (true alignment calibrated) or from indels as such.
- **D8. Conservative rule:** C1∧C3, MA(C1, C3), and C3 alone, each with and without F1, are compared by D1 and D2 in every null-like cell, by power per stratum, and by B7's false discovery proportion. A rule is recommended for conservative genome-wide scans only if it is *calibrated* in every null-like cell tested. Neither agreement nor MA becomes a default.
- **D9. Genome-like mixture:** after BH at q = 0.05 and 0.10, per procedure: the empirical false discovery proportion of this one mixture with an exact interval, power, and the false discoveries among corrupted genes. This is not an estimate of expected FDR in real genomes. With 80 selected genes it detects gross failures of FDR control only.

## Order

1. Commit the tooling: the new simulator cells (A2_p20, A3_p10, A5_p20, the four L25/L75 × divergence null cells, the two clustered cells, A(10, 0.1) with indels), F1fg and F3 in `phylo-windows.py`, the mixture builder, partial refits of changed alignments, and summaries with MA, noise floor, and consistency counts. Check that every existing cell's alignments are byte-identical to 2A's.
2. Smoke test in a separate root (4 genes per cell), then discard it.
3. Launch B1 and B3-B7.
4. Evaluate D4 on B1 and M4, record it in the change log, then launch B2 with 1 or 3 starts.
5. Report stage 2B, with the recommendation, before anything else is run.

## Reproducibility

- **The same dndsR and HyPhy as 2A:** a read-only library installed from the stage 2B tooling commit, whose dndsR package is identical to `5020a91`'s, and the image `sha256:e24b394df53597922308224fd2507956fb6f70bd146e470540b1c91179f1d17f`, both checked by every job (see the change log).
- **Independent replicates:** every 2B scenario has its own simulation seed (`phylo-stage2b-scenarios.tsv`, written by `phylo-stage2b-tables.py`), none shared with 2A. No 2A fit is reused. Within 2B, corrupted selection cells reuse their clean cell's seed, so the same genes are compared with and without corruption; and a filtered unit fits only the alignments its filter changed, taking the rest from the unit it was filtered from. The existing simulator cells reproduce 2A's alignments byte for byte.
- **Kept per gene:** calls with 1 and 3 starts, every start's likelihoods and AICc, the alignment checksum, the unit each fit came from, truth and corruption records (`genes.tsv.gz`), and a checksum of every simulated and filtered alignment (`checksums.md5`).
- One processor class (AMD EPYC 9654), resumable units with `DONE` markers, raw results outside the repository, summaries and the report in it.

## Limits

- **One phylogenetic configuration:** five taxa, one tree, one foreground clade (At+Al), κ = 2, P3A parameters. The 10-taxon simulation is postponed. Conclusions apply to the configurations tested and are not extended to other trees, taxon numbers, or foreground choices.
- The generating model is close to BUSTED's. Corruption is applied as one segment of empirical length per affected sequence. Indels are simulated at one rate.
- The genome-like mixture has one proportion of selected genes (10%) and no multinucleotide changes.
- Which real flagged segments are errors is unknown; only simulation knows the truth.

## Change log

| Change | Reason | Effect on results |
|---|---|---|
| Simulation seeds 8001-8150, none shared with 2A (7001-7135); no 2A fit is reused, so B1 and B2 no longer reuse 2A cells. To stay within the approved budget, selection cells (B3, B3b, B4 power cells, B5, B5b, B6's selection cell) have 120 genes instead of 150, as in 2A; the genome-like mixture has 800 genes (720 null, 80 selected, 20 per selected cell) instead of 1,000. Estimated: about 1,460 core-hours for the main blocks, plus B2 at 120 (one start) or 361 (three starts) | The maintainer required stage 2B's replicates to be independent of 2A's. Without reuse, the approved sizes would cost about 1,840, above the approved 1,566 | Precision: about ±9 points near 50% power at 120 genes. Non-inferiority depends mainly on discordant genes, which masking makes few |
| Fresh freeze at the tooling commit instead of 2A's library | No 2A fit is reused. The dndsR package is identical: `git diff 5020a91` over `R/`, `DESCRIPTION`, `NAMESPACE`, and `inst/` is empty. The image is the same by digest | None |
| B5 also fits foreground-only masking on the corrupted selection cells, and B6 fits all-sequence masking on A(10, 0.1) with indels, both only where the alignment changes | Masking-scope comparison wherever alignments differ; separates F1 from F3's gap step | A few core-hours |
| The maintainer approved 120 genes per selection cell and 800 in the mixture, and fixed the 5-point non-inferiority margin: results without the precision to establish it are reported as inconclusive, never by relaxing the margin | Clarification at approval | None |
| Budget enforced across the whole essential workflow, replacing the check before B2 alone. Every unit has an estimate (`phylo-stage2b-estimates.tsv`). Before a unit starts, its job adds the core-hours used (every unit's `meta.tsv`), the estimates of units running, and its own; if that exceeds the cap, the job stops, writes `BUDGET_STOP`, and the work pauses for the maintainer. The main blocks' cap is 1,686 (the ceiling of 1,806 less B2 with one start). B2 launches only if the core-hours used plus B2's estimate, scaled by the measured-to-estimated cost ratio (at least 1), stay within 1,806; otherwise it waits for approval. Simulation and summary jobs (minutes on 8 cores) are not counted | The maintainer's ceiling of about 1,806 core-hours, which the ±20% estimate does not extend | None |
| Smoke test (step 2) at `0aa137d`, root `p3s2b_smoke`: 101 of 101 units finished, 0 failed, the D4 audit and every summary table ran. It used 386 allocated core-hours, because each unit held 48 cores for 4 genes per cell; it is not part of the 1,806 ceiling, which covers the essential blocks. If the smoke counts against the 3,100 approved, about 2,079 remain, and 1,806 committed would leave about 270 uncommitted, against the reserve of at least 650 stated above. The main blocks are not launched until the maintainer decides how the smoke is counted. Later smoke tests should hold fewer cores | The study's budget | None on results; budget only |
