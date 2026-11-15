# Phylogenetic module: calibration stage 2B, and the stage 2 conclusions

Simulated data only. Design, decision criteria, and change log, all fixed before any 2B job: [phylo-calibration-stage2b.md](../dev/phylo-calibration-stage2b.md). Stage 2A: [phylo-calibration-stage2a.md](phylo-calibration-stage2a.md). Recommendations here are for the maintainer to decide on. **No default has changed, and nothing here is added to the package.**

Tables: `docs/validation/data/phylo_stage2b_*.tsv` (`units`, `rules`, `paired`, `starts`, `filtering`, `noise`, `mixture`, `audit`, `cost`). Unit names read `block__scenario__config__filter__seeds`. Per-gene results (`genes.tsv.gz`, with every start's likelihoods, AICc, and alignment checksums) and every HyPhy JSON are kept outside the repository, in the results root.

## Run

- 115 units, B1 and B3-B7 then B2: 115 units, all finished, **0 genes failed or excluded**. B2, the clean null grid, ran with three starts because the D4 audit failed in one cell (below).
- Frozen at `3603efc`, whose dndsR package is identical to 2A's `5020a91`, in a read-only library. Development image `sha256:e24b394df53597922308224fd2507956fb6f70bd146e470540b1c91179f1d17f` (HyPhy 2.5.103, not published), checked by every job. Every fit ran on AMD EPYC 9654.
- Inputs as 2A: `params.json` md5 `6a678fadefcbcf6b4d6d6377505997a9`, segment lengths md5 `30766db96c312d825e0c49e331db23bc`. Simulation seeds 8001-8150, none shared with 2A (7001-7135), so every 2B cell is an independent replicate. HyPhy seeds 1, 2, 3 (4, 5, 6 for the noise floor).
- Configurations: C0 BUSTED (dndsR's default), C1 `--multiple-hits Double`, C3 BUSTED-E with `Double`, C4 `Double+Triple`, all with `--srv Yes`. Filters: F1 masks dndsR's flagged dissimilar windows in every sequence (`mask`); F1fg only in foreground sequences (`mask_fg`); F3 is F1 plus codons within 3 of an internal gap (`gapflank`). Calls are at α = 0.05 on `lrt_best_fits`: the best alternative and best null likelihood over the starts, never the smallest p-value.
- Budget: 1,475.5 core-hours for stage 2B: 1,194.5 for B1 and B3-B7, 281.1 for B2 (estimated 361). The approved prospective ceiling was 1,806; the smoke test's 386 is recorded separately, as the maintainer ruled. Stage 2A used 635.

## Answers

### 1. False-positive control

Rate at α = 0.05, three starts, exact 95% interval; D1 class (*calibrated* if the upper bound is at most 5%). 300 genes per cell unless stated.

| Scenario | C1 | C3 | C1∧C3 | MA(C1, C3) |
|---|---|---|---|---|
| N0 L25 div ×0.5 | 0.3% (0-1.8) | 0% (0-1.2) | 0% | 0% |
| N0 L25 div ×2 | 0.7% (0.1-2.4) | 0% | 0% | 0.7% |
| N0 L50 div ×0.5 | 0.3% (0-1.8) | 0% | 0% | 0% |
| N0 L50 | 1.0% (0.2-2.9) | 0% | 0% | 0.7% |
| N0 L75 div ×0.5 | 0.3% (0-1.8) | 0% | 0% | 0.3% |
| N0 L75 | 1.3% (0.4-3.4) | 0% | 0% | 1.3% |
| N0 L75 div ×2 | 0% (0-1.2) | 0% | 0% | 0% |
| N0 L25 | 0.7% (0.1-2.4) | 0% (0-1.2) | 0% | 0.3% |
| N0 L50 div ×2 | 0.3% (0-1.8) | 0% | 0% | 0% |
| N1 L50 (stronger purifying selection) | 0% (0-1.2) | 0% | 0% | 0% |
| M1, multinucleotide, δ 0.2 | 2.0% (0.7-4.3) | 0% | 0% | 0.7% |
| M1, multinucleotide, δ 0.3 [150] | 0% (0-2.4); C4 0% | | | |
| M4, realigned indels | **8.3% (5.5-12.1), unsafe** | 0.3% (0-1.8) | 0.3% | 7.7% (4.9-11.3), inconclusive |
| M4, true alignment (diagnostic) | 2.0% (0.7-4.3) | | | |
| M4 + F1 | **8.0% (5.2-11.7), unsafe** | | | |
| M4 + F3 (gap flanks masked) | 1.7% (0.5-3.8) | | | |
| R1, one corrupted At segment, + F1 | 1.0% (0.2-2.9) | 0% | 0% | 0.7% |
| Mixture, 720 null genes, unmasked | 4.4% (3.1-6.2), inconclusive | 1.0% (0.4-2.0) | 1.0% | 4.0% (2.7-5.7), inconclusive |
| Mixture, 720 null genes, + F1 | 1.1% (0.5-2.2) | 0% (0-0.5) | 0% | 0.8% |

- **Correct nulls are calibrated under every configuration**: in the nine-cell grid of gene length (L25, L50, L75: the quartiles of codon counts in the real P3A orthogroups) by divergence (×0.5, 1, 2), plus N1, C1 ranges from 0% to 1.3% (largest upper bound 3.4%), and C3 makes no call in 2,700 genes. This replicates 2A (0-1.3%).
- **Multinucleotide changes are handled by double hits.** C1 is calibrated at δ 0.2 and 0.3, as in 2A (0.7%). Without double hits, 2A measured 97% (C0). Triple hits add nothing here (C4 0%).
- **Indels are the remaining failure of C1.** M4 gives 8.3%, replicating 2A's 8.7%. The true alignment brings it to 2.0%, so the excess comes from alignment error around indels, not from indels as such. F1 does not help: it flags only 33 of 300 M4 genes. F3 makes M4 calibrated (1.7%), but it costs power (below). C3 is calibrated on M4 (0.3%), as in 2A.
- **Alignment corruption.** One empirical-length corrupted segment in a foreground sequence gave 89-92% false positives under C1 and 24-28% under C3 in 2A, unmasked. Masking it is calibrated in 2B for C1 (1.0%) and C3 (0%), as in 2A.
- **In the genome-like mixture, every unmasked C1 false positive beyond the clean rate comes from foreground corruption.** 28 of the 720 null genes carry a corrupted foreground segment. C1 calls 24 of them (86%), against 1.4% of the 500 clean null genes and 0.5% of the 192 null genes corrupted only in the background. C3 calls 7 of the 28 (25%). Masking removes all of them under both C1 and C3.

### 2. Power

Three starts, 120 genes per cell, rate (exact 95% interval). Weak: ω+ ≤ 3; moderate: ω+ = 5; strong: ω+ = 10. A(ω+, p+): ω+ on a share p+ of foreground sites. Clustered: the same selected codons as one contiguous block.

| Stratum | Cell | C0 | C1 | C3 | C1∧C3 | MA |
|---|---|---|---|---|---|---|
| weak | A(2, 0.2) | | 2.5% (0.5-7.1) | 1.7% | 0.8% | 2.5% |
| weak | A(3, 0.05) | | 4.2% (1.4-9.5) | 0% | 0% | 4.2% |
| weak | A(3, 0.1) | | 8.3% (4.1-14.8) | 1.7% | 1.7% | 8.3% |
| weak | A(3, 0.2) | 10.8% (5.9-17.8) | 10.8% (5.9-17.8) | 5.8% | 5.8% | 10.8% |
| moderate | A(5, 0.1) | | 22.5% (15.4-31.0) | 8.3% | 6.7% | 23.3% |
| moderate | A(5, 0.2) | | 40.8% (32.0-50.2) | 26.7% | 25.8% | 40.8% |
| strong | A(10, 0.05) | 35.8% (27.3-45.1) | 32.5% (24.2-41.7) | 7.5% | 6.7% | 32.5% |
| strong | A(10, 0.1) | 65.8% (56.6-74.2) | 62.5% (53.2-71.2) | 40.8% | 39.2% | 62.5% |
| strong | A(10, 0.2) | 85.8% (78.3-91.5) | 77.5% (69.0-84.6) | 74.2% | 70.0% | 79.2% |
| clustered | AC(10, 0.05) | | 35.8% (27.3-45.1) | 12.5% | 12.5% | 35.0% |
| clustered | AC(5, 0.1) | | 33.3% (25.0-42.5) | 15.8% | 15.8% | 31.7% |
| indels | A(10, 0.1), realigned | | 69.2% (60.1-77.3) | | | |

Paired differences on the same genes (difference, Newcombe 95% interval, McNemar p, discordant genes). D3 calls a difference real only if p < 0.05 and the discordant genes exceed the noise threshold, 12 of 120 (below).

| Comparison | Cell | Difference | Discordant | p | D3 |
|---|---|---|---|---|---|
| C3 − C1 | A(5, 0.1) | −14.2 (−21.7, −7.2) | 21 | 0.0002 | material |
| C3 − C1 | A(5, 0.2) | −14.2 (−20.7, −7.5) | 19 | 0.0001 | material |
| C3 − C1 | A(10, 0.05) | −25.0 (−33.3, −16.9) | 32 | 2e-8 | material |
| C3 − C1 | A(10, 0.1) | −21.7 (−29.3, −13.4) | 30 | 9e-7 | material |
| C3 − C1 | A(10, 0.2) | −3.3 (−9.6, 2.8) | 14 | 0.42 | no |
| C3 − C1 | AC(10, 0.05) | −23.3 (−31.0, −15.8) | 28 | 7e-9 | material |
| C3 − C1 | AC(5, 0.1) | −17.5 (−24.5, −10.7) | 21 | 1e-6 | material |
| C3 − C1 | weak cells | −0.8 to −6.7 | 3-8 | 0.008-1 | no (below threshold) |
| C1 − C0 | A(3, 0.2), A(10, 0.05), A(10, 0.1) | 0 to −3.3 | 0-12 | ≥ 0.29 | no |
| C1 − C0 | A(10, 0.2) | −8.3 (−13.9, −3.4) | 10 | 0.002 | no (below threshold) |
| C4 − C1 | A(10, 0.1), A(10, 0.2) | −2.5, −3.3 | 11, 8 | 0.55, 0.29 | no |

- **Weak selection is barely detectable by any configuration**: 2.5-10.8% under C1, near the null rate under C3. That is a property of the test at this sample size (five taxa), not of a setting.
- **C3's protection costs material power at moderate and strong selection, and with clustered selection**: 14-25 points, each meeting D3. Only at A(10, 0.2) is the loss small. C1∧C3 is within 4.2 points of C3 everywhere.
- **Double hits cost little against BUSTED.** C1 − C0 is 0 to −3.3 points except at A(10, 0.2), −8.3 points, which is significant by McNemar but below the noise threshold. 2A measured a larger loss there with one start (91% to 70%); with three starts it is 86% to 78%.
- **Triple hits (C4) never meet D3 against C1** and change no D1 class (D5): Double is enough in the settings tested, as in 2A.
- **Model averaging gives BUSTED-E almost no weight** (mean Akaike weight 1.5-6%), so MA tracks C1 within 2 points in every cell. Three starts lowered the weight, rather than raising it: at A(10, 0.1), from 9.4% to 5.9%. Better BUSTED-E optimization does not make MA conservative here.

### 3. All-sequence versus foreground-only masking

D6: masked minus unmasked power on the same genes (for corrupted cells, masked corrupted minus clean), Newcombe 95% interval, 5-point non-inferiority margin.

| Selection | Cell | F1 | F1fg |
|---|---|---|---|
| clean, uniform (9 cells), C1 | weak to strong | 0 to +0.8, every interval within ±3.3: **non-inferior** | identical calls to F1 in all 9 |
| clean, uniform (9 cells), C3 | | 0 in every cell: **non-inferior** | |
| clean, indels | A(10, 0.1), realigned | +0.8 (−1.1, 2.9): non-inferior | |
| corrupted, C1 | A(3, 0.2) + R1 | +0.8 (−1.9, 3.9): non-inferior | same as F1 |
| corrupted, C1 | A(5, 0.1) + R1 | +1.7 (−3.6, 7.0): non-inferior | same |
| corrupted, C1 | A(10, 0.1) + R1 | +0.8 (−5.4, 7.1): inconclusive | same |
| corrupted, C1 | A(10, 0.2) + R1 | +3.3 (−0.8, 7.8): non-inferior | same |
| corrupted, C3 | A(10, 0.1) + R1 | +2.5 (−2.9, 7.8): non-inferior | |
| **clustered, C1** | **AC(10, 0.05)** | **−20.0 (−29.0, −10.7): inferior** | −20.0 (−29.2, −10.4): inferior |
| **clustered, C1** | **AC(5, 0.1)** | −9.2 (−17.4, −0.8): inconclusive | −11.7 (−20.0, −3.2): inconclusive |

- **On uniform selection, masking is free.** It changed 0 or 1 calls in every clean cell, replicating 2A (0 changes in A(10, 0.1) and A(10, 0.2)). It touches 0-8% of clean genes and removes under 0.05% of codons.
- **With corruption, masking restores clean power** within 0.8-3.3 points, replicating 2A. Unmasked, corruption raised C1's "power" to 92-99% by false calls: +81 points at A(3, 0.2).
- **Masking removes genuine clustered selection.** When the selected codons form one block, the dissimilar-window rule flags them as it would flag an error: it changed 119 of 120 genes and removed 1.8-2.0% of codons. AC(10, 0.05) loses 20 points, and AC(5, 0.1) 9 points. This is new in 2B; 2A had no clustered cell.
- **F1fg against F1.** Foreground-only masking removes fewer codons (mixture: 0.04% against 0.37%; AC(5, 0.1): 1.6% against 2.0%). It gave the same rates as F1 everywhere except AC(5, 0.1), where it lost 3 more genes (−2.5 points, p 0.25); in AC(10, 0.05) one call swapped each way. Neither avoids the clustered loss, because the selected block is in the foreground. **The evidence does not separate the two scopes.**
- **D6 outcome.** F1 moves corrupted nulls to *calibrated* and is non-inferior in clean and corrupted cells, but it is *inferior* in one clustered cell. By the prespecified rule, masking is recommended **only as a sensitivity analysis, with both results reported**, not as a package feature. One corrupted cell and one clustered cell are inconclusive: the precision is too low to call them either way.

### 4. One start against three, and optimization instability

| Cells | Calls changed, 1 vs 3 starts | Direction |
|---|---|---|
| Correct nulls (B1 and B2: 10 cells, C1 and C3) | 0-1 of 300 per cell | |
| M4 | C1: 3 (+1, −2); C3: 1; true alignment: 5 lost; F3: 6 (+1, −5) | mostly lost |
| Power, C1 | 0-12 of 120 (up to 10%) | mostly gained |
| Power, C3 | 0-18 of 120 (up to 15%) | mostly gained: A(10, 0.2) +18, 0 lost |
| Mixture, 800 genes | C1 6, C3 5 | gained |

Two independent three-start runs of C1 (seeds 1-3 against 4-6) on the same genes disagreed on 9 of 120 genes at A(10, 0.1) and 5 of 120 at A(10, 0.2). Pooled, that is 5.8%, with an upper 95% bound of 9.6%. This gives D3's threshold of 12 discordant genes in 120. The two runs' rates were within 1-4 points.

- **D4: one start fails the audit in one of 12 cells.** M4 under C1 changed 3 calls, against the limit of 2, so B2 ran with three starts, as prespecified. In correct nulls, extra starts changed at most one call in 300, as in 2A (0 of 600).
- **Under selection, three starts matter, most for BUSTED-E.** C3 at A(10, 0.2) went from 59.2% to 74.2%, and C3 at A(10, 0.1) from 35.8% to 40.8%. 2A's largest gain, C3 at A(10, 0.1) from 39% to 53%, replicates in direction but not in size.
- **Practical importance.** Rates in aggregate are stable within a few points once three starts are used. Individual calls near α are not: about 6% of genes in strong power cells flip between two stabilized runs. A per-gene call near the threshold should be read as uncertain, and a scan should report its seeds.
- HyPhy's reported LRT differed from twice its likelihood difference by more than 0.5 in at most 2 genes per unit. Inference uses the likelihoods.

### 5. False discoveries in the genome-like mixture (D9)

800 genes: 720 null, 80 selected (20 each of A(3, 0.2), A(5, 0.1), A(10, 0.05), A(10, 0.1)). Each sequence is corrupted independently at the empirical per-taxon rates. BH across the 800.

| Procedure | q | Discoveries | False | FDP (95% interval) | Power | False from corrupted foreground |
|---|---|---:|---:|---|---:|---:|
| C1 | 0.05 | 22 | 17 | **77% (55-92)** | 6.3% | 17 |
| C1 + F1 | 0.05 | 1 | 0 | 0% (0-97.5) | 1.3% | 0 |
| C1 + F1fg | 0.05 | 1 | 0 | 0% (0-97.5) | 1.3% | 0 |
| C3 | 0.05 | 7 | 5 | 71% (29-96) | 2.5% | 5 |
| C3 + F1 | 0.05 | 0 | 0 | | 0% | 0 |
| MA | 0.05 | 10 | 7 | 70% (35-93) | 3.8% | 7 |
| C1 | 0.10 | 25 | 19 | 76% (55-91) | 7.5% | 19 |
| C1 + F1 | 0.10 | 2 | 0 | 0% (0-84) | 2.5% | 0 |
| MA + F1 | 0.10 | 2 | 0 | 0% (0-84) | 2.5% | 0 |
| C3 + F1 | 0.10 | 0 | 0 | | 0% | 0 |

- **Unmasked, BH does not control false discoveries in this mixture.** C1 has 77% false discoveries at q = 0.05, and every one is a null gene with a corrupted foreground segment. BH assumes valid null p-values, and these are not. C3 and MA do not fix it (71%, 70%).
- **Masking removes every false discovery**, under C1, C3, and F1fg alike.
- **But very little survives BH.** Before correction, C1 called 29-30% of the 80 selected genes; after BH, 1-6%. With 10% of genes selected, mostly weakly, and a five-taxon test, the scan has almost no power after correction. This is one mixture: the masked FDP interval (0-97.5%) says nothing about FDR control, only that the unmasked failure is gross.

### 6. What replicates stage 2A, and what remains uncertain

| Finding | 2A | 2B (independent replicates) | |
|---|---|---|---|
| Correct nulls calibrated, all configurations | 0-1.3% | 0-1.3% in 10 cells, every upper bound ≤ 3.4% | replicates |
| Double hits fix multinucleotide changes | C1 0.7% | C1 2.0% (δ 0.2), 0% (δ 0.3) | replicates |
| M4 inflates C1; C3 does not | 8.7%; 0% | 8.3%; 0.3% | replicates |
| F1 does not fix M4 | 8.0% (C0) | 8.0% (C1) | replicates |
| Masking corrupted nulls calibrates them | 0-0.7% | 1.0% (C1), 0% (C3) | replicates |
| Masking costs nothing on clean uniform selection | 0 calls changed | 0-1 calls in 18 cell-configurations | replicates |
| Masking restores clean power under corruption | C0, C3 | within 0.8-3.3 points: C1 in 4 cells, C3 in 1 | replicates |
| C3 costs about half the power | 8-39% vs 34-61% (C1, one start) | material loss in 6 of 7 moderate, strong, and clustered cells | replicates |
| C4 equals C1 | no net difference | no D3 difference | replicates |
| MA tracks BUSTED, no corruption protection | one-start preview | three starts; Akaike weight of BUSTED-E falls | replicates, now with three starts |
| Three starts raise C3's power | +14 points at A(10, 0.1) | +5 at A(10, 0.1), +15 at A(10, 0.2) | direction replicates, size does not |
| Double hits' cost at A(10, 0.2) | 91% to 70%, one start | 86% to 78%, three starts | smaller, not material by D3 |
| C0 power at A(10, 0.2) | 91% (84-95) | 86% (78-92) | consistent within intervals |

New in 2B, not tested in 2A: masking removes clustered selection; F3 fixes M4 but costs power under selection with indels; masking prevents BH failure in a mixture; foreground-only masking behaves as all-sequence masking.

**Still uncertain:**
- Whether masking is non-inferior for moderately clustered selection (AC(5, 0.1)), and in one corrupted cell (A(10, 0.1) + R1). Both are inconclusive at 120 genes.
- F3's power cost under indels: −10.8 points (−17.2, −4.4). That meets D3 but is inconclusive by the D7 margin.
- C3 + F1 on realigned indels: not fitted. C3 alone is calibrated there.
- D3's noise threshold was measured only in power cells, and is applied to null cells as the prespecified rule says (29 of 300). Extra starts change at most one null call, so the threshold is probably conservative for nulls. That is why C3 − C1 on M4 (24 discordant, McNemar p 1e-7) does not count as a D3 difference.
- Whether any of this holds for other trees, more taxa, other foregrounds, or real data, where nobody knows which flagged windows are errors.

## Recommendations, where the evidence supports one

For the maintainer. None changes a default; masking is a validation script, not a dndsR feature.

1. **Use `--multiple-hits Double` (C1) rather than plain BUSTED (C0) whenever multinucleotide changes are possible.** C0 fails under them (97% in 2A). C1 is calibrated under every null tested except realigned indels, and its power cost against C0 is not material by D3. *Supported.*
2. **Use three starts combined by likelihood** (`--extra-seeds 2,3`, reading `lrt_best_fits`) for any inference under selection, and report the seeds. One start is enough for correct nulls, but not for M4 (D4) or under selection. *Supported.*
3. **Triple hits are not needed** in the settings tested (D5). *Supported.*
4. **Run masking as a sensitivity analysis, not a filter or a default.** Report results with and without it. A gene that is significant only unmasked, with a flagged foreground window, is a likely artifact. A gene that is significant only unmasked, with no corruption evidence but a block of divergent codons, may be clustered selection. The D6 rule gives this outcome: masking calibrates corrupted nulls but is inferior for clustered selection. *Supported.* Unmasked genome-wide scans are not trustworthy if the alignments may carry corrupted segments: 77% false discoveries in the mixture.
5. **For conservative genome-wide scans, C3 + F1 (or C1∧C3 + F1)** is calibrated in every null-like cell where it was run, in 2A and 2B. It costs 14-25 points of power at moderate and strong selection. *Supported, with the limit that C3 + F1 was not fitted on M4.*
6. **Model averaging is not a substitute for BUSTED-E or for masking** here: it tracks C1, and is not calibrated on M4 or in the unmasked mixture (D8). *Supported.*
7. **No recommendation** on all-sequence against foreground-only masking (the evidence does not separate them), or on F3 (not non-inferior, D7). Weak selection is out of reach of a five-taxon test under every configuration.

## Limits

- One phylogenetic configuration: five taxa, one tree, one foreground clade (At+Al), κ = 2, P3A parameters. Not extended to other trees, taxon numbers, or foregrounds.
- The generating model is close to BUSTED's. Corruption is one empirical-length segment per affected sequence, and indels are simulated at one rate.
- Selection cells have 120 genes: about ±9 points near 50% power. Several D6 comparisons are inconclusive for this reason, and inconclusive is not evidence of non-inferiority.
- The mixture is one realization, with 10% selected and no multinucleotide changes. Its FDP is not an expected FDR.
- Which real flagged windows are errors is unknown. Only simulation knows the truth.

## Cost

By block (core-hours): B1 207.9, B2 281.1, B3 281.3, B3b 30.0, B4 120.6, B5 154.1, B5b 89.8, B6 145.6, B7 165.0; total 1,475.5 (`phylo_stage2b_cost.tsv`). With 2A (635) and the smoke test (386), the study used 2,497 of the approved about 3,100.
