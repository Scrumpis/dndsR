# BUSTED-E as an alternative to masking: pilot results

Status: **pilot run once to the design fixed before running ([busted-e-pilot.md](../dev/busted-e-pilot.md), commit `dbe95e5`).** Outcome: **BUSTED-E does not improve on masking; every criterion fails.** Mechanism and integration checks: [busted-e-forensics.md](busted-e-forensics.md). No further BUSTED-E work is planned (maintainer, 2026-11-23).

## Runs

BUSTED-E (`phylo_busted(error_sink = "Yes")`, otherwise study B's C1 settings: double hits, three starting points, seed 1, HyPhy 2.5.103) on study B's exact inputs (8 taxa): E-E0 (300 error-free genes), E-E1 (the first 300 of study B's 600 genes with a corrupted segment in a foreground sequence), E-MIX (study B's 1,000-gene mixture, 116 selected). BUSTED and masked BUSTED results are study B's, on the same genes. 92.7 core-hours allocated (ceiling 150). Tables: `data/busted_e_pilot_*.tsv`.

## Results against the criteria

| Criterion | BUSTED | Masked BUSTED | BUSTED-E | Outcome |
|---|---|---|---|---|
| 1. Gate, E-E1 false positives (upper bound ≤ 7.5%) | 95.0% (91.9-97.2) | 0.7% (0.1-2.4) | **32.3% (27.1-37.9)** | fail |
| 1. Gate, E-E0 not anticonservative | 0.3% | 0.3% | 0% (0-1.2) | pass |
| 2. E-MIX false-discovery proportion (upper bound ≤ 7.5%) | 51.6% (62 discoveries) | 0% (16) | **55.6% (21.2-86.3), 9 discoveries** | fail |
| 3. E-MIX power above masked BUSTED (interval excluding 0) | 25.9% | 13.8% | **3.4%**; difference from masked −10.3 points (−16.4 to −5.1) | fail |

Gate 1 failed, so, as the design specified, the mixture was still run (within the ceiling) and is reported as descriptive. In the mixture BUSTED-E both made a majority of false discoveries and found fewer selected genes than masked BUSTED. The power figures after correction carry the BH-adaptivity caveat documented in the forensics (finding 5); the raw-scale comparison does not change the conclusion, since BUSTED-E's discoveries are few and mostly false.

## Reading

Under substantial corruption of one foreground sequence by a non-homologous segment, BUSTED-E is not a remedy: it neither controls false discoveries nor preserves power relative to masking. The forensics show why (the segment fits as moderate ω, which the error class, bounded at ω ≥ 100, cannot represent) and that dndsR reproduces HyPhy exactly. This is a limitation under this specific corruption, not a general verdict on BUSTED-E. No confirmation study is proposed.
