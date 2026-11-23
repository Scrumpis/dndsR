# BUSTED-E as an alternative to masking (pilot design)

Status: **bounded pilot, approved by the maintainer 2026-11-23 (`DECISIONS.md`); design fixed before any run.** Not authorization for a new stable method or a default change.

## Why

In study B's mixture (8 taxa, 1,000 genes, foreground corruption at stage 2B's empirical rates), BUSTED C1 had 51.6% false discoveries after BH; masking the flagged windows removed them all but cut power from 25.9% to 13.8% (−12.1 points, −20.0 to −5.3; [report](../validation/lineage-robustness-core.md#remaining-arms-16-taxa-and-busted-2026-11-23)). BUSTED-E (`--error-sink Yes`) adds an error-sink class meant to absorb residual alignment error. The question is whether it controls false discoveries in that mixture with less power loss than masking.

**Prior evidence (5 taxa, stage 2B, same HyPhy 2.5.103, so compatibility is established):** unmasked BUSTED-E gave 24-28% false positives with a corrupted foreground segment and 71% false discoveries in a genome-like mixture; it cost 14-25 points of power on clean data. The pilot tests whether this holds on 8 taxa, where the error sink may be better identified; the prior evidence makes failure likely, so the pilot stops early if its first gate fails.

## Arms (all study B's inputs, byte for byte)

| Arm | Source |
|---|---|
| BUSTED C1 (`--multiple-hits Double`, `--starting-points 3`, seed 1) | study B's results (`study-B3`, frozen `49fd39f`), not rerun |
| BUSTED C1, masked | study B's results |
| BUSTED-E: C1 with `--error-sink Yes` | run here, `phylo_busted(error_sink = "Yes")`, on the same alignments |

BUSTED-E is tested unmasked: it is offered as an alternative to masking, not as an addition to it (stage 2B already found masked BUSTED-E calibrated but less powerful than masked BUSTED).

## Data

| Unit | Source unit (study B) | Input | Genes |
|---|---|---|---:|
| E-E0 | B-E0-T8b-busted | true alignments | 300 |
| E-E1 | B-E1-T8b-busted | corrupted (foreground segment in t1) | 300 (the first 300 of 600, by name) |
| E-MIX | B-MIX-T8b-busted | corrupted mixture | 1,000 (116 selected) |

Runner `tools/dev/validation/busted-e-pilot.sb`, units `busted-e-pilot-units.tsv`, root `busted-e` frozen at the commit adding this design. Comparisons with study B's arms are on the same genes (E-E1: the same 300).

## Measures

False-positive rate at 0.05 on `p_value` (E-E0, E-E1), Clopper-Pearson; in E-MIX, BH-corrected `p_adj`: discoveries, false-discovery proportion (Clopper-Pearson, as study B), power over the 116 selected genes, and the power difference from masked BUSTED with a 95% gene-bootstrap interval (2,000 resamples, seed 20261122, as `phylo-study-B-exclude.R`); run time per gene against BUSTED's.

## Criteria (fixed now)

BUSTED-E **improves on masking** only if all hold:

1. **Gate 1 (E-E1, run first with E-E0):** false-positive rate with a corrupted foreground segment within tolerance (upper bound at most 7.5%). E-E0 not anticonservative (lower bound at most 5%).
2. **FDR (E-MIX):** FDP upper bound at most 7.5% (the multiple-testing study's "controls").
3. **Power (E-MIX):** power above masked BUSTED's (13.8%), with the bootstrap difference interval excluding 0.

If gate 1 fails, E-MIX is still run only if its cost is within the ceiling, because the maintainer asked specifically about the mixture's power and FDR; it is then reported as descriptive. If all three hold, a confirmation on an independent mixture (BMIX, new seed 381, 1,000 genes; BUSTED, masked BUSTED and BUSTED-E; about 3 × 63 core-hours) is proposed with its cost, not run without the maintainer's go-ahead, since it would exceed this pilot's budget.

## Cost

BUSTED C1 on 8 taxa measured at about 226 core-seconds per gene; BUSTED-E measured in a smoke test (8 mixture genes) before launch and recorded below. **Ceiling 150 core-hours.**

**Smoke test before launch (8 mixture genes, development only, not pooled):** BUSTED-E ran on HyPhy 2.5.103 through `phylo_busted(error_sink = "Yes")` without error (dndsR prints that BUSTED-E is experimental); 8 genes in 336 s on 8 cores, one gene per core, so at most about 336 core-seconds per gene (BUSTED C1: 226). Projection for 1,600 genes: about 110-150 core-hours, within the ceiling.
