# PRANK seed sensitivity (diagnostic)

Status: **diagnostic requested by the maintainer (2026-11-23), not part of the alignment confirmation's criteria.** Design: [alignment-confirmation.md](../dev/alignment-confirmation.md), "Seed sensitivity". Data: the alignment pilot's development genes only (BE3h, seed 304, genes 1-150, already reported); never the held-out or confirmatory data.

PRANK v.170427 codon mode (`-codon -once -seed=<s> -t=<species topology>`), seeds 1-5, then aBSREL (HyPhy 2.5.103, defaults) on each seed's alignments; runner `prank-seeds.sb`, about 13 core-hours.

| Measure | Result |
|---|---|
| Genes whose alignment differs between any two seeds | 48 of 150 (same lengths; gap placement differs) |
| aBSREL branch tests whose LRT changes | 18 of 1,950; largest change 0.137 |
| Branch `p_raw` crossing 0.05 between seeds | 0 (largest p-value range 0.004) |
| False-positive rate, `p_raw` | 3.59% at every seed |
| False-positive rate, gene `p_value` | 2.67% at every seed |
| Gene-level calls changing between seeds | 0 |

**Reading.** Without a fixed seed PRANK is not reproducible (smoke test: 3 of 10 reruns differed), so the seed must be fixed and recorded for exact reproducibility. On these genes, which seed is fixed does not change aBSREL's inference: the differences are confined to gap placement that moves few statistics, and by little. This is one indel process on 8 taxa; it does not establish seed-insensitivity elsewhere.
