# Alignment-method pilot: MAFFT (current), PRANK codon mode, MACSE (results)

Status: **bounded exploratory pilot, run once to the design fixed before running ([design](../dev/alignment-pilot.md), commit `6962b4a`); not confirmatory.** Outcome under the preregistered criteria: **no aligner is nominated**. PRANK and MACSE both pass the false-positive and power criteria but fail the cost criterion (at most 10 times A0's alignment cost). Current `phylo_prepare` behaviour stays, and nothing in dndsR changes because of this pilot. What follows is for the maintainer to decide.

## Runs

| | |
|---|---|
| Frozen library | commit `6962b4a`, image `dndsr_v0.4.0_sha256-f1b42d9d.sif`, HyPhy 2.5.103, PAML 4.10.10 |
| Aligners | MAFFT 7.526 L-INS-i via `phylo_prepare` (A0); PRANK v.170427 `-codon -once -t=<species topology>` (A1); MACSE 2.07 `alignSequences` defaults (A2); true alignment (T, diagnostic only). Separate conda environment; validation tools, not dndsR dependencies |
| Null data | cell BE3h (homologous codon indels, 8 taxa, no selection), seed 304, genes 1-150. Their true alignments are byte-identical to study B's BE3h unit (150 of 150 compared). Genes 151-300 are held out, unexamined |
| Selected data | cell BE3hS (BE3h with every gene selected on t1, t2 and their stem, classes 0.05/0.3/10), seed 321, 150 genes |
| Methods | aBSREL (defaults, every branch tested) and the branch models (M0, two-ratio, two-ratio null), unmasked, on each arm |
| Nodes | Intel Xeon Gold 6148 only; Slurm array 21580170 (2 tasks) |
| Cost | 31.7 core-hours allocated (ceiling 40) |
| Summary | `tools/dev/validation/aln-pilot-summary.R`; tables in `data/alignment_pilot_*.tsv` |

All 600 alignments reproduced their CDS (0 invalid, 0 missing). MACSE inserted no frameshift characters.

## Accuracy and usable positions (P1, P3)

Mean per gene; sum-of-pairs (SP) and total-column (TC) scores at the codon level against the true alignment.

| Arm | Set | SP | TC | Length (codons) | Columns with no gap | Columns gapped in at most half the taxa |
|---|---|---|---|---|---|---|
| A0 MAFFT | null | 0.9971 | 0.9848 | 315.0 | 285.4 | 302.0 |
| A1 PRANK | null | 0.9987 | 0.9934 | 315.7 | 285.2 | 302.0 |
| A2 MACSE | null | 0.9976 | 0.9868 | 315.2 | 285.4 | 302.0 |
| A0 MAFFT | selected | 0.9961 | 0.9795 | 315.7 | 284.6 | 302.0 |
| A1 PRANK | selected | 0.9983 | 0.9912 | 316.5 | 284.5 | 302.0 |
| A2 MACSE | selected | 0.9971 | 0.9844 | 315.9 | 284.6 | 302.0 |

Every aligner is close to the truth on average, but A0 is inexact (SP below 1) in 138 of the 150 null genes. Usable positions differ by under one codon between aligners.

## Size on null genes (P1)

At 0.05, Clopper-Pearson 95% intervals, categories as study B. 1,950 branch tests and 150 genes per arm.

| Output | A0 MAFFT | A1 PRANK | A2 MACSE | T (true) |
|---|---|---|---|---|
| aBSREL `p_raw` per branch | 7.3% (6.2-8.6), anticonservative | **3.6% (2.8-4.5), calibrated** | 5.4% (4.4-6.5), within tolerance | 2.7% (2.0-3.5) |
| aBSREL gene `p_value` | 18.7% (12.8-25.8), anticonservative | 2.7% (0.7-6.7) | 9.3% (5.2-15.2), anticonservative | 2.0% (0.4-5.7) |
| aBSREL gene `p_adj` (rejections) | 13 of 150 | 0 | 3 | 0 |
| aBSREL `p_bh_lineage` (rejections) | 18 of 1,950 | 0 | 4 | 0 |
| `p_two_vs_m0` | 4.7% (1.9-9.4) | 4.0% (1.5-8.5) | 4.0% (1.5-8.5) | 3.3% (1.1-7.6) |
| `p_fg_omega_gt1` | 0% | 0% | 0% | 0% |

Every null gene is null, so any `p_adj` or `p_bh_lineage` rejection is false. The A0 rates match study B's full BE3h cell (7.1% and 21%), so this subset is representative.

Gene-bootstrap difference from A0 in `p_raw` false-positive rate (2,000 resamples, seed 1): A1 −3.7 points (−4.6 to −2.9); A2 −1.9 (−2.7 to −1.2); T −4.7 (−5.6 to −3.7). PRANK removes about four-fifths of the excess over the true alignment, MACSE about two-fifths.

## Power on selected genes (P2)

Share of 150 selected genes detected at 0.05, unadjusted.

| Output | A0 | A1 | A2 | T |
|---|---|---|---|---|
| aBSREL `p_raw` on any foreground branch | 99.3% | 99.3% | 99.3% | 99.3% |
| aBSREL gene `p_value` | 98.0% | 94.0% | 96.0% | 94.7% |
| `p_two_vs_m0` | 98.0% | 97.3% | 98.0% | 99.3% |
| `p_fg_omega_gt1` | 0% | 0% | 0% | 0% |

Foreground power differs by 0 points from A0 for both alternatives (interval 0 to 0: the same 149 genes). A0's higher gene `p_value` power is the same inflation that gives it 18.7% false positives; A1 matches the true alignment. `p_fg_omega_gt1` has no power here because only 10% of sites are selected, so the branch's mean ω stays below 1.

## Cost (P3)

| Arm | Core-seconds per gene (null; selected) | Ratio to A0 |
|---|---|---|
| A0 MAFFT via `phylo_prepare` | 0.64; 0.75 (wall time on 8 threads × 8, including start-up: an upper bound) | 1 |
| A1 PRANK | 48.9; 62.5 | 76-83 |
| A2 MACSE | 15.9; 21.9 | 25-29 |

**Deviation from the design.** The design expected PRANK to take under 1 s per gene; it took about 49-62 s. During the run I reported to the maintainer that PRANK was "still within the cost rule", taking A0 to be about 8 core-seconds per gene. That figure was wrong: A0 measured 0.64-0.75. Neither alternative is within 10 times A0. For scale, aBSREL and the branch models together cost about 64 core-seconds per gene on 8 taxa, so PRANK would roughly double the cost of an 8-taxon analysis. The cost criterion is fixed and is not changed after results.

## Nomination (preregistered criteria)

| Criterion | A1 PRANK | A2 MACSE |
|---|---|---|
| 1. `p_raw` false-positive rate lower than A0 (difference interval excludes 0) and own upper bound at most 7.5% | pass (−3.7, interval −4.6 to −2.9; upper 4.5%) | pass (−1.9, −2.7 to −1.2; upper 6.5%) |
| 2. Foreground power non-inferior (loss under 5 points, lower bound above −10) | pass (0, 0 to 0) | pass (0, 0 to 0) |
| 3. Alignment cost at most 10 times A0's | **fail** (76-83 times) | **fail** (25-29 times) |
| Nominated | no | no |

## Alignment uncertainty (P4)

Flag = any disagreement between A0's homologous codon pairs and the alternative's (threshold-free), judged by study B's detectability rule: at least 80% of affected genes flagged and at most 5% of unaffected genes.

| Disagreement | Target | Affected flagged | Unaffected flagged | Detectable |
|---|---|---|---|---|
| A0 vs A1 | A0 inexact (null) | 124/138 (90%) | 1/12 (8%) | no (unaffected over 5%; only 12 unaffected genes) |
| A0 vs A1 | A0 inexact (selected) | 134/139 (96%) | 3/11 (27%) | no |
| A0 vs A1 | realignment-induced positive (null) | 67/68 (99%) | 58/82 (71%) | no |
| A0 vs A2 | A0 inexact (null) | 125/138 (91%) | 5/12 (42%) | no |
| A0 vs A2 | A0 inexact (selected) | 126/139 (91%) | 6/11 (55%) | no |
| A0 vs A2 | realignment-induced positive (null) | 64/68 (94%) | 66/82 (80%) | no |

Disagreement between aligners finds almost every gene with an alignment error, but nearly every gene with indels has one, so the flag does not separate genes with induced false positives from others (71-80% of genes without one are flagged). As a per-gene filter it would be close to excluding every gene with indels. It does not meet the rule and is not proposed as a diagnostic.

## Conclusions

- Realignment error in these simulations comes from MAFFT on translated proteins. A codon-aware, phylogeny-aware aligner (PRANK) reduces aBSREL's branch false-positive rate from 7.3% to 3.6% and its gene rate from 18.7% to 2.7%, near the true alignment's, with no loss of foreground power. MACSE helps less.
- By the preregistered rule no aligner is nominated, because both cost far more than 10 times the current step. This rule was set before results and stands; whether the cost criterion was the right one is for the maintainer to judge, and any change would be a new, labelled decision, not a reinterpretation of this pilot.
- The branch models were near nominal under every aligner (3.3-4.7%), consistent with study B.
- Aligner disagreement is not a useful per-gene uncertainty flag here.
- Current behaviour stays. The limitation is documented in [lineage-selection.md](../lineage-selection.md): aBSREL is anticonservative on alignments realigned with indels, and the dissimilar-window flag does not catch it.

## Limits

One indel process (homologous codon indels from study B's simulator), 8 taxa, 300 codons, one selection pattern, 150 genes per set; the simulator's indels may favour or disfavour particular aligners. PRANK was given the topology without branch lengths, as dndsR has only the topology. Held-out genes 151-300 are untouched and remain available.

## Options for the maintainer (recommendations, not decisions)

1. **Keep current behaviour** and the documented limitation (the default under the preregistered outcome).
2. **Offer PRANK codon mode as an optional aligner in `phylo_prepare`**, off by default, documented as reducing aBSREL's false positives in this pilot at about 50-60 s per gene. This would be a new interface and a new optional dependency, and should first be confirmed on held-out genes 151-300, new seeds and the 16-taxon tree, with its cost stated (about 8 core-hours for 300 genes × 2 aligners on 8 taxa, more on 16).
3. **Recommend in the documentation** that users who care about aBSREL results on indel-rich data align with a codon-aware aligner themselves and pass the alignments in, without changing dndsR, citing this pilot as exploratory.

My recommendation is option 3 now, with option 2 considered only after the confirmatory run. Option 3 changes no software. I have not implemented any of these.
