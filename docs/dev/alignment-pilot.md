# Alignment-method pilot (design)

Status: **bounded exploratory pilot, approved by the maintainer 2026-11-22 (`DECISIONS.md`, item 6); this design is fixed before any run.** It is not confirmatory. Its only possible outcomes are a report and, at most, a proposal for a confirmatory study on held-out data. No change to `phylo_prepare` follows from it alone.

## Why

In study B's BE3h cell (homologous codon indels, 8 taxa, no selection), aBSREL's raw branch p-values gave 7.1% false positives at 0.05 and its gene p-value 21% after `phylo_prepare` realigned the sequences, against 2.3% and 3.0% on the true alignment ([lineage-robustness-core.md](../validation/lineage-robustness-core.md)). So the excess comes from realignment. This pilot asks whether a practical, codon-aware alternative removes it without losing power, and whether a measure of alignment uncertainty, computed without the truth, marks the genes concerned.

## Questions

- **P1.** For each aligner, how close is the alignment to the truth, and what are aBSREL's and the branch models' false-positive rates on null genes?
- **P2.** What is the power on genes selected on the foreground, under the same indel process?
- **P3.** What does each aligner cost, and how many usable codon positions does it keep?
- **P4 (alignment uncertainty).** Does disagreement between two aligners, computed without the true alignment, mark the genes with alignment error or with realignment-induced false positives?

## Aligners (inputs as `phylo_prepare` receives them: unaligned CDS and the species tree)

| Arm | Method | Settings |
|---|---|---|
| A0 | `phylo_prepare` as it is | MAFFT 7.526 L-INS-i on the translations, back-translated (the default) |
| A1 | PRANK v.170427, codon mode | `prank -codon -once -t=<species topology>`; internal node labels removed (PRANK rejects them), no branch lengths, as dndsR has only the topology |
| A2 | MACSE 2.07 | `macse -prog alignSequences` with defaults; any frameshift character `!` becomes a gap and is counted |
| T | true alignment | diagnostic only, never used to align, choose settings, or tune |

The settings above are the tools' documented defaults for coding sequences, chosen before running. They are not changed after results. A1 and A2 alignments are written as `phylo_prepare` writes its own (same file names and FASTA headers), so every method downstream sees the same kind of input.

## Data

- **Null:** cell BE3h, T8b, seed 304, genes 1-150, regenerated in the pilot's root; their true alignments are checked against study B's BE3h unit (`study-B3`) byte for byte. Genes 151-300 of BE3h are **held out**: not aligned or examined here, kept for a confirmatory comparison.
- **Selected:** new cell BE3hS (BE3h with every gene selected on the foreground, site classes 0.05/0.3/10 there, as the mixture's selected genes), T8b, new seed 321, 150 genes. Added to `phylo-lineage-sim2.R` without changing any other cell.
- 300 codons, site classes 0.05/0.3/1 at 0.6/0.3/0.1 elsewhere, the stage 2 parameters (`params.json`), as study B.

## Methods run on every arm

aBSREL (defaults) and the branch models (M0, two-ratio, two-ratio null), as study B, through dndsR, on each arm's alignments. No masking (masking is a separate sensitivity analysis and did not change BE3h).

## Measures

- **Accuracy (against T):** codon-level sum-of-pairs score (share of the true homologous codon pairs recovered) and total-column score (share of true columns reproduced exactly), per gene, averaged.
- **Usable positions:** aligned length in codons; codon columns with no gap; columns with gaps in at most half of the taxa.
- **Size (null):** aBSREL `p_raw` per branch, `p_bh_lineage` per branch, gene `p_value`; `p_two_vs_m0` and `p_fg_omega_gt1` per gene; at 0.05 with Clopper-Pearson intervals; categories as study B.
- **Power (selected):** share of selected genes with aBSREL `p_raw` ≤ 0.05 on any foreground branch (t1, t2 or their stem), gene `p_value` ≤ 0.05, `p_two_vs_m0` ≤ 0.05, `p_fg_omega_gt1` ≤ 0.05, unadjusted (size and power on the same scale).
- **Cost:** wall seconds per gene for alignment, on one core.
- **P4 uncertainty:** per gene, the share of A0's homologous codon pairs that A1 does not have (A0-A1 disagreement), and the same for A0-A2. Flag = any disagreement (above 0), a threshold-free rule. Targets: (a) genes whose A0 alignment differs from T (sum-of-pairs below 1); (b) genes with a realignment-induced positive (a branch with `p_raw` ≤ 0.05 under A0 and above 0.05 under T). Sensitivity, and flags on genes with an exact A0 alignment, reported against study B's *detectable* rule (at least 80% of affected genes, at most 5% of unaffected). No other threshold is tried.

## Criteria (fixed now)

An alternative aligner is **nominated** for a confirmatory study only if all hold:

1. Its aBSREL `p_raw` false-positive rate on the null genes is lower than A0's, with the 95% gene-bootstrap interval of the difference (2,000 resamples of genes, seed 1) excluding 0, and its own upper bound is at most 7.5%.
2. Its power on the selected genes (aBSREL foreground `p_raw`) is non-inferior to A0's by the stage 2B rule: loss under 5 points and the interval's lower bound above −10.
3. Its alignment cost is at most 10 times A0's per gene.

Gene `p_value`, `p_bh_lineage` and the branch models are reported with the same intervals but do not decide nomination. If no aligner is nominated, current behaviour stays and the limitation is documented. A nominated aligner is proposed to the maintainer for confirmation on the held-out genes 151-300, new seeds and the 16-taxon tree, with its cost, before anything runs.

## Precision, cost, and budget

150 null genes give 1,950 branch tests (interval about ±1.2 points at 7%); 150 selected genes give an interval of about ±8 points on power. aBSREL and branch models about 64 core-seconds per gene per arm on 8 taxa (study B's pilot): 300 genes × 4 arms ≈ 21 core-hours; alignment ≈ 2 core-hours (MACSE about 14 s per gene, PRANK under 1 s with a guide tree). **Ceiling: 40 core-hours.** Intel nodes only (`--constraint=intel18`), because HyPhy's results differed between AMD and Intel CPUs in study B.

## Implementation

`tools/dev/validation/aln-pilot.sb` (one array task per data set) and `aln-pilot-summary.R`, from a frozen library at the commit that adds this design, root `aln-pilot` beside the study roots. Aligners from a separate conda environment (bioconda `prank=v.170427`, `macse=2.07`, `mafft=7.526`); they are validation tools, not dndsR dependencies.
