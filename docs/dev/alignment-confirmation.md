# Alignment confirmation: MAFFT (current) against PRANK codon mode (design)

Status: **approved by the maintainer 2026-11-23 (`DECISIONS.md`); this design and its criteria are fixed before any confirmatory result is seen.** It follows the exploratory [alignment pilot](../validation/alignment-pilot.md), whose BE3h genes 1-150 are the *development set*: they are used here only diagnostically (already reported), never pooled with the confirmatory data. The question is whether an optional PRANK alignment mode in `phylo_prepare` is justified. No default changes from this study; a positive result leads to a proposal for an optional mode, with its own tests and documentation, for the maintainer's approval.

## Arms

| Arm | Method | Settings |
|---|---|---|
| A0 | `phylo_prepare` as it is | MAFFT 7.526 L-INS-i on translations, back-translated (the default) |
| A1 | PRANK v.170427, codon mode | `prank -codon -once -seed=1 -t=<species topology without node labels>`: the pilot's settings plus a fixed seed (amendment below) |
| T | true alignment | diagnostic only: never used to choose settings |

MACSE is not carried forward: in the pilot it removed less than half of the excess at 25-29 times MAFFT's cost. No other aligner is added: the pilot gives no evidence that another inexpensive codon-aware method would do better, and a third arm would add cost without a reason fixed in advance.

## Data (all simulated, homologous indels as BE3h; held out from the pilot)

| Unit | Cell | Tree | Genes | Seed | Purpose |
|---|---|---|---:|---:|---|
| C-N8h | BE3h | T8b | 150 (genes 151-300) | 304 | null, the pilot's held-out genes |
| C-N8 | BE3h | T8b | 300 | 371 | null, new seed |
| C-S8 | BE3hS | T8b | 150 | 372 | every gene selected on the foreground (t1, t2, stem; classes 0.05/0.3/10) |
| C-M8 | BMIXh | T8b | 1,000 | 373 | 10% selected, mixture for false discoveries |
| C-N16 | BE3h | T16 | 300 | 374 | null |
| C-S16 | BE3hS | T16 | 150 | 375 | selected |
| C-M16 | BMIXh | T16 | 500 | 376 | mixture |

300 codons, site classes 0.05/0.3/1 at 0.6/0.3/0.1 elsewhere, `params.json` as study B. Runner `tools/dev/validation/aln-confirm.sb`, units `aln-confirm-units.tsv`, root `aln-confirm` frozen at the commit adding this design. aBSREL with defaults (every branch tested) on each arm. The branch models were near nominal under every aligner in study B and the pilot and are not rerun.

## Measures

- **Accuracy (against T):** codon sum-of-pairs and total-column scores per gene, as the pilot.
- **Completeness and validity:** share of genes with an alignment that reproduces the input CDS; aligned length, ungapped columns.
- **Type I error (null units):** aBSREL `p_raw` per branch and gene `p_value`, at 0.05, Clopper-Pearson intervals; gene-bootstrap difference A1 − A0 (2,000 resamples, seed 1).
- **False discoveries (mixtures):** false-discovery proportion with 95% gene-bootstrap interval (`lineage-mt-summary.R`'s method) for gene `p_adj`, `p_bh_lineage` per lineage family and pooled, and the number of lineages without selection that reject anything.
- **Power:** selected units: share of genes with `p_raw` ≤ 0.05 on any foreground branch, and gene `p_value` ≤ 0.05; mixtures: power of gene `p_adj` and `p_bh_lineage`. Differences A1 − A0 with gene-bootstrap intervals.
- **Cost:** per-gene wall seconds and peak resident memory of PRANK (`/usr/bin/time`, one core per gene); A0's whole `phylo_prepare` step; aBSREL's run time on the same genes for scale.
- **Reproducibility:** both aligners rerun on the first 20 genes of every unit; byte-identical output required.

## Criteria (fixed now)

A1 is **proposed as an optional `phylo_prepare` mode** only if all hold, at **both** 8 and 16 taxa:

1. **Branch Type I error.** On each null unit (C-N8h, C-N8, C-N16), A1's `p_raw` rate is lower than A0's, with the bootstrap difference interval excluding 0, and A1's own upper bound is at most 7.5%.
2. **Gene Type I error.** On the same units, A1's gene `p_value` rate is lower than A0's with the difference interval excluding 0, and A1 is not anticonservative (lower bound at most 5%).
3. **False discoveries.** In each mixture (C-M8, C-M16), A1 does not *fail* for gene `p_adj` or pooled `p_bh_lineage` (FDP lower bound at most 5%, the multiple-testing study's rule), and its FDP point estimate is below A0's.
4. **Power.** On C-S8 and C-S16, A1's foreground `p_raw` power is non-inferior to A0's (loss under 5 points, interval lower bound above −10). Gene `p_value` power is reported against T's, not A0's, because A0's is inflated by false positives.
5. **Completeness and reproducibility.** At least 99% of genes aligned and valid; all 20 reruns per unit byte-identical.
6. **Cost.** PRANK's per-gene time at most the per-gene time of the aBSREL run on the same genes (so the optional mode at most doubles an aBSREL analysis), at 16 taxa; peak memory at most 4 GB per gene. *This replaces the pilot's rule (at most 10 times MAFFT), which PRANK failed: MAFFT's cost is negligible beside the tests, and the maintainer asked whether the computational penalty is acceptable for the whole analysis. The new rule was set after the pilot showed PRANK at about 1.75 times aBSREL's cost on 8 taxa; it is stated here as a new criterion for a new study, not a reinterpretation of the pilot.*

If 1-5 hold but 6 fails, the report gives the cost and the maintainer decides. If 1 or 2 fail on either tree, current behaviour stays and the limitation remains documented. Whatever the outcome, external codon alignments stay supported, and the documentation says that an in-frame alignment is not evidence of correct homology.

**Diagnostic (not deciding):** T on every unit, to separate alignment error from other causes.

## Precision

Null: 150-300 genes per unit, 1,950-8,700 branch tests (branch-rate interval about ±1 point at 5%); gene rate about ±3 points. Selected: 150 genes, about ±4 points near 95%. Mixtures: about 100 and 50 selected genes; FDP intervals about ±5-10 points.

## Cost and budget

Measured in the pilot on 8 taxa: PRANK about 55, aBSREL about 28 core-seconds per gene. 16 taxa measured in a smoke test on development seeds before launch (recorded below). Projection: 8 taxa 1,600 genes × (PRANK 55 + aBSREL 3 × 28) ≈ 62 core-hours; 16 taxa 950 genes × (PRANK + 3 × aBSREL). **Ceiling 200 core-hours**; if the smoke test projects above, C-M16 is reduced to 300 genes first, and that is reported. Intel nodes only.

## Amendments before the confirmatory runs (2026-11-23, from a smoke test on development seeds 391-392, never pooled)

- **PRANK needs a fixed seed to be reproducible.** Without one, PRANK reran 3 of 10 smoke genes to a different alignment (its default seed varies between runs); MAFFT reran all 10 identically. With `-seed=1`, 4 of 4 16-taxon genes reran byte-identically. A1 therefore uses `-seed=1`, and any optional mode would set a seed and record it. The pilot ran without a seed; its results stand as one realization.
- **Measured cost.** PRANK on 16 taxa: about 115 core-seconds per gene (8 taxa: about 49), peak memory under 10 MB; `phylo_prepare` (MAFFT) peak about 260 MB for the whole step. Against aBSREL's measured 28 (8 taxa) and 64 (16 taxa) core-seconds per gene, PRANK is about 1.75-1.8 times aBSREL, so criterion 6 (at most 1 times) is projected to fail; as fixed above, the cost then goes to the maintainer for decision rather than ending the proposal. Projection for the whole study: about 60 core-hours on 8 taxa and 80 on 16, about 145 in all (ceiling 200).

## Maintainer's rulings of 2026-11-23, applied before any confirmatory result is read

**Cost criterion (6) unchanged, reported, not decisive.** Its outcome will be reported exactly as preregistered. Following the maintainer's ruling that runtime is not a scientific acceptance criterion (`DECISIONS.md`), the **decision rule for proposing an optional mode is criteria 1-5**; criterion 6 is reported beside it. This is a documented revision of the decision rule, made before the confirmatory results, not a reinterpretation of them.

**Practical utility assessment (separate from the gate; fixed now).** Reported from the confirmation's own measurements:

- per-gene and genome-scale cost: PRANK's core-seconds per gene on 8 and 16 taxa, projected to 10,000 orthogroups, beside the cost of the lineage workflow that follows (aBSREL, branch models, BUSTED C1, RELAX on the same trees, from the measured times of the studies), so the alignment's share of an end-to-end analysis is explicit;
- reuse: `phylo_prepare` writes alignments once and every method reads them, so the alignment cost is paid once per data set, not per model;
- reliability gained per unit of compute: false positives removed per 1,000 null genes (A0 − A1, `p_raw` and gene `p_value`), false discoveries removed in the mixtures, against PRANK's extra core-hours;
- power: any genuine selection signal lost (foreground `p_raw` and gene `p_value` in the selected units, gene `p_adj` and `p_bh_lineage` power in the mixtures), A1 against A0 and against T.

If A1 meets criteria 1-5, it is proposed as an **optional high-confidence alignment mode**, with this assessment; the default stays MAFFT until the maintainer decides and the mode has its own tests and documentation.

**Seed sensitivity (diagnostic, not part of the criteria).** On the pilot's development genes only (BE3h genes 1-150, already reported; not the held-out or confirmatory data): PRANK at seeds 1-5 and aBSREL on each (`prank-seeds.sb`, about 13 core-hours). Reported: share of genes whose alignment differs between seeds, sum-of-pairs against the truth per seed, the `p_raw` and gene `p_value` false-positive rates per seed, and how many gene-level calls change between seeds.
