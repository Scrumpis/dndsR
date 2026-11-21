# Study B: lineage tests under real-data errors (design)

Status: **approved by the maintainer 2026-11-21 (pilot first, then the full study if the pilot passes its gates, within about 250 core-hours); amended below before the pilot ran.** Gap 2 of the [audit](validation-audit.md): BUSTED's dominant real-data failure (a corrupted foreground segment: 86% of such null genes called, 77% false discoveries after BH in a genome-like mixture, [stage 2B](../validation/phylo-calibration-stage2b.md)) was never measured for aBSREL or the branch models, and the lineage studies simulated perfect alignments on trees of 5-16 taxa. BUSTED itself was tested only on 5 taxa.

## Questions

- **B1.** Under each error below, what is the false-positive rate of aBSREL (`p_raw` per branch, `p_bh_lineage`, gene `p_adj`), `p_two_vs_m0`, `p_fg_omega_gt1`, and BUSTED (C1: double hits, three starts), with no selection anywhere?
- **B2.** In a genome-like mixture (null and selected genes, errors at empirical rates), what share of discoveries after the documented correction is false (`p_bh_lineage` for aBSREL; BH across genes for the gene-level tests)?
- **B3.** Does masking flagged windows (`phylo-windows.py`, the stage 2 sensitivity analysis) restore calibration, and at what cost in power, for these methods and trees? Masking stays a sensitivity analysis unless this study and the maintainer decide otherwise; no automatic exclusion is introduced.
- **B4.** Which existing QC flags (`dissimilar_window`, `high_ds`, `few_synonymous_changes`, dS above 3 times the median) mark the false positives, and how many true positives they would remove if used as filters (reported, not adopted).

## Errors (each a null cell; foreground = the first two taxa and their stem, as in the lineage studies)

| Cell | Error | Source of the generator |
|---|---|---|
| E0 | none (reference, same seeds' structure) | `phylo-lineage-sim.R` |
| E1 | one corrupted segment of empirical length in one foreground sequence | stage 2 corruption |
| E2 | one corrupted segment in one background sequence | stage 2 corruption |
| E3 | codon indels at rate 0.02, realigned by `phylo_prepare` (true alignment discarded) | stage 2 M4 |
| E4 | a misassigned ortholog: one foreground taxon's sequence replaced by a paralog that diverged at the root of the tree (dS of that tip several times the median) | new, small addition to the simulator |
| E5 | gene-tree discordance (the foreground taxa's gene tree differs from the species tree) | stage 2 M5, extended to T8b |
| E6 | 4 of 16 taxa missing per gene | `T16m` |

Trees: T8b for every cell; T16 for E0, E1, E3 (the errors that inflated BUSTED). Site classes 0.05/0.3/1 at 0.6/0.3/0.1 (B0 of the lineage study), 300 codons.

## Criteria (fixed before running)

- Size at 0.05 judged by the lineage studies' categories and Bradley's tolerance (2.5-7.5% for the whole interval). A method is **robust to an error** if E-cell size is within tolerance or calibrated; **not robust** if outside or anticonservative; otherwise inconclusive.
- B2: false discovery proportion with its 95% interval, compared with the nominal q.
- B3: masking is **non-inferior** in power if the loss is under 5 points with the interval's lower bound above -10 (the stage 2B D6 rule).

## Sizes and precision

300 genes per null cell (an interval of about plus or minus 2.5 points at 5%); 600 for E1 on T8b, where the decisive comparison is. Mixture: 1,000 genes (900 null, 100 selected on the foreground with omega 10 on 10% of sites), errors at the empirical per-taxon rates of stage 2.

## Seeds and cost

New seeds 301-312. aBSREL about 20 core-seconds per gene on 8 taxa and 64 on 16; branch models (M0, two-ratio, null) about 36 on 8 taxa; BUSTED C1 with three starts measured in the pilot. About 4,500 gene-cells, each run unmasked and masked: **about 250 core-hours**. Pilot: 20 genes per cell (about 10 core-hours), checking that each error is generated as intended (for example that E4's paralog has the planned dS) before the full runs.

## Amendments before the pilot (2026-11-21, before any run of this study)

The maintainer asked for defensible input-quality safeguards: evaluate incorrect orthology, alignment errors, indels, missing taxa, discordance, and recombination; decide which problems can be detected automatically, which warrant warnings, and which require exclusion; and separate upstream model limitations from dndsR defects. Nothing below was chosen from results of this study.

**Implementation.** `phylo-lineage-sim2.R` (cells `BE0`, `BE3`, `BE4`, `BE5`, `BE7`, `BMIX`), `phylo-lineage-study2.sb` (corruption by `phylo-windows.py`, realignment by `phylo_prepare`, masked reruns), units `phylo-study-B-units.tsv` (12 units, seeds 301-312) and `phylo-study-B-pilot-units.tsv` (20 genes each, seeds 351-362, development only).

**Added cell E7, recombination** (seed 308): the second half of each gene's codons evolves on E5's discordant gene tree (t2 and t3 exchanged), the first half on the species tree; one breakpoint at the middle.

**Generators, made precise.** E4: a paralog of t1 that split from the species lineage at the root (a new root, stem 10% of t1's depth; the paralog's branch is t1's depth plus that stem), replacing t1. E5: t2 and t3 exchanged, which is one nearest-neighbour interchange at the stem of the foreground clade (t1, t2), so the foreground clade is not monophyletic in the gene tree. E3: indels as stage 2's M4 (rate 0.02, mean 2 codons), realigned by `phylo_prepare` from unaligned CDS; the true alignment is kept for diagnosis. E1 and E2: `phylo-windows.py corrupt internal` on t1 or t8. Mixture: each taxon's corruption rate is one of stage 2B's five empirical per-taxon rates (t1 0.032, t2 0.020, t3 0.056, t4 0.083, t5 0.154, t6 0.056, t7 0.083, t8 0.032); 10% of genes selected on the foreground with classes 0.05/0.3/10 there.

**Methods and settings.** aBSREL (defaults), branch models (M0, two-ratio, two-ratio null), BUSTED as stage 2B's C1 (double hits, three starting points), each on the alignment as given and on its `mask` copy.

**Classification of each error, per method (fixed now).**

- *Robust*: size within tolerance (criteria above) or calibrated.
- *Detectable*: a signal computable from the input marks at least 80% of affected genes while marking at most 5% of the error-free genes of E0. Candidates, all existing: the window profile (`dissimilar_window`), `high_ds`, `few_synonymous_changes`, a tip's dS more than 3 times the gene's median tip dS, and, for E5 and E7, `gene_tree_mode = "report"`'s topology comparison where a gene tree is available. Sensitivity and specificity are reported for every pair.
- *Warn*: not robust, detectable, and masking does not restore calibration, or costs more power than the non-inferiority rule allows: a warning naming the genes is justified.
- *Exclude*: not robust, detectable, and size returns within tolerance when the flagged genes are removed, at a power cost within the non-inferiority rule. Exclusion is recommended to the maintainer, never applied by default without approval.
- *Limitation*: not robust and not detectable: documented as a limit of the method.

**Upstream limitation or dndsR defect (fixed now).** For each method and error cell where size is outside tolerance, the 20 genes with the smallest p-values are rerun directly in HyPhy or codeml, without dndsR, from the alignment and tree dndsR wrote (as `phylo-lineage-relax-direct.sb` and `phylo-direct-lineage.sh` did). If the direct statistics agree with dndsR's (LRT within 1e-3 or a relative 1e-4), the failure belongs to the model or engine on that input; if not, it is a dndsR defect and is fixed before anything else.

**Pilot gates (fixed now).** (1) every unit completes; (2) errors are generated as planned: one logged segment per alignment in E1 and E2; at least 95% of E3 orthogroups realigned and tested; in E4, t1's M0 branch length at least 3 times the median of the other terminal branches in at least 90% of genes; corruption rates in the mixture's log consistent with the planned rates (reported); (3) measured cost projects the full study at no more than 250 core-hours. If it projects above, arms are dropped in this order, keeping every aBSREL and branch-model arm: BUSTED masked runs on E2, E4, E5, E6, E7; then BUSTED unmasked on E5 and E7; and what was dropped is reported as not run.
