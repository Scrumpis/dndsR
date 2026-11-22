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

**Pilot arms.** The pilot runs aBSREL and the branch models on every unit, and BUSTED only on the two error-free units, to measure its cost (the smoke test showed BUSTED with double hits and three starts at about 100-500 core-seconds per gene on 8 taxa, far above aBSREL's 20).

**Pilot gates (fixed now).** (1) every unit completes; (2) errors are generated as planned: one logged segment per alignment in E1 and E2; at least 95% of E3 orthogroups realigned and tested; in E4, t1's M0 branch length at least 3 times the median of the other terminal branches in at least 90% of genes; corruption rates in the mixture's log consistent with the planned rates (reported); (3) measured cost projects the full study at no more than 250 core-hours. If it projects above, arms are dropped in this order, keeping every aBSREL and branch-model arm: BUSTED masked runs on E2, E4, E5, E6, E7; then BUSTED unmasked on E5 and E7; and what was dropped is reported as not run.

## Pilot outcome (2026-11-21, development data; gates as fixed above)

`phylo-study-pilot-gates.R B` on the 12 pilot units (seeds 351-362), 9bea53a:

- **Passed:** all units complete; one segment per alignment in E1 and E2 (T8b, T16); 100% of E3 orthogroups realigned and tested (T8b, T16).
- **Failed, generator: E4.** The paralog's branch reached 3 times the median terminal branch in 75% of 20 genes (gate 90%). The split was too recent (stem 10% of t1's depth). Changed to half the depth before the main runs, and the E4 pilot is rerun (`phylo-study-B-pilot2-units.tsv`, the same seed) to check the gate again. No other cell's generator is affected.
- **Failed, cost:** 1,090 core-hours projected, against 250. Measured per gene (unmasked; masked about the same): aBSREL 28 core-seconds on 8 taxa and 89 on 16; branch models 36 and 228; BUSTED C1 226 and 536. The design had not projected branch models on 16 taxa, and BUSTED is 8-10 times aBSREL. After the drop order fixed above (BUSTED masked on E2, E4-E7, then unmasked on E5 and E7), it is still about 930; with no BUSTED at all, about 303.
- Mixture (20 genes): no gene selected (2 expected), corruption shares t3 0.05, t4 0.05, t6 0.15, t7 0.05: too few genes to compare with the planned rates; reported, not gated.

**Staged execution, pending the maintainer's decision on the rest.** Run now, within budget: every aBSREL and branch-model arm, unmasked and masked, on the 8-taxon and T16m units except E4 (`phylo-study-B-units-core.tsv`, the main seeds; about 144 core-hours). Held: E4 until its generator passes; the 16-taxon units (about 160 core-hours for aBSREL and branch models); every BUSTED arm (about 790 core-hours as designed, about 630 after the drop order). The held arms keep their seeds, so running them later reproduces the same simulated genes.

## Core arms outcome and an E3 generator defect (2026-11-22, after the core results)

The core arms and the E4 unit (re-piloted generator, gate passed 20 of 20) ran in full; report: [lineage-robustness-core.md](../validation/lineage-robustness-core.md). The direct-rerun rule above found no dndsR defect: E1 and E2 agree to 1e-12, and E3's few differences come from HyPhy behaving differently on AMD and Intel nodes.

**Generator defect found after the results, recorded rather than silently fixed.** E3's aBSREL inflation persisted on the true alignment (a diagnostic run, `p_raw` 11.2%, mostly terminal branches), which led to `add_indels()` (`phylo-sim.R`, shared with stage 2's M4). An insertion gave each descendant of the inserting branch its own random codons, so an insertion on an internal branch was a non-homologous segment in its descendants: an alignment error built into the "true" alignment, not the indel process the design names. The E3 results are kept and reported as measuring that mixed error. A homologous option (`add_indels(homologous = TRUE)`: one inserted sequence copied to every descendant) was added at `49fd39f` without changing the default, so stage 2 and E3 still reproduce. Cell BE3h reruns E3 with it, the same seed (304), the same methods, criteria and arms (`phylo-study-B-units-E3h.tsv`, root `study-B3`, about 11 core-hours, within the approved budget). No criterion or threshold changed.

**BE3h outcome.** Ran in full (aBSREL and branch models, as given and masked). Realigned: `p_raw` 7.1%, gene `p_value` 21%, `p_bh_lineage` 0.8%, `p_two_vs_m0` 5.3%. A diagnostic run of aBSREL on BE3h's true alignments (not a prespecified arm, the same as the E3 diagnostic) gave 2.3% and 3.0%, matching E0, so the remaining excess is alignment error from realignment. Direct HyPhy reruns of the 20 most significant genes agree 260 of 260. Reported in [lineage-robustness-core.md](../validation/lineage-robustness-core.md).

## Amendment before the remaining arms (2026-11-22, before any of them ran)

The maintainer approved the remaining arms with a ceiling of about 950 core-hours (all arms, the drop order not applied) and ruled that the unrun indel arms use the homologous-insertion generator. So **BE3h replaces BE3** for E3 on 16 taxa (seed 311, as planned) and for BUSTED on 8 taxa (seed 304, BE3h's seed). Every other unit keeps its cell, seed and settings. BUSTED runs as C1 (double hits, three starts), as given and masked, on E0, E1, E2, BE3h, E4, E5, E6, E7 and the mixture (8 taxa and T16m), and with aBSREL and the branch models on the 16-taxon units E0, E1 and BE3h. Code: the frozen library at `49fd39f` (root `study-B3`). Units: `phylo-study-B-units-rest.tsv`; BUSTED-only units are new rows named `-busted`, so they regenerate the same genes from the same seeds, which is checked against the earlier units' simulated alignments before results are read. Criteria, classification rules, mixture analysis and direct-rerun rule are unchanged.
