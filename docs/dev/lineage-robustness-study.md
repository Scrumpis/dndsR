# Study B: lineage tests under real-data errors (design)

Status: **draft for the maintainer's review; not run.** Gap 2 of the [audit](validation-audit.md): BUSTED's dominant real-data failure (a corrupted foreground segment: 86% of such null genes called, 77% false discoveries after BH in a genome-like mixture, [stage 2B](../validation/phylo-calibration-stage2b.md)) was never measured for aBSREL or the branch models, and the lineage studies simulated perfect alignments on trees of 5-16 taxa. BUSTED itself was tested only on 5 taxa.

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
