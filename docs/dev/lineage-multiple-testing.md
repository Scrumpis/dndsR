# Multiple testing in the lineage outputs: definitions and a mixture study (design)

Status: **the audit of definitions is done (below). The mixture study answers the maintainer's item 3 of 2026-11-22 (`DECISIONS.md`); its design is fixed before any run.** Results: [validation/lineage-multiple-testing.md](../validation/lineage-multiple-testing.md). It is a follow-up to study B, labelled as such. Its preregistered criteria are not changed after results.

## Why

The study B core report called `p_bh_lineage` calibrated in BE3h because 0.8% of its branch tests were below 0.05. For a BH-adjusted output that rate guarantees nothing. In all-null data, BH's false-discovery rate is the probability that a family makes any rejection, and `p_bh_lineage` made false discoveries in 11 of 13 lineage families. aBSREL's gene `p_adj` made 30 in 300 null genes ([correction](../validation/lineage-robustness-core.md)). No study has measured these outputs' false-discovery control or power with null and selected genes mixed. The earlier power figure for `p_bh_lineage` (93%, [lineage calibration](../validation/phylo-lineage-calibration.md)) came from cells where every gene was selected, which favours BH.

## Audit of the definitions (from the code: `R/phylo_hyphy.R`, HyPhy 2.5.103 `aBSREL.bf` and `libv3/convenience/math.bf`)

| Output | Computed by | Hypothesis | Family | Error controlled |
|---|---|---|---|---|
| `p_raw` (per branch and gene) | HyPhy: LRT of the branch's adaptive model against the same with ω ≤ 1 for that branch; p from the mixture 0.5 χ²₀ + 0.2 χ²₁ + 0.3 χ²₂ | this branch, in this gene, has some sites with ω > 1 | none | one test's type I error |
| `p_holm` (per branch and gene) | HyPhy `math.HolmBonferroniCorrection`: the k-th smallest of m tested p-values times (m − k + 1), capped at 1, **without** the running maximum of standard Holm | as `p_raw` | the m branches tested in that gene (all branches by default) | rejecting each branch with `p_holm` ≤ α lies between Holm and Hochberg's step-up procedure; the family-wise error within the gene is controlled when the branches' p-values are independent or positively dependent (Hochberg's condition), not under arbitrary dependence |
| gene `p_value` | dndsR: the smallest `p_holm` | some tested branch of this gene has sites with ω > 1 (global null: none has) | one gene | as a global test it is Hochberg's (the minimum over k of (m − k + 1) p₍ₖ₎), not Bonferroni's m·p₍₁₎ as the documentation implies; they differed in 1-4 genes per unit of study B and never across 0.05 |
| gene `p_adj` | dndsR: BH on `p_value` across tested orthogroups | as `p_value` | all tested genes | false-discovery rate over genes, if each `p_value` is valid |
| `p_bh_lineage` (per branch and gene) | dndsR: BH on `p_raw` across the tested genes in which that species-tree branch was tested | as `p_raw` | one family per species-tree branch (13 on an 8-taxon tree); families are separate | false-discovery rate **within one lineage**, if each `p_raw` is valid. Nothing is controlled across lineages: pooling discoveries from all lineages is a separate question this design measures |
| branch models `p_adj_*` | dndsR: BH across orthogroups of each test's p | per test (two-ratio vs M0, foreground ω > 1, free vs M0) | all tested genes, per test | false-discovery rate over genes, per test |

Consequences for documentation (to be written after the study, and reviewed): `p_raw` is a per-branch test with no correction; `p_holm` and gene `p_value` answer "is any branch of this gene selected?" within one gene; `p_adj` is the gene-level output when many genes are tested; `p_bh_lineage` answers "in which genes is this particular lineage selected?" for a lineage chosen in advance. None of these is valid when `p_raw` itself is not, as under realignment error (BE3h).

## Mixture study

**Questions.**

- **M1.** With null and selected genes mixed, do aBSREL's `p_bh_lineage` (per lineage family and pooled over families), gene `p_adj`, `p_holm` per branch, and the branch models' `p_adj_two_vs_m0` and `p_adj_fg_omega_gt1` keep the false-discovery proportion at or below 0.05, with error-free alignments and with realigned indels?
- **M2.** What power does each have on the selected genes and branches?
- **M3 (exploratory, computed afterwards from the same p-values, no code change).** One BH family over every branch test of every gene, as an alternative for exploratory all-branch analyses.

**Cells.** 8 taxa (T8b), 1,000 genes each, 300 codons, site classes 0.05/0.3/1 at 0.6/0.3/0.1. Ten percent of genes are selected on the foreground (t1, t2 and their stem, classes 0.05/0.3/10 there), as study B's mixture.

| Unit | Cell | Errors | Seed |
|---|---|---|---|
| Q-MIX-T8b | BMIX, no corruption | none | 331 |
| Q-MIXh-T8b | BMIXh (new: BMIX's selection with BE3h's homologous indels), realigned by `phylo_prepare` | indels | 332 |
| Q-MIXh-T8b, true alignment (diagnostic arm) | same genes, aBSREL on the true alignments | none | 332 |

**Methods.** aBSREL with defaults (every branch tested) and the branch models (M0, two-ratio, two-ratio null) through dndsR, unmasked; the runner and settings of study B (`phylo-lineage-study2.sb`).

**Truth.** A gene is selected or not (`genes.tsv`). A branch test is true when its branch is a foreground branch of a selected gene; every other branch test is null. In selected genes the background branches are null.

**Measures.** For each output at 0.05: discoveries, false-discovery proportion with a 95% gene-bootstrap interval (2,000 resamples of genes, seed 1), and power (share of true genes, or true foreground branch tests, discovered). For `p_bh_lineage`: per family; over all families pooled; and, for the families with no true alternative (every branch outside the foreground), the number with any rejection.

**Criteria (fixed now).** An output **controls** false discoveries in a cell when the upper bound of its FDP interval is at most 0.075 (1.5 times q, the liberal end of Bradley's rule used throughout); it **fails** when the lower bound exceeds 0.05; otherwise **inconclusive**. Families with no alternative are judged by their count with any rejection against 5% of the families (binomial, reported). Power is reported, not judged.

**Precision and cost.** About 100 selected genes and 300 true foreground branch tests per cell; an FDP interval of about ±3-5 points at 5%. aBSREL about 28 and branch models about 36 core-seconds per gene on 8 taxa: 2,000 genes ≈ 36 core-hours, plus about 8 for the diagnostic arm; **ceiling 50 core-hours**. Intel nodes only.

**Implementation.** Units `tools/dev/validation/lineage-mt-units.tsv`, the study B runner, root `study-Q` from a frozen library at the commit adding this design; summary `lineage-mt-summary.R`. Report: `docs/validation/lineage-multiple-testing.md`.

## Stage 2: hypothesis families and candidate corrections (design; approved 2026-11-23, fixed before running)

The maintainer ruled (2026-11-23, `DECISIONS.md`): no global correction is added yet; define the hypothesis families, decide which procedure suits each, and validate empirical FDR and power in mixtures with realistic dependence, then hold any new optional correction pending validation. Correcting more p-values does not fix false positives from alignment error (stage 1: every aBSREL output failed after realignment), so the candidates are judged on error-free alignments, and realigned data are reported but not judged.

### The three families

| Family | Hypotheses | Typical question | Existing output |
|---|---|---|---|
| F-gene | per gene g, H_g: no tested branch of g has sites with ω > 1 | which genes show episodic selection anywhere | gene `p_adj` (BH over genes on HyPhy's per-gene minimum of `p_holm`) |
| F-lineage(b) | for one species-tree branch b chosen in advance, H_{g,b} over genes g | in which genes is lineage b selected | `p_bh_lineage` (BH within b's family) |
| F-all | every gene-by-branch H_{g,b} | which branches of which genes, no lineage chosen in advance | none (stage 1's M3 was exploratory) |

The within-gene family (branches of one gene) is HyPhy's `p_holm`; it answers a question about one gene and is not a cross-gene correction.

### Candidate procedures for F-all (computed from existing p-values; no change to dndsR)

- **P-BH:** one BH over all gene-by-branch `p_raw` (stage 1's M3).
- **P-BY:** Benjamini-Yekutieli over the same, valid under arbitrary dependence, more conservative.
- **P-BB (hierarchical):** Benjamini-Bogomolov two-stage selection. Stage 1: BH at q over genes on the gene `p_value` (the F-gene screen), R genes selected of m. Stage 2: within each selected gene, BH at q·R/m over its branch `p_raw`. Controls the average FDR over selected genes; its overall FDP over all branch discoveries is measured too.
- **P-lineage-pooled (reference):** existing `p_bh_lineage` discoveries summed over lineages (not designed for F-all; stage 1 held only because lineages without selection were silent).

### Dependence

Within a gene, branch tests share one alignment and tree, so they are dependent (realistic, present in every cell). Across genes, the simulator is independent. Selection on varying lineages (more realistic than one fixed foreground) is added by cell **BMIXv**: 10% of genes selected, each on one random branch (probability 0.7) or two (0.3), any lineage; truth per branch from `gene_branch_truth.tsv` (ω mean above 1). Across-gene dependence (for example co-selected pathways) is not simulated and is a stated limit.

### Evaluation data

| Data | Selection | Taxa | Genes | Status |
|---|---|---|---:|---|
| Q-MIX-T8b (stage 1) | fixed foreground | 8 | 1,000 | done |
| Q-MIXh-T8b, true-alignment arm (stage 1) | fixed foreground | 8 | 1,000 | done |
| C-M8, C-M16, arm T (alignment confirmation) | fixed foreground | 8, 16 | 1,000, 500 | running in that study |
| **Q2-V8, Q2-V16** (new, seeds 377, 378) | random lineages (BMIXv) | 8, 16 | 1,000, 500 | to run: aBSREL on true alignments, `lineage-mt2-units.tsv`, root `study-Q2` |
| Realigned: Q-MIXh-T8b, C-M8 and C-M16 arm A0 | fixed foreground | 8, 16 | | reported, not judged |

### Criteria (fixed now)

- Per data set and procedure, at q = 0.05: discoveries, FDP with 95% gene-bootstrap interval (2,000 resamples, seed 1), power (share of true gene-by-branch hypotheses discovered), and for P-BB the average within-selected-gene FDP.
- A procedure **controls** F-all in a data set if its FDP upper bound is at most 0.075; **fails** if the lower bound exceeds 0.05.
- **Proposal rule:** among P-BH, P-BY and P-BB, those that control F-all on every error-free data set; of these, the one with the highest power averaged over the data sets is proposed to the maintainer as a candidate optional column. If none controls everywhere, nothing is proposed and the documentation says F-all has no validated procedure. Nothing is implemented in this stage.
- F-gene (`p_adj`) and F-lineage (`p_bh_lineage` per family) are re-evaluated on the new data with stage 1's criteria.

### Cost

Q2-V8 and Q2-V16: aBSREL about 28 and 64 core-seconds per gene: about 17 core-hours. Everything else reuses existing p-values. **Ceiling 25 core-hours.** Summary: `lineage-mt2-summary.R`, written before the new results are read.
