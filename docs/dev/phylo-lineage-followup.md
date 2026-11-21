# Lineage-level methods: follow-up study design

Fixed on 2026-11-21, after the results of the first study ([phylo-lineage-calibration.md](phylo-lineage-calibration.md), report [validation/phylo-lineage-calibration.md](../validation/phylo-lineage-calibration.md)) had been read, and before any simulation below was run. It is a separate study: its questions come from those results, so it uses new seeds and states its hypotheses, nulls, and criteria here. The maintainer asked for it (decisions in `DECISIONS.md`, 2026-11-21). Changes after this point go in the change log at the end, with reasons.

## Audit of the first report

Before new work, each claimed verdict was checked against the prespecified criterion (upper 95% Clopper-Pearson bound at most 0.05 for *calibrated*) and the replication counts. One claim was wrong: aBSREL's gene-level `p_value` meets it only on T5 (5/300, upper 3.8%); on T8b, T8p, and T16 it is *consistent with nominal* (12/300, 8/300, 9/300; upper bounds 6.9%, 5.2%, 5.6%), not calibrated. The report, methods page, NEWS, and skill said "calibrated on every tree" and are corrected. *Consistent with nominal* means only that the data do not show a rate above 0.05; it is not evidence that the rate is at most 0.05. `p_fg_omega_gt1` met the criterion, but only at a foreground omega of 0.2, far below its null boundary.

## Tests and their nulls (checked in `R/phylo_branch_omega.R`)

| Output | Alternative model | Null model | Hypotheses | Reference distribution |
|---|---|---|---|---|
| `p_two_vs_m0` | two-ratio (omega_bg, omega_fg free) | M0 (one omega) | H0: omega_fg = omega_bg; H1: they differ (two-sided) | chi-square, 1 df (interior null) |
| `p_fg_omega_gt1` | two-ratio | two-ratio with omega_fg fixed at 1 (`fix_omega = 1`) | H0: omega_fg <= 1; H1: omega_fg > 1 (one-sided) | signed root `sign(w_fg - 1) * sqrt(LRT)` is N(0, 1) at omega_fg = 1, an interior point of (0, inf); p = P(Z > root). Equal to PAML's convention of halving the chi-square(1) p-value when w_fg > 1. Exact at the boundary, conservative below it |
| `p_free_vs_m0` | free-ratio (one omega per branch) | M0 | H0: one omega on all branches | chi-square, branches minus 1 df |

The two branch-model tests answer different questions: a foreground omega of 1 against a background of 0.2 is a large *difference* (`p_two_vs_m0`) and no evidence of omega above 1 (`p_fg_omega_gt1`). The fixed-omega null is fitted from one start; an under-optimized null would inflate its LRT, so its convergence is examined (part D).

## Criteria

Size at 0.05 for each null cell, with a 95% Clopper-Pearson interval, judged two ways:

1. The first study's categories (*calibrated*: upper bound at most 0.05; *consistent with nominal*: estimate at most 0.05, upper bound above; *anticonservative*: lower bound above 0.05; otherwise *above nominal, not significantly*). Under this rule a test exactly at its level is *consistent with nominal* at best, so it cannot alone show calibration at a boundary.
2. Bradley's liberal criterion, added for this study because of that: **within tolerance** if the whole interval lies in [0.025, 0.075]; **outside tolerance** if the interval lies wholly above 0.075 or below 0.025; otherwise **inconclusive**. With 600 genes a rate of 0.05 gives an interval of about 3.4-7.1%, so cells meant to show calibration at the boundary have 600 genes.

Rates at 0.01 are reported, not judged.

## Part D: foreground omega exactly 1

- **Cell D0**: foreground (the first two taxa and their stem) omega 1 at every site; background 0.2; no site variation. Trees T5 (300 genes), T8b (600), T16 (600), T16m (300; foreground taxa can be pruned, and genes where the foreground is not identifiable are excluded by dndsR). Models M0, two-ratio, two-ratio null; seeds 101-104.
- **Hypothesis D**: at omega_fg = 1, `p_fg_omega_gt1` has size 0.05. Criterion: within Bradley's tolerance on T8b and T16. Also reported: the two-sided size of `lrt_fg_omega_eq1` against chi-square(1); the fraction of genes with estimated omega_fg above 1 (expected near one half); mean LRT (expected 1); how often `p_two_vs_m0` is significant (the difference test, which should be powerful here); and, for genes with p at most 0.05, whether the null fit's likelihood is consistent with the alternative (`lrt_negative`, `codeml_warning`).

## Part R: RELAX's false-positive rate

First study: 9.0% (27/300) on T16, and LRTs above chi-square(1) on every tree (mean 1.15-1.28, expected 1), each gene fitted from one starting point (HyPhy's default, dndsR's default). Four candidate causes, each with its own check:

| Part | Question | Data | Criterion |
|---|---|---|---|
| R1 | Does HyPhy, run without dndsR, give the same result? | The first study's 300 C0 T16 alignments; tree written by hand; HyPhy with dndsR's settings and seed (`tools/dev/validation/phylo-lineage-relax-direct.sb`) | LRT equal to dndsR's within 1e-3 for every gene |
| R2 | Is it optimization? | The same 300 genes with 10 starting points (direct); and new genes (R3) refitted with `starting_points = 10` through dndsR | Optimization is implicated if more starts lower the LRT of the false positives (change in null and alternative log-likelihood reported) and bring size within Bradley's tolerance while one start does not |
| R3 | Does it survive a new seed? | Cell C0 on T16, 600 new genes, seed 105, one start | Reproduced if the size is outside Bradley's tolerance or anticonservative |
| R4 | Is it the simulator? The first study's simulator draws each site's omega class once for all branches; HyPhy's RELAX assumes a class drawn independently on each branch | Cell **C0i**: C0 with classes drawn per branch, T16, 600 genes, seed 107, one start | Misspecification is implicated if C0i is within tolerance while R3 is not |

Sample size is a fifth: with 3 test branches of 29, the test set carries few substitutions, and the chi-square approximation may be poor whatever the cause; the R3 and R4 cells measure the size directly at 600 genes. The thresholds of the first study are not changed for any outcome.

## Part F: free-ratio test with missing taxa

First study: `p_free_vs_m0` 8.0% (24/300) on T16m, 5.0% on T16; mean LRT/df 1.05 on T16m, variance 1.17 times the chi-square's. Not associated with codeml warnings, disagreeing starts, merged branches, or flagged branches. A subgroup (genes missing a foreground taxon) had a higher rate, but on the balanced T16 with one omega every tip is equivalent, so it is chance among several comparisons; it is not tested further.

- **Cell A0 on T16m, replicate**: 300 new genes, seed 108, models M0 and free-ratio.
- **Hypothesis F**: the excess is reproducible. Judged on the replicate alone (first study's categories and Bradley's), and on the 600 genes pooled. Diagnostics reported: mean LRT/df, its variance against 2 df, and the size by number of merged branches (the one prespecified subgroup).
- Whatever the outcome, `p_free_vs_m0` is a test of any variation among branches, with many parameters; its per-branch omegas are descriptive estimates, not tests.

## Budget

A smoke test of 10 genes per unit first. The whole follow-up is limited to 120 core-hours, projected from the smoke test; over that, cells are trimmed in the order: D0 T5, D0 T16m, the dndsR 10-start refit (R2 keeps the direct 300 genes).

## Outputs and labels

Results root outside the repository. Summary: `tools/dev/validation/phylo-lineage-followup-summary.R`. Report: appended to [validation/phylo-lineage-calibration.md](../validation/phylo-lineage-calibration.md) as a separate section. Whatever the results, outputs whose calibration failed are flagged in the tables dndsR writes, and the documentation says where.

## Change log

- 2026-11-21, smoke test (10 genes per unit, at `c9123a5`; read for run time and failures only, its results discarded). Two changes to part R, before the full runs:
  - HyPhy's `--starting-points 10` stopped with an error ("Matrix/List LHS expected") on 2 of 10 genes, run directly and through dndsR alike, and on others reached a *lower* likelihood than one starting point (for example the null of gene 1: -4600.15 against -4599.37). It is not usable as an optimization check in HyPhy 2.5.103. R2 instead fits each gene with three random seeds (1, 2, 3), one starting point each, and takes the best likelihood per model, as dndsR does for BUSTED with `extra_seeds`.
  - Run directly on the same alignments with the hand-written tree, HyPhy's LRTs differed from dndsR's by up to 0.1 (gene 8: 1.86 against 1.76): RELAX's optimum depends on how the same data and tree are written. R1 is therefore split: **R1**, HyPhy directly on the alignment and tree files dndsR wrote (archived), must equal dndsR's result within 1e-3, checking dndsR's handling; **R1'**, the hand-written inputs, measures the sensitivity to input representation, reported, not judged.
- Budget: projected from the smoke test (wall time times 8 cores, scaled to the gene count), 192 core-hours with the dndsR 10-start refit, over the 120 limit. Trimmed in the prespecified order: D0 T5, D0 T16m, the dndsR 10-start refit (which is also unusable, above). The remainder projects to about 160 core-hours by that method, but it overestimates: for A0 T16m with M0 and free-ratio only, it projects 20,220 s of wall time for 300 genes, while the first study's A0 T16m took 13,952 s for 300 genes with all four models, a factor of at least 1.45 (10 genes on 8 cores run in two uneven waves). Corrected by 1.45, about 110. A watcher cancels the runs if the actual cost reaches 120 core-hours. Without D0 T16m, the boundary is not tested with missing taxa.
