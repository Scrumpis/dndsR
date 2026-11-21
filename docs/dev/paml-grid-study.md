# Study P: PAML branch models across a parameter grid (design draft)

Status: **draft for the maintainer's approval; not run.** Gap 4 of the [audit](validation-audit.md). The first lineage studies tested the branch models at one divergence, one alignment length (300 codons), and five trees; the omega = 1 boundary at two trees; and starting values only as dndsR's three defaults.

## Hypotheses tested by dndsR (verified in `R/phylo_branch_omega.R`)

| Output | Null | Alternative | Reference distribution |
|---|---|---|---|
| `p_two_vs_m0` | one omega for the tree (M0) | foreground and background omegas differ (two-ratio) | chi-square, df 1 |
| `p_free_vs_m0` | M0 | one omega per branch (free-ratio) | chi-square, df = branches - 1 |
| `p_fg_omega_gt1` | foreground omega = 1 (two-ratio, fixed) | foreground omega > 1 | one-sided signed root: P(Z > sign(omega_fg - 1) sqrt(LRT)), equivalent to a 50:50 mixture of 0 and chi-square(1) |

The nulls are nested in the alternatives with the parameter interior, except omega_fg = 1 under the one-sided alternative (handled by the signed root) and omega estimates at the upper bound 999 (flagged `omega_at_upper_bound`).

## Questions

- **P1, size of `p_two_vs_m0`** with no foreground difference, across divergence, alignment length, and tree.
- **P2, size of `p_fg_omega_gt1`** at omega_fg exactly 1, the boundary where the follow-up gave 5.5% and 6.8% (inconclusive), across divergence and length.
- **P3, power** for omega_fg 0.5 and 2 against omega_bg 0.2 (two-ratio) and omega_fg 2 (above 1).
- **P4, estimation**: bias and spread of omega_fg and omega_bg (two-ratio) and per-branch free-ratio omega, by divergence and length, against the simulated truth; how often omega is at a bound.
- **P5, `p_free_vs_m0` with missing taxa** on 8 taxa (2 of 8 absent per gene), where the follow-up could not establish calibration on 16 taxa (8% then 4%).
- **P6, optimizer**: with 5 random starting omegas (0.01-5, log-uniform) in addition to dndsR's 3, how often a better optimum exists (log-likelihood higher by more than 0.01), and whether it changes a test's decision; starting kappa varied the same way.
- **P7, misspecification**: P1 and P2 with omega varying among sites (classes 0.05/0.3/1), which branch models do not model.

## Cells

Divergence: the tree's branch lengths times 0.25, 1, or 4. Length: 100, 300, 1,000 codons. Trees: T8b, T16, T8p.

| Block | Truth | Grid | Genes per cell | Cells |
|---|---|---|---|---|
| P1 | omega 0.2 everywhere | 3 divergences x 3 lengths x 3 trees | 300 | 27 |
| P2 | omega_bg 0.2, omega_fg 1 | 3 x 3, T8b | 300 | 9 |
| P3 | omega_fg 0.5 or 2 | 3 x 3, T8b | 100 | 18 |
| P5 | omega 0.2, 2 of 8 taxa absent | 3 lengths, T8b, divergence 1 | 300 | 3 |
| P7 | classes 0.05/0.3/1, no foreground difference; and omega_fg classes scaled to mean 1 | 3 lengths, T8b, divergence 1 | 300 | 6 |

P4 uses P1-P3's fits. P6 refits 300 genes from P1 (divergence 4 and length 100, where flat likelihoods are likeliest, and the reference cell) with 5 extra starts.

## Criteria (fixed before running)

As the lineage studies: size at 0.05 by Clopper-Pearson interval, **calibrated** (upper bound at most 5%), **within Bradley's tolerance** (interval inside 2.5-7.5%), **anticonservative** (lower bound above 5%), otherwise inconclusive. A test is validated for a cell's conditions if calibrated or within tolerance there. Power reported with no threshold. P4: median relative bias within plus or minus 10% for omega_bg and omega_fg at 300 codons and divergence 1 (the operating point the documentation describes), reported everywhere else. P6: a better optimum in at most 1% of genes, and no change of decision at 0.05 in more than 0.5%; otherwise the default number of starts is reported to the maintainer as insufficient (defaults are not changed without approval).

## Precision

300 genes per null cell: an interval of about plus or minus 2.5 points at 5%, and P2 pooled over the 9 cells (2,700 genes) about plus or minus 0.8 points, enough to tell 6% from 5%.

## Seeds and cost

New seeds 501-563, one per cell. Measured: M0 + two-ratio + null about 36 core-seconds per gene on 8 taxa at 300 codons; free-ratio adds about 100. Scaled by length (roughly linear) and taxa (T16 about 3 times): P1 about 120 core-hours, P2 about 35, P3 about 25, P5 about 35, P7 about 20, P6 about 15: **about 250 core-hours**. A pilot of 20 genes per cell (about 20 core-hours) checks the scaling of divergence and length in the simulator, and the timings, before the full runs. Simulator: `phylo-lineage-sim2.R` with a divergence factor and the codons argument, as in study R-D.

## Not covered

Branch-site models (not in dndsR), codon-frequency model choice (dndsR fixes F3x4), alignment and orthology error (study B).
