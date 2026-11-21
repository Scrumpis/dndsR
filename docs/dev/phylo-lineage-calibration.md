# Lineage-level methods: calibration design

Fixed on 2026-11-20, before any simulation was analysed. Changes after that are listed in the change log at the end, with reasons. Methods: `phylo_branch_omega` (PAML codeml branch models), `phylo_absrel` (HyPhy aBSREL), `phylo_relax` (HyPhy RELAX), and the pooled omega of `phylo_summarize`. Software correctness (dndsR reproduces the engines run directly) is checked separately, in tests and `docs/validation/phylo-lineage.md`; this study asks whether the outputs mean what they say on data with a known truth.

## Questions

1. **Size.** Under each method's null, how often is a test significant at 0.05, per test and after dndsR's corrections?
2. **Power.** Under alternatives of the kind the methods target, how often are true effects found?
3. **Estimation.** How close are per-branch omega estimates, and the pooled lineage omega, to the truth, and do dndsR's reliability flags mark the bad ones?
4. **Robustness.** What changes with tree shape, branch length (very short and very long branches), denser sampling, missing taxa per gene, and site-to-site omega variation that branch models ignore?

## Simulator

`tools/dev/validation/phylo-lineage-sim.R`, built on the codon simulator of the BUSTED calibration (`phylo-sim.R`: HKY mutation, F3x4 codon frequencies and kappa estimated from 9,084 real Brassicaceae orthogroups, `params.json`), extended so that every branch can have its own omega classes. Independent of dndsR's package code. Branch lengths are given in the units of the BUSTED calibration and multiplied by its factor `c`; 300 codons per gene unless stated.

## Trees

| Name | Shape | Purpose |
|---|---|---|
| T5 | The empirical five-taxon Brassicaceae tree of the BUSTED calibration | Continuity with that study |
| T8b | Balanced, 8 taxa, terminal 0.03, internal 0.02 | Typical case |
| T8p | Pectinate, 8 taxa, terminal lengths from 0.005 to 0.25 | Very short and very long (saturated) branches |
| T16 | Balanced, 16 taxa, terminal 0.03, internal 0.015 | Denser sampling, more parameters |
| T16m | T16 with a random quarter of taxa removed from each gene (at least 6 kept) | Missing taxa and pruning |

## Cells

Null cells have 300 genes, alternative cells 100. Foreground for two-ratio and RELAX: the clade of the first two taxa and its stem. aBSREL tests every branch.

| Cell | Truth | Trees | Measures |
|---|---|---|---|
| A0 | One omega, 0.2, on every branch and site | T5, T8b, T8p, T16, T16m | Size of `p_two_vs_m0`, `p_free_vs_m0`, `p_fg_omega_gt1`; per-branch omega error by flag |
| A1 | Foreground 0.6, background 0.2 | T5, T8b, T16 | Power of `p_two_vs_m0`; error of `omega_fg`, `omega_bg` |
| A2 | Each branch its own omega, drawn log-normal around 0.2 (sd of log 0.7), fixed per tree | T8b, T8p, T16 | Per-branch omega error; pooled lineage omega against the branch's true omega |
| A3 | Sites in three classes (omega 0.05, 0.3, 1; 60/30/10%), the same on every branch | T5, T8b | Size of the branch-model tests when sites vary |
| A4 | Foreground 1.5 at every site, background 0.2 | T8b | Power of `p_fg_omega_gt1` |
| B0 | Sites 0.05/0.3/1 (60/30/10%), no omega above 1 anywhere | T5, T8b, T8p, T16 | Per-branch size of `p_raw`; per-gene size of `p_value` (Holm minimum) |
| B1 | As B0, but 10% of sites at omega 10 on one terminal branch | T8b, T16 | Power on that branch (`p_holm`, `p_bh_lineage`) |
| C0 | Same omega classes on test and reference (k = 1) | T5, T8b, T16 | Size of RELAX |
| C1 | Test classes are reference classes raised to k = 0.5 (relaxation) | T8b, T16 | Power, and the sign of k |
| C2 | k = 2 (intensification) | T8b, T16 | Power, and the sign of k |

Reference classes for RELAX: 0.05, 0.3, 1.5 (60/30/10%), so that both relaxation and intensification change them.

## Criteria (fixed in advance)

For each null cell and test, the false-positive rate at 0.05 with a 95% Clopper-Pearson interval:

- **Calibrated**: upper bound at most 0.05.
- **Consistent with nominal**: estimate at most 0.05, upper bound above it.
- **Anticonservative**: lower bound above 0.05. Reported as a limitation in the user documentation, with the conditions.

Estimation (A cells): median absolute log2 error of omega, overall and by dndsR's reliability flag. The flags are useful if flagged estimates have clearly larger errors than unflagged ones (median absolute log2 error at least twice as large). The pooled lineage omega (A2) is judged by its median absolute log2 error against the branch's true omega, compared with the median of per-gene estimates.

Power is reported, not judged: there is no fixed acceptable power.

A method is described in the documentation as calibrated only for the conditions where it met the size criterion, and only in simulation.

## Budget and order

A smoke test of 10 genes per cell first, to measure run time. The full study is launched only if its projected cost stays under 250 core-hours; otherwise cells are trimmed in this order: T16m, then T8p for B0, then A4. Results root outside the repository; summaries and the report come into `docs/validation/phylo-lineage-calibration.md`.

## Change log

- 2026-11-20: smoke test, 10 genes per cell, all 27 units at `f7aeaea`, no failures. Projected cost of the full study from its timings: about 246 core-hours, dominated by free-ratio fits on 16 taxa (about 15 core-minutes per gene), under the 250 core-hour limit, so the full study was launched as designed (array 21467783), at the same commit, with no cells trimmed.
- 2026-11-20: the full study filled the scratch file quota within about 20 minutes (1.0 million files): `phylo_branch_omega` kept about 140 files per orthogroup (each start's copies of the inputs, and codeml's auxiliary outputs). The array was cancelled, codeml's auxiliary files from this work were deleted (about 140,000), and the partial units, at `f7aeaea`, were discarded unread. The package was changed to keep about 25 files per orthogroup, with `keep_runs = FALSE` keeping none; the driver uses `keep_runs = FALSE` and archives HyPhy's folders per unit. No numerical code changed. The study restarts at the commit that makes this change; the smoke test's timings stand.
- 2026-11-21: the full study completed at `dcf3231` (array 21469993): 27 of 27 units, 0 genes failed, 234.6 core-hours. The summary script, committed at `f7aeaea` before any result was read, has a fourth verdict this design did not define, *above nominal, not significantly* (estimate above 0.05, interval including it); the report uses it as such. The design does not measure the size of `p_fg_omega_gt1` at its boundary (foreground omega 1); the report says so rather than adding a cell after seeing results. Report: [validation/phylo-lineage-calibration.md](../validation/phylo-lineage-calibration.md).
