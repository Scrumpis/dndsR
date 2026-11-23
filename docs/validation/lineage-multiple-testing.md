# Multiple testing in the lineage outputs: mixture study (results)

Status: **follow-up to study B, run once to the design fixed before running ([design and audit of definitions](../dev/lineage-multiple-testing.md), commit `48aa75d`).** It answers the maintainer's item 3 of 2026-11-22. Criteria were not changed after results. Recommendations at the end are for the maintainer to decide; no code or default has changed because of them.

## Runs

| | |
|---|---|
| Frozen library | commit `48aa75d`, image `dndsr_v0.4.0_sha256-f1b42d9d.sif`, HyPhy 2.5.103, PAML 4.10.10 |
| Q-MIX-T8b | cell BMIX without corruption, 8 taxa, 1,000 genes, seed 331; 111 genes selected on t1, t2 and their stem (classes 0.05/0.3/10); aBSREL and branch models on the true alignments |
| Q-MIXh-T8b | cell BMIXh (BMIX's selection with BE3h's homologous indels), seed 332, 1,000 genes, 87 selected; realigned by `phylo_prepare`; aBSREL and branch models |
| Q-MIXh-T8b, true alignment | the same 1,000 genes, aBSREL on their true alignments (diagnostic arm) |
| Nodes | Intel Xeon Gold 6148 only; Slurm 21579083 (2 tasks), 21581106 |
| Cost | 47.2 core-hours allocated (ceiling 50) |
| Summary | `tools/dev/validation/lineage-mt-summary.R`; tables `data/lineage_multiple_testing_{fdr,families}.tsv` |

Truth: a gene is selected or not (`genes.tsv`); a branch test is true only on a foreground branch (from `sim/branch_truth.tsv`) of a selected gene. Each selected gene has 3 true branch tests and 10 null ones.

## False-discovery proportion and power at 0.05

FDP with a 95% gene-bootstrap interval (2,000 resamples, seed 1). **Controls**: upper bound at most 0.075; **fails**: lower bound above 0.05; otherwise inconclusive. Power is the share of true genes, or true foreground branch tests, discovered.

| Output | Error-free (Q-MIX) | Realigned indels (Q-MIXh) | Same genes, true alignment |
|---|---|---|---|
| aBSREL `p_bh_lineage`, all families pooled | 181 disc., FDP 1.1% (0-2.9), **controls**; power 53.8% | 248, 46.8% (40.0-54.1), **fails**; 50.6% | 119, 0.8% (0-2.9), **controls**; 45.2% |
| aBSREL gene `p_adj` | 86, 2.3% (0-6.0), **controls**; 75.7% | 170, 58.2% (50.9-65.4), **fails**; 81.6% | 59, 3.4% (0-8.6), inconclusive; 65.5% |
| M3 (exploratory): one BH family over all branch tests | 129, 2.3% (0-5.4), controls; 37.8% | 231, 55.0% (48.0-62.3), fails; 39.8% | 84, 2.4% (0-6.1), controls; 31.4% |
| aBSREL `p_holm` per branch (a within-gene correction, for context) | 202, 6.9% (3.6-10.9), inconclusive; 56.5% | 365, 60.0%, fails; 55.9% | 161, 11.2% (6.3-16.9), fails; 54.8% |
| aBSREL gene `p_value`, unadjusted (for context) | 118, 10.2% (5.1-16.2), fails; 95.5% | 262, 68.7%, fails; 94.3% | 96, 15.6%, fails; 93.1% |
| aBSREL `p_raw`, unadjusted (for context) | 520, 46.0%, fails; 84.4% | 996, 78.3%, fails; 82.8% | 465, 54.4%, fails; 81.2% |
| branch models `p_adj_two_vs_m0` | 109, 3.7% (0.8-7.7), inconclusive; 94.6% | 69, 1.4% (0-4.9), **controls**; 78.2% | (not run) |
| branch models `p_adj_fg_omega_gt1` | 0 discoveries; power 0% | 0; 0% | (not run) |

The unadjusted outputs and `p_holm` are not designed to control the false-discovery rate across genes; they are shown so the effect of each correction is visible.

### `p_bh_lineage` per family

| Lineage family | Error-free | Realigned indels | True alignment |
|---|---|---|---|
| t1 (true alternative) | 71 disc., 1 false, FDP 1.4% (0-4.6), controls; power 63% | 67, 13 false, 19.4% (10.1-30.0), fails; 62% | 46, 1 false, 2.2% (0-7.5), inconclusive; 52% |
| t2 (true alternative) | 68, 0 false, controls; 61% | 68, 10 false, 14.7% (6.8-23.8), fails; 67% | 55, 0 false, controls; 63% |
| n3, their stem (true alternative) | 42, 1 false, 2.4% (0-8.1), inconclusive; 37% | 23, 3 false, 13.0% (0-27.8), inconclusive; 23% | 18, 0 false, controls; 21% |
| 10 families with no alternative: any rejection | 0 of 10 (95% CI 0-31%) | **10 of 10** (69-100%), 4-17 false discoveries each | 0 of 10 (0-31%) |

## Answers to the design's questions

- **M1, error-free alignments.** `p_bh_lineage` controls false discoveries per lineage family and when its discoveries are pooled across families. Pooling separate BH families carries no general guarantee; it held here because the null lineages made no rejections. Gene `p_adj` controls false discoveries over genes. The branch models' `p_adj_two_vs_m0` is inconclusive (point estimate 3.7%, upper bound 7.7%, just over the 7.5% limit). `p_adj_fg_omega_gt1` made no discoveries, because selection on 10% of sites does not raise the branch's mean ω above 1.
- **M1, realigned indels.** Every aBSREL output fails. About half of `p_bh_lineage`'s and `p_adj`'s discoveries are false, and every lineage without selection rejects something. On the true alignments of the same genes, both corrected outputs return to near 1-3%. The failure comes from realignment error making aBSREL's `p_raw` invalid, not from the corrections. BH cannot correct invalid p-values. The branch models' `p_adj_two_vs_m0` controls under the same error (1.4%), with power reduced from 95% to 78%.
- **M2, power.** At 8 taxa, `p_bh_lineage` finds 61-63% of true terminal-branch tests and 37% of the stem's. Gene `p_adj` finds 76% of selected genes, and `p_adj_two_vs_m0` 95%. The branch models are told the foreground in advance; aBSREL searches every branch.
- **M3 (exploratory).** One BH family over all 13,000 branch tests controls false discoveries on error-free data at 2.3%, but finds 38% of true branch tests against 54% for the per-lineage families. That is the price of not choosing lineages in advance.

## Which quantity answers which question

Proposed wording for the user documentation, conditional on correct alignments:

| Question | Output | Correction across genes |
|---|---|---|
| Is this gene under episodic selection on some branch? (gene level) | aBSREL gene `p_adj` | included (BH over tested genes); `p_value` is uncorrected across genes |
| Is the pre-specified foreground lineage faster-evolving in this gene? (gene level, foreground chosen in advance) | branch models `p_adj_two_vs_m0` | included |
| In which genes is lineage X selected, for a lineage chosen in advance? (branch level) | `p_bh_lineage` within lineage X's family | included (BH within the family) |
| Which branches of which genes, with no lineage chosen in advance? (exploratory) | no current output guarantees it. Pooled `p_bh_lineage` held here only because null lineages were silent; M3, one BH family over all branch tests, is the guaranteed option (not implemented) | — |
| Is this branch of this gene selected, as a single planned test? | `p_raw` | none; use only for one pre-specified test |

`p_holm` corrects only across branches within one gene; it is not a cross-gene correction. The counts in `absrel_n_sig_lineage` use `p_bh_lineage`, so they are valid per lineage on error-free alignments, but summing them across lineages is exploratory.

## Limits

One cell design per condition: 8 taxa, 300 codons, 10% of genes selected on 3 of 13 branches, one indel process, 1,000 genes. Other selection rates, tree sizes or patterns of dependence between genes may differ. The branch models were not run on the true alignments of Q-MIXh. Masking was not run (it was a separate sensitivity analysis in study B). This is one run per cell.

## Recommendations for the maintainer (not implemented)

1. **Documentation only (no behaviour change).** Replace "no output is recommended for lineage counts yet" with the table above; state that every aBSREL output, corrected or not, lost false-discovery control on indel-realigned alignments, while the branch models kept it; keep the alignment-check warning's wording that an unflagged gene is not shown to be aligned correctly.
2. **Optional new output (a behaviour change, needs approval).** Add M3, one BH family over all branch tests of all genes, as a new column (for example `p_bh_all_branches`) for exploratory analyses without a pre-chosen lineage. It would change no existing column.
3. **Rename nothing, change no default.** The definitions of `p_holm`, `p_value` and `p_adj` stay as they are. Their documentation is corrected to say that HyPhy's `p_holm` omits Holm's running maximum, and that the gene `p_value` is Hochberg's global test rather than Bonferroni's.
