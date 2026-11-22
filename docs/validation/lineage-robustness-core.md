# Lineage tests under input errors: aBSREL and branch models (study B, core arms)

**Question.** Which input errors break aBSREL (`phylo_absrel`) and the branch models (`phylo_branch_omega`)? Which errors can be detected from the input? For each, should dndsR warn, exclude, mask, or document a limitation? And is a failure upstream (HyPhy or codeml on that input) or a dndsR defect? Design, criteria and classification rules, fixed before running: [docs/dev/lineage-robustness-study.md](../dev/lineage-robustness-study.md).

**Scope of this report.** These are the core arms only: every aBSREL and branch-model arm, as given and masked, on the 8-taxon tree and on 16 taxa with 4 missing. The maintainer has not yet decided on the 16-taxon units and the BUSTED arms (cost, below), so nothing here applies to BUSTED, or to 16 complete taxa.

**Answer.**
- **A corrupted segment in one sequence** (E1 foreground, E2 background) breaks both methods. It is detected in every affected gene by the window profile dndsR already reports (`dissimilar_windows`), which flags none of the error-free genes. Masking the flagged windows restores calibration.
- **Realigned indels** inflate aBSREL's raw branch p-values and its gene p-value. `p_bh_lineage` does not control false discoveries here either: see the correction below. Masking does not fix this, and no signal detects it. With correctly simulated (homologous) insertions, BE3h, the excess is smaller (`p_raw` 7.1%, gene `p_value` 21%) and disappears on the true alignment (2.3%, 3.0%): it comes from alignment error introduced by realignment. E3's larger excess also included a simulator defect (below).
- **Paralogs, gene-tree discordance, missing taxa and recombination** (E4-E7) leave both methods within or below 5%. For discordance and paralogs this holds partly because dndsR's existing zero-length-foreground exclusion removes many of the affected genes.
- **Engine or dndsR:** wherever size was outside tolerance, HyPhy and codeml run directly on dndsR's files reproduce dndsR's statistics. These failures belong to the models on that input, not to dndsR.

## Runs

300 genes per cell (600 for E1), 300 codons, site classes 0.05/0.3/1 at 0.6/0.3/0.1, seeds 301-312 (E4 at its re-piloted generator). Code at `37aed18` (frozen library in `study-B2/lib`), image `dndsr_v0.4.0_sha256-f1b42d9d`, HyPhy 2.5.103, codeml 4.10.10. Units: `tools/dev/validation/phylo-study-B-units-core.tsv` and `phylo-study-B-units-core-E4.tsv`; BE3h at `49fd39f`, `phylo-study-B-units-E3h.tsv`. About 130 core-hours, plus about 13 for BE3h, its direct reruns and the true-alignment diagnostic. Summary: `phylo-study-B-summary.R`; mixture arms: `phylo-study-B-exclude.R`; direct reruns: `phylo-study-B-direct.R`. Tables: `data/lineage_robustness_core_*.tsv`.

| Cell | Error |
|---|---|
| E0 | none |
| E1 | one corrupted segment (stage 2's empirical lengths) in a foreground sequence (t1); 600 genes |
| E2 | the same in a background sequence (t8) |
| E3 | codon indels (rate 0.02, mean 2 codons), realigned by `phylo_prepare`; insertions on internal branches were not homologous across descendants (simulator defect, below) |
| BE3h | E3 rerun with homologous insertions (`add_indels(homologous = TRUE)`, code `49fd39f`, same seed 304), aBSREL and branch models |
| E4 | t1 replaced by a paralog that split at half the tree's depth |
| E5 | discordance: t2 and t3 exchanged in the gene tree |
| E6 | 16 taxa, 4 missing per gene |
| E7 | recombination: the second half of each gene on E5's gene tree |
| MIX | 1,000 genes, 10% selected on the foreground (classes 0.05/0.3/10), each taxon corrupted at stage 2B's empirical rate |

## Size at 0.05, as given (masked in brackets where it differs)

Categories as the lineage studies: *calibrated* (upper bound ≤ 5%), *within* Bradley's 2.5-7.5%, *anticonservative* (lower bound > 5%), otherwise *inconclusive*.

| Cell | aBSREL `p_raw`, per branch | aBSREL `p_bh_lineage`, per branch | aBSREL gene `p_value` | `p_two_vs_m0` | `p_fg_omega_gt1` |
|---|---|---|---|---|---|
| E0 | 2.1% calibrated | 0% calibrated | 2.3% calibrated | 4.3% inconclusive | 0% calibrated |
| E1 | **9.5% anticons.** [2.1% calibrated] | **7.7% anticons.** [0%] | **99.7% anticons.** [1.0% calibrated] | **45% anticons.** [3.3% inconclusive] | 0% calibrated |
| E2 | **9.9% anticons.** [2.2% calibrated] | **7.7% anticons.** [0%] | **100% anticons.** [1.7% calibrated] | **9.3% anticons.** [5.3% inconclusive] | 0% calibrated |
| E3 | **14.4% anticons.** [13.8%] | 5.1% within [4.6% within] | **48% anticons.** [44%] | 5.7% inconclusive | 0% calibrated |
| BE3h | **7.1% anticons.** [7.1%] | 0.8% calibrated [0.8%] | **21% anticons.** [21%] | 5.3% inconclusive [5.3%] | 0% calibrated |
| E4 | 2.7% calibrated | 0% calibrated | 2.3% inconclusive (213 tested) | 4.3% inconclusive | 0% calibrated |
| E5 | 1.4% calibrated | 0.2% calibrated | 4.1% inconclusive (49 tested) | 2.7% inconclusive | 0% calibrated |
| E6 | 2.1% calibrated | 0% calibrated | 1.3% calibrated | 4.3% inconclusive (279) | 0% calibrated |
| E7 | 2.7% calibrated | 0% calibrated | 3.7% inconclusive | 4.3% inconclusive | 0% calibrated |

Intervals are in `data/lineage_robustness_core_size.tsv`. E0's `p_two_vs_m0` (4.3%, 2.3-7.3) cannot be separated from 5%, so for the branch models E3-E7 are judged against that same baseline. None is anticonservative.

**Exclusions that already exist.** dndsR excludes an orthogroup whose foreground branch has zero length on the fitted tree (`foreground_branch_zero_length`). That rule removed 251 of 300 E5 genes from aBSREL: when the gene tree disagrees, the foreground stem collapses. It also removed 87 of 300 E4 genes. The branch models tested all of them. These exclusions are listed with their reason in `absrel_results.tsv`, as documented.

## Detection

Share of genes each signal flags, against the error-free E0 (criterion for *detectable*: at least 80% of affected genes, at most 5% of E0):

| Signal | E1 | E2 | E3 | BE3h | E4 | E5 | E7 | E0 | Mixture: sensitivity / false flags on unaffected |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| window profile (`dissimilar_windows`) | **100%** | **100%** | 7% | 0.3% | 0% | 0% | 0% | 0% | **99.5% / 0%** |
| any `qc_flags` | 100% | 100% | 17% | 9% | 3% | 6% | 2% | 0.7% | 99.5% / 2.3% |
| a tip's dS > 3× the gene's median | 17% | 19% | 0% | 0% | **100%** | 0% | 0% | 0% | 26% / 0% |
| `few_synonymous` | 0% | 0% | 0% | 0% | 67% | **99.7%** | 2% | 0% | 0% / 0% |
| `high_ds` | 0.7% | 0% | 0% | 0% | 0% | 0% | 0% | 0% | 0.2% / 0% |

E1 and E2 are detectable by the window profile; E4 by the tip-dS ratio; E5 by `few_synonymous` (in the branch-model fit on the species tree). E3, BE3h and E7 are not detectable by any signal. The topology comparison of `gene_tree_mode = "report"` was not run in these arms, which take the species tree.

## Masking and exclusion in the mixture

1,000 genes; 428 corrupted, 116 selected. False discovery proportion after each output's documented correction (BH across genes; for aBSREL branches, BH per species-tree branch, `p_bh_lineage`), and power on the selected tests. In the exclusion arm, genes flagged by the window profile are removed, and power still counts every selected test, so an excluded gene counts as not found. The power difference is against the as-given arm, with a gene-bootstrap 95% interval. Criterion for masking or exclusion: loss under 5 points, lower bound above −10.

| Test | Arm | Discoveries | FDP (95% CI) | Power | Power difference (95% CI) |
|---|---|---:|---:|---:|---:|
| aBSREL per branch | as given | 722 | **74.4%** (71.0-77.5) | 53.2% | |
| | masked | 192 | 5.2% (2.5-9.4) | 52.3% | −0.9 (−4.5 to +1.1), **non-inferior** |
| | flagged excluded | 117 | 5.1% (1.9-10.8) | 31.9% | −21.3 (−28.5 to −14.4), not |
| aBSREL gene `p_value` | as given | 494 | **78.3%** (74.4-81.9) | 92.2% | |
| | masked | 89 | 7.9% (3.2-15.5) | 70.7% | −21.6 (−30.1 to −10.1), not |
| | flagged excluded | 53 | 5.7% (1.2-15.7) | 43.1% | −49.1 (−58.9 to −36.4), not |
| `p_two_vs_m0` | as given | 116 | **25.0%** (17.4-33.9) | 75.0% | |
| | masked | 99 | 4.0% (1.1-10.0) | 81.9% | +6.9 (+2.6 to +11.9), **non-inferior** |
| | flagged excluded | 62 | 4.8% (1.0-13.5) | 50.9% | −24.1 (−32.0 to −16.3), not |
| `p_fg_omega_gt1` | as given, masked | 0 | — | 0% | |

*Corrected 2026-11-22:* the first version of this table counted every branch of a selected gene as a true test for aBSREL per branch (the branches table's `foreground` column is TRUE for every tested branch when all branches are tested). The true tests are the foreground branches (t1, t2, their stem) of selected genes, 348 rather than 1,508. That raised the as-given FDP from 66.3% to 74.4%, and power from 16.1% to 53.2% (masked 12.1% to 52.3%; excluded 7.4% to 31.9%). The classifications did not change. The gene-level and branch-model rows were unaffected. Fixed in `phylo-study-B-summary.R` and `phylo-study-B-exclude.R`, which now read the truth from `sim/branch_truth.tsv`.

The as-given arm's power is inflated by corruption, which also produces true discoveries in corrupted selected genes. So the power differences overstate what masking costs on clean data. Read them as the non-inferiority rule specifies, not as clean-data power. `p_fg_omega_gt1` found nothing: with 10% of sites at ω = 10, the foreground's whole-gene ω is about 1.1, too close to 1 for 300 codons.

## Upstream model or dndsR

For every method and cell outside tolerance, the 20 most significant genes were rerun directly, without dndsR code. aBSREL ran on the alignment and tree dndsR wrote. codeml ran on a PHYLIP file written from that alignment, a hand-written foreground-labelled tree, and a control file with the provenance's settings, keeping the best of three ω starts.

| Cell, run | Statistics compared | Agree (LRT within 1e-3 or relative 1e-4) | Largest difference |
|---|---:|---:|---:|
| E1 aBSREL | 260 branch LRTs | 260 | 3e-13 |
| E2 aBSREL | 260 | 260 | 5e-13 |
| E1 `p_two_vs_m0` | 20 | 20 | 5e-14 |
| E2 `p_two_vs_m0` | 20 | 20 | 4e-15 |
| E3 aBSREL | 260 | 253 | 0.08 |
| E3 aBSREL masked | 260 | 249 | 0.08 |
| BE3h aBSREL | 260 | 260 | 5e-14 |

On E3, 7 and 11 branch LRTs differed, in 3 and 4 genes. Rerun three times on the same node, HyPhy gave identical results each time, so the differences are not run-to-run noise. The E3 unit had run on an AMD node (amr-191); the direct reruns ran on Intel Skylake (skl-076). Rerun directly on an AMD node (amr-215), HyPhy reproduced dndsR's values for all four genes (for example t1 of `BE3_T8b_0282`: 10.3212 against dndsR's 10.3212, and 10.4020 on Skylake). **The differences come from HyPhy behaving differently on the two CPU types, not from dndsR.** They appear only in gapped alignments here, and are at most 0.08 in LRT. It is a reproducibility caveat for HyPhy results across clusters: dndsR records the HyPhy build and seed, but not the CPU type. The earlier finding that RELAX's optimum is not reproducible to 1e-3 in 14 of 300 genes may have the same cause; that was not tested.

## Classification, under the rules fixed in advance

| Error | aBSREL per branch (`p_bh_lineage`) | aBSREL gene `p_value` | Branch models (`p_two_vs_m0`) |
|---|---|---|---|
| E1, E2 corrupted segment | not robust; detectable; masking restores at a cost within the rule → **mask** (warn when flagged) | not robust; detectable; masking restores size but loses more power than the rule allows → **warn** | not robust; detectable; masking restores and gains power → **mask** (warn when flagged) |
| E3, BE3h realigned indels | robust by the preregistered per-test rate (`p_bh_lineage` 5.1% in E3, 0.8% in BE3h), but false discoveries in 11 of 13 lineage families (correction below); `p_raw` not robust; not detectable → **limitation** | not robust, not detectable → **limitation** | robust (inconclusive, not anticonservative) |
| E4 paralog | robust | robust (inconclusive, not anticonservative) | robust |
| E5 discordance | robust; most genes excluded by the existing rule | robust on the 49 tested | robust |
| E6 missing taxa | robust | robust | robust |
| E7 recombination | robust | robust (inconclusive) | robust |

All of these are upstream properties of the models on such inputs (previous section).

### E3 and the homologous-insertion rerun (BE3h)

After the core results, aBSREL on E3's *true* alignments (no realignment; a diagnostic, not a prespecified arm) still gave `p_raw` 11.2% (terminal branches 16.2%, internal 3.1%) and gene `p_value` 36%. That traced to the simulator: `add_indels()` gave each descendant of an inserting branch its own random codons, so an internal insertion was non-homologous in its descendants, an alignment error built into the "true" alignment. The E3 results are kept, and measure that mixed error. The fix is an option, `homologous = TRUE`, with the default unchanged so that stage 2 and E3 reproduce; the change and the rerun are recorded in the design ([lineage-robustness-study.md](../dev/lineage-robustness-study.md)). BE3h reran E3 with it: same seed, methods, criteria and arms, on Intel Skylake (skl-078), about 11 core-hours.

| aBSREL on BE3h | `p_raw` per branch | terminal | internal | gene `p_value` |
|---|---:|---:|---:|---:|
| realigned by `phylo_prepare` (the arm) | 7.1% (6.3-8.0) | 6.8% (5.8-7.9) | 7.6% (6.3-9.1) | 21% (16.5-26.1) |
| true alignment (diagnostic) | 2.3% (1.9-2.8) | 1.8% (1.3-2.4) | 3.1% (2.3-4.1) | 3.0% (1.4-5.6) |
| E0, no error | 2.1% (1.6-2.5) | 1.8% | 2.4% | 2.3% |

On the true alignment BE3h matches E0: gaps themselves do not inflate aBSREL. The excess after realignment is therefore alignment error, spread over terminal and internal branches alike, and HyPhy run directly on dndsR's files reproduces it exactly (previous section). Masking does not remove it, because the window profile flags almost none of these genes (0.3%). `p_bh_lineage`'s per-test rate is 0.8%, but it made false discoveries in 11 of 13 lineage families (correction below); the branch models stay at baseline (`p_two_vs_m0` 5.3%, E0 4.3%). Which aligner step causes the error, and whether another aligner would avoid it, was not tested.

## Correction (2026-11-22, after the maintainer's review): what the BH outputs' rates mean

The size tables above apply the preregistered per-test rule to every output, including `p_bh_lineage` and the BH-adjusted gene `p_adj`. For an output already corrected by Benjamini-Hochberg, that rate is not an error-rate guarantee. When every gene is null, every rejection is false, so the false-discovery rate of a BH family equals the probability that the family has any rejection, which BH holds at or below 0.05. The right check is how many families made any rejection. The preregistered numbers stay as they are; this check is added after the results and is labelled as such.

| Cell | `p_bh_lineage`: families (one per species-tree branch) with any rejection | rejections | aBSREL gene `p_adj` (one family): rejections | branch models' `p_adj_*`: rejections |
|---|---:|---:|---:|---:|
| E0 | 0 of 13 | 0 | 0 | 0 |
| E1 | 1 of 13 (t1) | 600 | | |
| E2 | 1 of 13 (t8) | 300 | | |
| E3 | 11 of 13 | 200 | 120 | 0 |
| BE3h | **11 of 13** | 33 | **30** | 0 |
| E4 | 0 of 13 | 0 | 0 | 0 |
| E5 | 1 of 13 | 1 | 0 | 0 |
| E6 | 0 of 65 | 0 | 0 | 0 |
| E7 | 0 of 13 | 0 | 0 | 0 |

So in BE3h `p_bh_lineage` made false discoveries in 11 of the 13 lineages tested, and aBSREL's gene-level `p_adj` made 30 false discoveries among 300 null genes: neither controls false discoveries under realignment error. The statement in the first version of this report that `p_bh_lineage` "stays within tolerance" and is "the output to use" was wrong. The branch models' adjusted p-values made no rejection in any null cell. The definitions of each p-value, and a mixture study of false-discovery control and power, are in [the multiple-testing design](../dev/lineage-multiple-testing.md); its report will be `lineage-multiple-testing.md`.

## Recommendations, for the maintainer's decision

None of these is implemented; each changes behaviour or output and needs approval.

1. **Warn on the window profile.** dndsR already reports `dissimilar_windows` per orthogroup. A gene with any flagged window should carry a warning in aBSREL and branch-model results, because in the mixture 74-78% of unmasked discoveries were false and the profile flagged 99.5% of corrupted genes with no false flags. This is a warning, not an exclusion: excluding flagged genes failed the power rule for every test.
2. **Offer masking as the documented safeguard** for aBSREL branch tests and the branch models (non-inferior power, FDP 4-5% against 25-74%). It should stay a sensitivity analysis, not a default, unless you decide otherwise. For aBSREL's gene-level `p_value`, masking restores size but costs 22 points of power, so present it with the warning, not as a fix.
3. **Document E3 as a limitation of aBSREL's raw branch and gene p-values after realignment.** *Withdrawn in part, see the correction below: `p_bh_lineage` is not a safe alternative here.* With correctly simulated indels (BE3h) the excess is about 7% per branch and 21% per gene, and comes from realignment, not from gaps: aBSREL on the true alignments is calibrated. The documentation should say that aBSREL's unadjusted branch and gene p-values are sensitive to alignment error in gapped genes, that dndsR's window profile does not detect this. Testing other aligners or alignment filters would be a separate study.
4. **No new safeguard for paralogs, discordance, missing taxa or recombination** on these trees. The existing zero-length-foreground exclusion is doing real work for discordance and should be described as such in the outputs documentation. The tip-dS ratio detects the paralog in every case and could be reported as an informational flag, but the tests stayed calibrated without it.

## Not run (held for the maintainer's decision)

- The 16-taxon units E0, E1 and E3 (about 160 core-hours for aBSREL and branch models).
- Every BUSTED arm: about 790 core-hours as designed, about 630 after the drop order fixed in advance. BUSTED costs 226 core-seconds per gene on 8 taxa and 536 on 16, 8-10 times aBSREL.

## Limits

One simulator; one 8-taxon tree for most cells; 300 codons; one omega distribution. The mixture's power figures are confounded by corruption (above). Gene-tree topology comparison was not evaluated. E5 and E7 use a single nearest-neighbour interchange; deeper discordance was not simulated.
