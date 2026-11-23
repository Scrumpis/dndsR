# Lineage-level selection: methods

Four commands compare selection among the branches of a species tree, from the same `phylo_prepare` folder: `phylo_branch_omega` (PAML codeml branch models), `phylo_absrel` (HyPhy aBSREL), `phylo_relax` (HyPhy RELAX), and `phylo_summarize`, which turns their per-gene results into one row per lineage. Tutorial: [tutorials/lineage-selection.md](tutorials/lineage-selection.md). BUSTED, for episodic selection on a fixed foreground, is in [phylo-methods.md](phylo-methods.md).

**Status: experimental, new in this release.** dndsR reproduces each engine run directly ([evidence](validation/phylo-lineage.md)). Their error rates and estimates were measured only in simulation, on five designed trees ([calibration](#calibration-in-simulation)): aBSREL's tests held their level on simulated error-free alignments but not after sequences with indels were realigned ([input errors](#input-errors)), the branch-model tests were close to it, RELAX on 16 taxa was anticonservative, reproducibly, and the calibration of `p_fg_omega_gt1` at a foreground omega of 1 and of `p_free_vs_m0` with missing taxa is not established. Those outputs carry a `calibration` flag in every results table. The methods have run on one small set of public Brassicaceae genes, which shows they work on real data, not that their conclusions are right. They are not validated for biological inference across arbitrary phylogenies.

## Which question each answers

| Command | Question | Unit | What a positive means | What it does not mean |
|---|---|---|---|---|
| `phylo_branch_omega` | How large is dN/dS on each branch, or on the foreground against the background? | Branch (free-ratio) or set of branches (two-ratio), averaged over all sites of a gene | The foreground's gene-wide omega differs from the background's, or exceeds 1 | Positive selection: a gene-wide omega above 1 is rare, and short branches, few synonymous changes, and alignment or orthology errors produce high values |
| `phylo_absrel` | Did some sites on this branch evolve with omega above 1? | Branch | Evidence, under the model, of episodic positive selection at some sites of that branch | Selection across the gene, or adaptation |
| `phylo_relax` | Is selection on test branches weaker (relaxed) or stronger (intensified) than on reference branches? | Gene, test set against reference set | `k < 1`: omega classes pulled towards 1; `k > 1`: pushed away from 1 | Positive selection. Intensification can be stronger purifying selection; relaxation can come from a smaller population size |
| `phylo_summarize` | Across genes, how do lineages compare? | Species-tree branch | A description of the genes tested on that branch | Evidence of adaptation, or of a trait association |

A single number rarely represents selection on a lineage. Gene-wide omega mixes purifying selection at most sites with whatever happens at a few; episodic selection (aBSREL) and relaxation (RELAX) are separate questions with separate models, and they can disagree for good reasons.

## Three kinds of output

Keep them apart when reporting:

| Kind | Outputs | What it is | Uncertainty given |
|---|---|---|---|
| **Model-based estimates**, per gene | `omega_m0`, `omega_bg`, `omega_fg` (codeml); `omega_free` per branch; aBSREL's `baseline_omega`; RELAX's `k` and omega classes | A maximum-likelihood value under a model that is certainly simplified. Not a test | None per gene. Free-ratio omegas are noisy (see calibration) and are for exploration |
| **Tests**, per gene or per branch in a gene | `p_two_vs_m0`, `p_free_vs_m0`, `p_fg_omega_gt1` and their `p_adj_*`; aBSREL `p_raw`, `p_holm`, `p_value`, `p_adj`, `p_bh_lineage`; RELAX `p_value`, `p_adj` | A likelihood-ratio test of one stated null, with a nominal error rate that holds only as far as the calibration shows (some do not hold it: `calibration` column) | p-values, corrected as stated |
| **Descriptive summaries**, per lineage | Everything in `phylo_summarize`'s tables: `pooled_omega` and its interval, medians, counts of significant genes, matched pairs | Summaries over genes of the estimates and tests above. Not tests of lineage hypotheses: genes are not independent, and a lineage difference is not evidence of adaptation or of a trait association | Bootstrap intervals over orthogroups; sign-test p-values described as descriptive |

## Branch identity

All four commands name branches as `phylo_busted` does: by the species-tree branches they correspond to (`species_tree_branches`). Terminal branches are named by taxon; internal branches by the species tree's own node labels if it has them, otherwise `mrca(x,y)`. **Label the internal nodes of your species tree** with plain names (letters, digits, `_`, `.`, `-`): they become the branch identifiers in every table and in the exported tree.

When taxa are missing from an orthogroup, branches of the species tree join; such a branch is named by all of them (`mrca(spA,spD)+spF`), flagged `merged_branch`, and left out of lineage summaries, because its estimate cannot be attributed to one lineage. The two branches at the root always join, because the analyses use unrooted trees: **include an outgroup**, so that no lineage you care about sits next to the root.

## phylo_branch_omega

codeml branch models (Yang 1998), with `NSsites = 0`, F3x4 codon frequencies (`codon_freq`), one transition/transversion ratio, and the tree topology given:

- `M0`: one omega for every branch.
- `free` (`model = 1`): one omega per branch. Exploratory: with one parameter per branch, each estimate is noisy.
- `two_ratio` (`model = 2`): one omega for foreground branches (the foreground file), one for the rest.
- `two_ratio_null`: as `two_ratio`, with the foreground omega fixed at 1.

Tests, each corrected by Benjamini-Hochberg across tested orthogroups: `p_two_vs_m0` (does the foreground differ, 1 df), `p_free_vs_m0` (does omega vary among branches), and `p_fg_omega_gt1` (is the foreground omega above 1; one-sided, from the signed likelihood root). Without a foreground file, only `M0` and `free` run.

Each model is fitted from three starting omegas (0.4, 0.1, 2) and the best likelihood kept; `lnl_spread_*`, `starts_at_best_*`, and the flag `starts_disagree:<model>` show when codeml's optimizer found different optima. A nested model that fits better than its alternative gives a negative likelihood ratio; it is set to 0 and flagged `lrt_negative:<test>`. `codeml_warning:<model>` passes on codeml's own "check convergence" message, which it prints when the optimizer may have stopped short (on the public Brassicaceae genes, 56 of 300 genes); look at those fits before relying on them.

Per-branch flags in `_branches.tsv`, from the free-ratio fit: `zero_length`, `few_synonymous_changes` (fewer than one expected synonymous substitution, `S*dS < 1`, so omega is unstable or undefined), `omega_at_upper_bound` (codeml's 999), `high_ds` (dS above 1, where saturation makes estimates unreliable), `merged_branch`, `not_in_species_tree`. `reliable_free` is true when none of the first four applies. These thresholds are rules of thumb; their value is measured in the calibration below.

**Files.** codeml writes many files. `phylo_branch_omega` keeps, per orthogroup, one alignment, one tree per model, and codeml's control file and main output for each model and start: about 25 files, enough to rerun only what changed. On genome-scale runs, where file systems often limit the number of files, use `--keep-runs false`: each orthogroup's codeml folder is removed once its results are in the tables. HyPhy analyses keep about five files per orthogroup.

## phylo_absrel

aBSREL (Smith et al. 2015) fits a variable number of omega classes per branch and tests each branch for a class above 1. `branches = "all"` (default) tests every branch, the exploratory analysis for comparing lineages; `"foreground"` tests only the foreground, with more power for a hypothesis fixed in advance. `baseline_omega`, from aBSREL's first stage, is one MG94 omega per branch for every branch; it is on a different scale from codeml's GY94 omega and is not interchangeable with it.

Multiple testing. Each output answers a different question, over a different family:

| Output | Question | Family corrected over |
|---|---|---|
| `p_raw` | Was this branch, in this gene, selected? | none |
| `p_holm` | the same, within this gene | the branches tested in the gene. HyPhy multiplies the k-th smallest p by (m − k + 1) without Holm's running maximum, so it lies between Holm's and Hochberg's procedures, valid when the branches' p-values are independent or positively dependent |
| `p_value` (per gene) | Was any tested branch of this gene selected? | the smallest `p_holm`: in effect Hochberg's global test, under the same condition |
| `p_adj` (per gene) | the same, across many genes | Benjamini-Hochberg across tested genes |
| `p_bh_lineage` | In which genes was this particular branch selected? | Benjamini-Hochberg across the genes in which that species-tree branch was tested, one family per branch. Nothing corrects across branches: counting discoveries over several lineages is not covered by it |

Each correction is valid only if `p_raw` is. In a simulated mixture of 1,000 genes (8 taxa, 10% selected on two terminal branches and their stem; [results](validation/lineage-multiple-testing.md)), with error-free alignments, `p_bh_lineage` kept the false-discovery proportion at or below 2.4% within each selected lineage. Lineages without selection made no rejections, and gene `p_adj` kept them at 2.3%. After realigning sequences with indels, about half of the discoveries from every aBSREL output were false, and every lineage without selection made some. The branch models' `p_adj_two_vs_m0` stayed at 1.4%. Which output to use, on correct alignments:

| Question | Output |
|---|---|
| Is this gene under episodic selection on some branch? | gene `p_adj` (`p_value` is not corrected across genes) |
| In which genes is lineage X selected, for a lineage chosen in advance? | `p_bh_lineage` within that lineage |
| Which branches of which genes, with no lineage chosen in advance? | no output guarantees this. Summing `p_bh_lineage` discoveries over lineages held in that simulation only because lineages without selection made no rejections; treat such counts as exploratory |
| One branch of one gene, planned in advance | `p_raw` |

`p_holm` corrects only across the branches of one gene, not across genes. One simulated design is limited evidence: other tree sizes, selection rates, or alignment problems may behave differently ([input errors](#input-errors)).

## phylo_relax

RELAX (Wertheim et al. 2015) fits omega classes on reference branches and the same classes raised to `k` on test branches. The foreground file's rows are the test set; rows with a column `set` = `reference` give the reference set, and without them every other branch is the reference. Both sets must stay identifiable after pruning and may not overlap; otherwise the orthogroup is excluded (`reference_*`, `test_reference_overlap`). `k_at_bound` marks fits at HyPhy's limits for `k` (0 or 50), which are extreme fits rather than measured sizes. `models = "Minimal"` fits only the null and alternative, which is much faster.

## phylo_summarize

One row per species-tree branch in `_lineages.tsv`, plus matched comparisons, foreground and RELAX summaries, and identifiers for joining other data.

- **Pooled omega**, not a mean of ratios: `(sum N*dN / sum N) / (sum S*dS / sum S)` over the genes in which the branch exists, unmerged, with dS at most 1. Genes contribute in proportion to their substitutions, so genes with almost no synonymous changes, whose individual ratios are unstable or undefined, do not dominate. A mean of per-gene ratios is dominated by exactly those. The interval resamples orthogroups. **Its limits, from exploratory reanalysis of the calibration outputs ([audit](dev/validation-audit.md), E1; a confirmatory study is designed):** with one omega at every site it covered the truth about 91% of the time at 25 genes and 94% at 50-100, short of 95% for small gene sets; when omega varies among sites, as it does in real genes, the free-ratio model's pooled omega was 4-11% too low, more so on long branches, so its interval missed the truth more often as genes were added (70% coverage at 100 genes). Treat the interval as a lower bound on the uncertainty, compare lineages of similar branch length or include branch length (`median_branch_length`, `median_ds`) as a covariate in a comparative model, and do not read small differences between lineages as real.
- **Per-gene medians and quartiles** of omega over reliable estimates, beside it: a different summary, of the typical gene rather than of all substitutions.
- **Coverage and diagnostics**: genes tested, estimated, merged or absent, reliable; median dS and branch length; the fractions flagged.
- **Matched pairs** (`_lineage_pairs.tsv`): two lineages estimated from different genes differ partly because the genes differ. Each pair is compared on the genes where both have reliable estimates: pooled omega of each on that shared set, the median per-gene log2 ratio, and a sign test, Benjamini-Hochberg corrected across pairs. Genes share ancestry and function, so these p-values are descriptive.

Everything in this table is descriptive. Formal tests stay per gene, in the analyses summarized.

### Joining lineage traits

`_branch_ids.tsv` has one row per species-tree branch: `branch_id`, `tree_label` (its Newick-safe form), `terminal`, `ape_node` (the node below the branch as `ape` numbers it, so ancestral states from `ape::ace` or `phytools` join directly), `parent_id`, descendant taxa, and the species tree's branch length. `_species_tree_ids.nwk` is the species tree labelled with `tree_label`. Terminal branches are species: join species-level traits (LTR content, genome size) by taxon name. Internal branches join by `branch_id`, for traits reconstructed on ancestors, such as genome-size change along a branch.

Lineage-level hypotheses ("lineages with more LTR activity have higher omega") need a phylogenetic comparative model that accounts for shared ancestry, unequal branch lengths, and the number of genes behind each estimate: for example `nlme::gls` with `ape::corBrownian`, `phylolm`, or `caper::pgls`, with the species tree written here. dndsR does not fit these models. Weighting by `n_pooled` or `sum_syn_changes`, and checking that results do not depend on lineages with few genes or saturated dS, are part of such an analysis.

## Limits

- **Single-copy genes only.** `phylo_prepare` excludes orthogroups with paralogs. Gene families that expand or contract, the ones genome expansion is made of, are therefore under-represented, and more so in lineages with more duplication. Conclusions about genome expansion from these tables apply to the single-copy genes tested.
- **Branch length.** Short branches carry few substitutions, so their omega is noisy; long branches saturate. The flags and the dS columns show both; neither is corrected.
- **Model assumptions.** Branch models assume one omega per branch for all sites; site-to-site variation is ignored (it is measured in the calibration). aBSREL and RELAX assume the tree, the codon model, and correct alignments. Alignment errors inflate positives ([input errors](#input-errors)); multinucleotide substitutions do too, as in the BUSTED calibration ([phylo-methods.md](phylo-methods.md#calibration-in-simulation)).
- **Calibration.** RELAX's `p_value` was anticonservative on 16 taxa (8.3-9.3% at 0.05 in every check); `p_fg_omega_gt1` at omega exactly 1 (5.5-6.8%) and `p_free_vs_m0` with missing taxa (8%, then 4% in a replicate) are not established as calibrated. Rows carry the `calibration` flag; treat single-gene positives from them with caution. Single-gene free-ratio estimates are noisy even when unflagged; compare lineages through pooled values.
- **Gene-tree discordance** is reported per orthogroup (`gene_tree_discordant`), not corrected.
- **Ascertainment.** Only orthogroups that pass QC are analysed; counts by reason are in each run's provenance.

## Calibration in simulation

Design fixed before any simulation: [dev/phylo-lineage-calibration.md](dev/phylo-lineage-calibration.md). Full results: [validation/phylo-lineage-calibration.md](validation/phylo-lineage-calibration.md). Simulated codon data, five trees of 5 to 16 taxa, 300 codons per gene; false-positive rates at 0.05 with 95% intervals:

- **aBSREL**: `p_raw` 1.5-2.6% per branch, gene `p_value` 1.7-4.0%, `p_bh_lineage` under 0.1%. `p_raw` and `p_bh_lineage` met the size criterion on every tree; the gene-level `p_value` on T5 only, and was consistent with 5% elsewhere. Power on a branch with 10% of sites at omega 10: 72-79% by `p_holm`, 93-96% by `p_bh_lineage`.
- **Branch models**: `p_two_vs_m0` 3.3-7.0% and `p_free_vs_m0` 4.0-6.7% with all taxa, none significantly above 5%, also when omega varies among sites; `p_free_vs_m0` 8.0% (5.2-11.7) when taxa are missing from genes, not reproduced in a follow-up (4.0%, 2.1-6.9; 6.0% of the 600), so calibration there is not established. Benjamini-Hochberg-corrected tests: no false positive. `p_fg_omega_gt1` never fired on purifying genes, and at a foreground omega of exactly 1, its boundary, it gave 6.8% (4.9-9.2) on 8 taxa and 5.5% (3.8-7.6) on 16: close to 5%, calibration not established (not tested with missing taxa). Power of `p_two_vs_m0` for foreground 0.6 against 0.2: 98-100%; of `p_fg_omega_gt1` for 1.5: 78%.
- **RELAX** (`Minimal`): 5.0% on 8 taxa, 7.0% on 5, **9.0% (6.0-12.8) on 16 taxa: anticonservative**, reproduced in a follow-up with new genes (8.3%), with sites' classes drawn per branch as RELAX assumes (8.5%), with HyPhy run directly (9.0%), and with the best of three seeds (9.3%). Its likelihood-ratio statistic is larger than chi-square(1) assumes (mean 1.28); dndsR does not correct it. Power: 59-62% for k = 0.5, 35-38% for k = 2, with the direction of `k` right in 95-100% of detections.
- **Estimates**: a single unflagged free-ratio omega is off by a median of 1.3- to 1.4-fold, and more than 2-fold in about one case in five; flagged ones by 8-fold or more, but the flags fire on at most 2.5% of branches. The **pooled omega** is within 3-6% of the truth on a typical branch (at most 15%), better than the median of per-gene values and far better than their mean, which reached 5-fold on a short branch.

These are properties under the simulation model, with correct alignments. Input errors are the next section.

## Input errors

A second simulation study (study B, criteria fixed before running; [report](validation/lineage-robustness-core.md)) added one error at a time to null genes on an 8-taxon tree, and one 16-taxon tree with 4 taxa missing per gene, and on a complete 16-taxon tree, for aBSREL, the branch models (M0, two-ratio) and BUSTED. These findings are for the configurations named, not a general guarantee.

- **A corrupted or misaligned segment in one sequence** gave aBSREL 9.5-9.9% false positives per branch and 99.3-100% per gene (8 and 16 taxa), the two-ratio test 9-51%, and BUSTED 95-96% when the segment was in a foreground sequence. The `dissimilar_window` flag marked every affected gene and none of the error-free ones.
- **Realigning sequences with indels** (`phylo_prepare`'s MAFFT step), with insertions simulated correctly, gave aBSREL 7.1% per branch (`p_raw`) and 21% per gene (`p_value`); on the true alignments the same genes gave 2.3% and 3.0%. So the excess is alignment error. After correction, `p_adj` made 30 false discoveries among 300 null genes, and `p_bh_lineage` made at least one in 11 of 13 lineages. **The `dissimilar_window` flag missed these errors (0.3% of genes).** On 16 taxa the gene `p_value` gave 26% (0.7% on the true alignments), and `p_bh_lineage` made false discoveries in 24 of 29 lineages. The branch models were unaffected (`p_two_vs_m0` 2.7-5.3%), and BUSTED nearly so (6.0%, with 0 and 3 false discoveries after correction on 8 and 16 taxa). HyPhy and codeml run directly on dndsR's files gave the same statistics. A codon-aware aligner (PRANK) removed most of aBSREL's excess in an exploratory pilot, at about 80 times MAFFT's cost; dndsR's aligner is unchanged ([pilot](validation/alignment-pilot.md)).
- **A paralog in place of an ortholog, a gene tree differing from the species tree, missing taxa, and recombination** did not raise either method above its error-free rate in the configurations tested: one paralog split at half the tree's depth, one exchange of two taxa, four of 16 taxa missing, one breakpoint. For discordance and paralogs this held partly because dndsR already excludes genes whose foreground branch has zero length (251 of 300 discordant genes for aBSREL). These configurations do not establish robustness to other forms or degrees of these errors, which are not screened.

**Supplied alignments.** Codon alignments made outside dndsR are accepted with `alignment_dir`. Being in frame is not evidence that the aligned codons are homologous, and the alignment check does not establish it either; the same limits apply as to `phylo_prepare`'s alignments.

**The alignment check.** When `dissimilar_window` flags tested orthogroups, aBSREL, BUSTED, RELAX and the branch models end with a warning giving the number flagged and how many of them are significant, and the run's `_provenance.json` records the same under `alignment_check`. The flag is worth acting on: review those alignments. An unflagged orthogroup is not shown to be correctly aligned.

**Masking** the flagged windows is an optional sensitivity analysis, not a default and not a supported dndsR command (the studies used `tools/dev/validation/phylo-windows.py`). In study B it restored aBSREL's, the branch models' and BUSTED's false-positive rates for corrupted segments, at a power cost of under 1 point for aBSREL branch tests (`p_bh_lineage`), 22 points for aBSREL's gene `p_value`, and 12 points for BUSTED's `p_adj`. It did nothing for realignment error, because the flag does not find those errors. It does not protect against false inference of selection in general: compare masked and unmasked results, and treat a gene significant only unmasked with a flagged window as a likely artifact.

## References

Yang Z (1998) Mol Biol Evol 15:568. Yang Z (2007) Mol Biol Evol 24:1586. Smith MD et al. (2015) Mol Biol Evol 32:1342. Wertheim JO et al. (2015) Mol Biol Evol 32:820.
