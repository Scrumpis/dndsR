# Lineage-level selection: methods

Four commands compare selection among the branches of a species tree, from the same `phylo_prepare` folder: `phylo_branch_omega` (PAML codeml branch models), `phylo_absrel` (HyPhy aBSREL), `phylo_relax` (HyPhy RELAX), and `phylo_summarize`, which turns their per-gene results into one row per lineage. Tutorial: [tutorials/lineage-selection.md](tutorials/lineage-selection.md). BUSTED, for episodic selection on a fixed foreground, is in [phylo-methods.md](phylo-methods.md).

**Status: experimental, new in this release.** dndsR reproduces each engine run directly ([evidence](validation/phylo-lineage.md)). Their error rates and estimates with dndsR's inputs were measured only in simulation, on a few designed trees ([calibration](#calibration-in-simulation)), and the methods have run on one small set of public Brassicaceae genes, which shows they work on real data, not that their conclusions are right. They are not validated for biological inference across arbitrary phylogenies.

## Which question each answers

| Command | Question | Unit | What a positive means | What it does not mean |
|---|---|---|---|---|
| `phylo_branch_omega` | How large is dN/dS on each branch, or on the foreground against the background? | Branch (free-ratio) or set of branches (two-ratio), averaged over all sites of a gene | The foreground's gene-wide omega differs from the background's, or exceeds 1 | Positive selection: a gene-wide omega above 1 is rare, and short branches, few synonymous changes, and alignment or orthology errors produce high values |
| `phylo_absrel` | Did some sites on this branch evolve with omega above 1? | Branch | Evidence, under the model, of episodic positive selection at some sites of that branch | Selection across the gene, or adaptation |
| `phylo_relax` | Is selection on test branches weaker (relaxed) or stronger (intensified) than on reference branches? | Gene, test set against reference set | `k < 1`: omega classes pulled towards 1; `k > 1`: pushed away from 1 | Positive selection. Intensification can be stronger purifying selection; relaxation can come from a smaller population size |
| `phylo_summarize` | Across genes, how do lineages compare? | Species-tree branch | A description of the genes tested on that branch | Evidence of adaptation, or of a trait association |

A single number rarely represents selection on a lineage. Gene-wide omega mixes purifying selection at most sites with whatever happens at a few; episodic selection (aBSREL) and relaxation (RELAX) are separate questions with separate models, and they can disagree for good reasons.

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

Multiple testing, three levels: `p_raw` (per branch), `p_holm` (HyPhy's Holm-Bonferroni across the branches tested in one gene), and per gene `p_value` (the smallest `p_holm`, a valid test that some branch was selected) with `p_adj` (Benjamini-Hochberg across genes). For lineages, `p_bh_lineage` applies Benjamini-Hochberg to `p_raw` across the orthogroups in which that species-tree branch was tested: the count of genes with `p_bh_lineage <= 0.05` on a branch is the lineage-level count to report.

## phylo_relax

RELAX (Wertheim et al. 2015) fits omega classes on reference branches and the same classes raised to `k` on test branches. The foreground file's rows are the test set; rows with a column `set` = `reference` give the reference set, and without them every other branch is the reference. Both sets must stay identifiable after pruning and may not overlap; otherwise the orthogroup is excluded (`reference_*`, `test_reference_overlap`). `k_at_bound` marks fits at HyPhy's limits for `k` (0 or 50), which are extreme fits rather than measured sizes. `models = "Minimal"` fits only the null and alternative, which is much faster.

## phylo_summarize

One row per species-tree branch in `_lineages.tsv`, plus matched comparisons, foreground and RELAX summaries, and identifiers for joining other data.

- **Pooled omega**, not a mean of ratios: `(sum N*dN / sum N) / (sum S*dS / sum S)` over the genes in which the branch exists, unmerged, with dS at most 1. Genes contribute in proportion to their substitutions, so genes with almost no synonymous changes, whose individual ratios are unstable or undefined, do not dominate. A mean of per-gene ratios is dominated by exactly those. The interval resamples orthogroups.
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
- **Model assumptions.** Branch models assume one omega per branch for all sites; site-to-site variation is ignored (it is measured in the calibration). aBSREL and RELAX assume the tree, the codon model, and correct alignments. Multinucleotide substitutions, alignment errors, and misassigned orthologs inflate positives, as in the BUSTED calibration ([phylo-methods.md](phylo-methods.md#calibration-in-simulation)).
- **Gene-tree discordance** is reported per orthogroup (`gene_tree_discordant`), not corrected.
- **Ascertainment.** Only orthogroups that pass QC are analysed; counts by reason are in each run's provenance.

## Calibration in simulation

Design fixed before any simulation: [dev/phylo-lineage-calibration.md](dev/phylo-lineage-calibration.md). Results: [validation/phylo-lineage-calibration.md](validation/phylo-lineage-calibration.md).

## References

Yang Z (1998) Mol Biol Evol 15:568. Yang Z (2007) Mol Biol Evol 24:1586. Smith MD et al. (2015) Mol Biol Evol 32:1342. Wertheim JO et al. (2015) Mol Biol Evol 32:820.
