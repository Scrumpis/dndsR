# Output files

All outputs for a comparison go in `<output_dir>/<comparison_name>/`. Contrasts go in `<output_dir>/contrast/`. Tables are tab-separated.

## extract_cds

| File | Contents |
|---|---|
| `<genome>_CDS.fasta` | One CDS per gene (or per transcript with `--group-by tx`), named by its GFF3 `ID` without a type prefix such as `gene-`, `rna-`, or `gene:`, whatever its `Name` (from v0.5.0; before, a file whose gene symbols were all unique was named by symbol) |
| `<genome>_AA.fasta` | Translated proteins, with `--export-proteins` |
| `<genome>_CDS_qc.tsv` | One row per gene or transcript, including those not written |
| `*_extract_cds.log` | Log for each genome |

QC columns: `id`, `transcript` (the isoform used), `n_isoforms`, `length_nt`, `multiple_of_3`, `internal_stop` (in-frame stop; NA if the length is not a multiple of 3), `written`.

## calculate_dnds

`<comparison>_dnds.tsv` has one row per orthologous pair.

| Column | Meaning |
|---|---|
| `query_id`, `subject_id` | Reciprocal best hit pair |
| `dN` | Nonsynonymous substitutions per nonsynonymous site |
| `dS` | Synonymous substitutions per synonymous site |
| `dNdS` | dN divided by dS. NA where dS is 0 or the estimate is unavailable |
| `perc_identity` to `score_raw` | DIAMOND statistics for the protein alignment, including `q_len`, `s_len`, `qcov`, and `evalue` |

`<comparison>_dnds_provenance.json` records versions, parameters, and MD5 hashes of the CDS inputs.

`<comparison>_dnds_reliability.tsv` (v0.5.0 and later) has one row per pair, in the same order. It flags estimates to treat with care and changes no value: the dN/dS table is written exactly as without it.

| Column | Meaning |
|---|---|
| `query_id`, `subject_id` | As in the dN/dS table |
| `codons_compared` | Codons with A, C, G, T in both sequences and a stop in neither |
| `syn_sites`, `syn_differences`, `nonsyn_differences` | Nei-Gojobori (1986) counts on those codons, whatever `--dnds-method` was |
| `flags` | Semicolon-separated, empty when none |

| Flag | Criterion |
|---|---|
| `few_codons` | Fewer than 100 codons compared |
| `few_synonymous_differences` | Fewer than 3 synonymous differences, so dS rests on very few events |
| `ds_below_observed` | dS below half the synonymous p-distance, with at least one synonymous difference: a failed estimate or a large difference between methods |
| `extreme_ratio` | dN/dS above 10 |
| `low_ds`, `high_ds`, `outlier_ds` | dS below 0.02, above 0.8, or above 3 times the table's median, as counted in the `calculate_dnds` message |
| `no_alignment` | No codon alignment was found for the pair |

Evidence for the criteria: [validation/estimators.md](validation/estimators.md), "Reliability flags".

### McDonald-Kreitman (experimental, v0.5.0 and later)

Written when the comparison has `snp_file` and `snp_side` columns. Method, inputs, and checks: [mcdonald-kreitman.md](mcdonald-kreitman.md).

`<comparison>_mk.tsv` has one row per 1:1 pair that was counted.

| Column | Meaning |
|---|---|
| `query_id`, `subject_id` | The pair |
| `reference_side` | `query` or `subject`: whose SNPs were counted |
| `pn`, `ps` | Polymorphic nonsynonymous and synonymous sites, by the amino-acid effect of the ALT allele |
| `dn`, `ds` | Divergent nonsynonymous and synonymous differences, path-averaged and so possibly fractional |
| `dn_int`, `ds_int` | `round()` of the above, the counts Fisher's test uses |
| `codons_used` | Aligned codons counted for divergence |
| `codons_excluded_polymorphic` | Codons left out because they carry a SNP |
| `alpha` | `1 - (ds * pn) / (dn * ps)` |
| `neutrality_index` | `(pn * ds) / (ps * dn)` |
| `dos` | Direction of selection: `pn / (pn + ps) - dn / (dn + ds)` |
| `fisher_p` | Fisher's exact test on `[pn, ps; dn_int, ds_int]`, NA when a margin is zero |
| `fisher_q` | Benjamini-Hochberg across the finite p-values in the file |

`<comparison>_mk_pooled.tsv` has one row: `comparison`, `reference_side`, `n_pairs`, the summed `pn`, `ps`, `dn`, `ds`, `dn_int`, `ds_int`, `codons_used`, `codons_excluded_polymorphic`, and `alpha`, `neutrality_index`, `dos`, `fisher_p` recomputed from the sums.

## append_annotations

`<comparison>_dnds_annot.tsv` adds to the above: `q_attributes`, `s_attributes`, `q_<term>`, and `s_<term>` for terms found (`ipr`, `go`, `kegg`, `panther`, `pfam`, `tigrfam`, `cog`, plus any `--custom` patterns), and `q_gff_seqname`, `q_gff_start`, `q_gff_end` with their `s_` equivalents. `q_` means query and `s_` means subject.

## Enrichment

| Command | Table |
|---|---|
| `ipr_enrichment` | `<comparison>_<q or s>_IPR_<ENTRY_TYPE>_enrichment.tsv`, one per entry type |
| `term_enrichment` | `<comparison>_<q or s>_<TERM>_enrichment.tsv` |
| `go_enrichment` | `<comparison>_<q or s>_GO_<ontology>_topGO_<algorithm>_<statistic>.tsv` |
| `covariate_enrichment` | `<comparison>_<q or s>_<GO or IPR>_covariate_enrichment.tsv` |

Each command also writes a plot, and a log for each side. `go_enrichment` and `term_enrichment` also write `<comparison>_go_enrichment_notice.json` or `<comparison>_term_enrichment_notice.json`, with the calibration notice, the dndsR version, and the settings. Their tables do not change.

The tested groups: positive genes have `dNdS > pos_threshold` (default 1), and the rest are other genes passing the filter (`dNdS` not NA and below `max_dnds`, default 10).

- `ipr_enrichment` and `term_enrichment`: `pos_count`, `nonpos_count`, `pos_total`, `nonpos_total`, `total_count`, `odds_ratio`, `ci_lower`, `ci_upper`, `p_value` (one-sided Fisher's exact test, not corrected for multiple testing), `p_adj` (corrected by `p_adjust_method`), `p_adjust_method`, `n_tested`, `enrichment`, and in `ipr_enrichment` `adjust_scope`: `global` (one family across all eligible entry types, the default), `per_type`, `pooled` (with `stratify_by_type = FALSE`), or `none`.
- Correction spans every term meeting `min_total` and, for `ipr_enrichment`, `max_prop`. `n_tested` is that number. `min_pos` then chooses which rows are written, so a table usually has fewer rows than `n_tested`, and `p_adj` cannot be recalculated from the rows alone.
- `covariate_enrichment`: `term`, `term_type`, `ontology`, `term_name`, `threshold` (the dN/dS above which a gene is positive), `universe` and `pos_total` (annotated and positive genes in that universe), `total_count`, `pos_count`, `expected` (positive genes predicted from length and dS), `obs_over_expected`, `fisher_p` (without adjustment, for comparison, and not corrected), `p_value` (adjusted, not corrected), `p_adj` (BH across `n_tested` terms, in all ontologies of that table), `n_tested`, and `model`. `min_pos` chooses which rows are written, after correction.
- `go_enrichment`: `annotated`, `significant`, `raw_p` (topGO's p-value for the chosen `algorithm`, not corrected for multiple testing), `p_adj`, `p_adjust_method`, `adjust_scope`, `n_tested`, `enrichment`. By default `p_adj` is BH across every tested term in all requested ontologies, for one comparison and side, and `n_tested` counts them, so it exceeds a table's rows. With `adjust_scope` `per_ontology`, each ontology has its own family. With `p_adjust_method` `none`, **`p_adj` equals `raw_p`**, `adjust_scope` is `none`, and `n_tested` is still the family's size.
- `enrichment` everywhere is the proportion of positive genes with the term divided by the proportion of other genes with it.

## dnds_ideogram

`*_ideogram_input.tsv` (input with coordinates), and one SVG and PNG per side.

## dnds_summary

`<comparison>_dnds_summary.tsv` has `scope` (`global` or `regional`), `side`, `region_status` (`all`, `region`, or `background`), counts of rows removed (`n_removed_na`, `n_removed_ge_max`), `n`, `mean_dnds`, `median_dnds`, `sd_dnds`, `se_dnds`, `ci_lower`, `ci_upper`, and `p_value`. The `normal` CI uses mean ± 1.96 SE, which is too narrow for small `n`, so prefer `--ci-method bootstrap` there.

## dnds_contrast

`contrast/<contrast>_<mode>_<sides>_dnds_contrast.tsv` has `n` paired genes, `mean_dnds_A`, `mean_dnds_B`, `mean_delta` (A minus B), `delta_ci_lower`, `delta_ci_upper`, `wilcox_p_value` (paired signed-rank), `n_pos_A`, `n_pos_B`, `discordant_Apos_Bnon`, `discordant_Anon_Bpos`, `mcnemar_p_two_sided`, and `binom_p_A_gt_B_on_discordant`. No multiple-testing adjustment is applied across contrasts.

## dnds_state_contrast

`--level gene` is disabled, because its results were not meaningful. The default, `--level pairwise`, classifies each row as Positive, Neutral, or Purifying and compares groups given by a `group` column in `--dnds-table-file`. `<label>_dnds_state_contrasts.tsv` has `state`, `groupA`, `groupB`, `model`, `log_or`, `se`, `p_value`, `p_adj`, `odds_ratio` (A relative to B), `or_ci_lower`, and `or_ci_upper`.

## phylo_prepare

Experimental. Outputs go in `<output_dir>/<run_name>/` (default `prepared`). Methods: [phylo-methods.md](phylo-methods.md).

| File | Contents |
|---|---|
| `alignments/<orthogroup>.fasta` | Codon alignment of each prepared orthogroup, sequences named `taxon|gene_id` |
| `proteins/<orthogroup>.faa` | The MAFFT protein alignment it was built from |
| `trees/<orthogroup>.nwk` | The pruned, unrooted topology with `{FG}` labels, as HyPhy will receive it |
| `gene_trees/<orthogroup>.nwk` | Supplied gene trees, with tips renamed to taxa (when given) |
| `cds/<taxon>/` | CDS extracted by `extract_cds` (genome input) |
| `species_tree.nwk`, `foreground.tsv` | Copies of the inputs, which `phylo_busted` checks against the provenance |
| `<run>_summary.tsv` | One row per orthogroup: `prepared`, `excluded`, or `failed` |
| `<run>_members.tsv` | One row per gene of every orthogroup: `copies` of its taxon, `cds_status`, and whether it was `used` |
| `<run>_sequences.tsv` | Per sequence of each orthogroup that passed membership and CDS checks: codons, `terminal_stop_removed`, `starts_with_atg`, `ambiguous_codons` |
| `<run>_provenance.json` | dndsR, R, ape, Biostrings, and MAFFT versions and arguments; settings; inputs with MD5s; counts by status, reason, and basis; prepared orthogroups per taxon |
| `prepare_manifest.json` | Where `phylo_busted --prepared-dir` finds everything |

Summary columns: `orthogroup`, `status`, `reason`, `exclusion_basis`, `detail`, `n_taxa`, `taxa`, `taxa_absent`, `n_codons`, `terminal_stops_removed`, `fg_requested_branches` and `fg_absent_branches` (species-tree branches requested, and those with no taxa present), the alignment QC metrics listed below, `n_fg_branches`, `qc_flags`, `alignment_md5`, `tree_md5`; and with gene trees, the `gt_*` columns below.

## phylo_busted

Phylogenetic results are separate from pairwise ones. They go in `<output_dir>/<run_name>/` (default `busted`), and contain no dN/dS estimates.

| File | Contents |
|---|---|
| `<run>_results.tsv` | One row per orthogroup, tested, excluded, or failed; with `--prepared-dir`, preparation's exclusions too |
| `<run>_branches.tsv` | Every branch given to HyPhy, for each orthogroup that reached HyPhy |
| `<run>_sequences.tsv` | Per-sequence QC |
| `<run>_provenance.json` | dndsR, R, ape, and HyPhy versions; HyPhy path; settings; seed, extra seeds, and HyPhy threads; inputs and MD5s; the prepared folder and its provenance's MD5; the foreground hypothesis; the correction family; counts by status, reason, and basis; tested orthogroups per taxon |
| `hyphy/<orthogroup>/` | The exact `alignment.fas` and `tree.nwk` given to HyPhy, `hyphy.log`, `busted.json`, and `record.json` (what a rerun must match to be reused); extra starts in `start_seed<k>/` |

Results columns:

| Column | Meaning |
|---|---|
| `orthogroup`, `method` | Orthogroup; `busted` |
| `status`, `stage`, `reason`, `exclusion_basis`, `detail` | `tested`, `excluded`, or `failed`; in `prepare` or `analysis`; why; `mathematical`, `no_data`, `policy`, `data_quality`, `input`, or `software`; a readable explanation |
| `n_taxa`, `n_codons`, `taxa_absent` | Size of the alignment, and species-tree taxa not in it |
| `n_fg_branches`, `fg_branches` | Foreground branches given to HyPhy, by its names |
| `fg_requested_branches`, `fg_absent_branches` | Species-tree branches requested, and those with no taxa present |
| `lrt`, `p_value` | HyPhy's likelihood ratio and p-value, from `seed`. p is 0.5 when the ratio is 0 |
| `p_adj`, `p_adj_method`, `family_size` | Benjamini-Hochberg across tested orthogroups in this run |
| `constrained_fitted`, `lnl_unconstrained`, `lnl_constrained` | Log-likelihoods. HyPhy fits the constrained model only when the unconstrained one shows evidence |
| `fg_omega_<k>`, `fg_prop_<k>`, `bg_omega_<k>`, `bg_prop_<k>` | ω and weight of HyPhy's rate class k (its class k-1) on foreground and background branches, unconstrained model. Classes are not sorted |
| `tree_length` | Sum of branch lengths in the unconstrained model |
| `n_starts` | Starts run: 1, or more with `--extra-seeds` |
| `start_seeds`, `lrt_min_starts`, `lrt_max_starts`, `lnl_unconstrained_best`, `lnl_null_best`, `lrt_best_fits`, `p_best_fits`, `primary_unconstrained_is_best`, `primary_null_is_best` | Across starts, with `--extra-seeds` only. `lrt_best_fits` and `p_best_fits` compare the best alternative and best null fits. They never change `p_value` or `p_adj` |
| `p_best_fits_adj` | Benjamini-Hochberg over `p_best_fits` across the tested orthogroups (the `p_adj` family), with `--extra-seeds` only. The value to report under the recommended three starts (v0.5.0 and later) |
| `gt_missing_species_splits`, `gt_extra_splits`, `gt_missing_splits`, `gt_discordant_near_foreground`, `gene_tree_mode` | With gene trees: species-tree branches absent from the gene tree, gene-tree branches absent from the species tree, which are missing, and whether one is or touches a foreground branch |
| `terminal_stops_masked`, `max_gap_fraction`, `max_ambiguous_codon_fraction`, `gappy_column_fraction`, `median_p_distance`, `max_p_distance`, `dissimilar_windows` | Alignment QC |
| `qc_flags` | Flags separated by `;`. See [phylo-methods.md](phylo-methods.md#alignment-checks-reported) |
| `alignment_md5`, `json`, `json_md5` | The input alignment's MD5, and HyPhy's JSON with its MD5 |

Branches columns: `orthogroup`, `branch` (HyPhy's name), `terminal`, `taxa_side` (the smaller side of the branch's split of taxa), `foreground`, `species_tree_branches` (species-tree branches joined into this one, `+`-separated; internal ones are named by the tree's node labels or `mrca(x,y)`), `in_species_tree` (whether the branch exists in the pruned species tree), `topology_source` (`species_tree` or `gene_tree`).

## phylo_branch_omega

Experimental, development code. In `<output_dir>/<run_name>/` (default `branch_omega`). Methods and flags: [lineage-selection.md](lineage-selection.md).

`_results.tsv`, one row per orthogroup:

| Column | Meaning |
|---|---|
| `orthogroup`, `status`, `stage`, `reason`, `exclusion_basis`, `detail`, `n_taxa`, `n_codons`, `taxa_absent`, `qc_flags` | As for `phylo_busted`; flags add `starts_disagree:<model>`, `codeml_warning:<model>`, `lrt_negative:<test>`, `<model>_failed` |
| `lnl_m0`, `np_m0`, `kappa_m0`, `omega_m0`, `tree_length_m0` | One-ratio model: log-likelihood, parameters, transition/transversion ratio, omega, tree length |
| `lnl_free`, `np_free`, `lrt_free_vs_m0`, `df_free_vs_m0`, `p_free_vs_m0`, `p_adj_free_vs_m0` | Free-ratio model and its test against M0, BH across tested orthogroups |
| `two_ratio_status`, `two_ratio_reason`, `n_fg_branches` | Whether the two-ratio model ran (`tested`, `excluded` with a reason when the foreground is not identifiable, `not_run` without a foreground) |
| `lnl_two_ratio`, `omega_bg`, `omega_fg`, `lrt_two_vs_m0`, `p_two_vs_m0`, `p_adj_two_vs_m0` | Two-ratio model and its test against M0 |
| `lnl_two_ratio_null`, `lrt_fg_omega_eq1`, `p_fg_omega_gt1`, `p_adj_fg_omega_gt1` | Foreground omega fixed at 1, and the one-sided test of foreground omega above 1 |
| `lnl_spread_<model>`, `starts_at_best_<model>`, `starts_ok_<model>` | Optimizer diagnostics across starting omegas |
| `n_branches`, `n_branches_reliable` | Branches, and those with reliable free-ratio estimates |
| `calibration` | On tested rows, the tests in this table whose calibration in simulation is a known problem, as `<output>:<code>` separated by `;`: `p_free_vs_m0:not_established_in_simulation`, `p_fg_omega_gt1:not_established_in_simulation`. Present only when those tests were computed. Explanation in the run's warning and `_provenance.json` (`calibration`), evidence in [validation/phylo-lineage-calibration.md](validation/phylo-lineage-calibration.md) |

`_branches.tsv`, one row per orthogroup and branch: `branch`, `terminal`, `taxa_side`, `species_tree_branches`, `in_species_tree`, `topology_source`, `foreground`; for each fitted model (`_m0`, `_free`, `_two`): `t` (branch length), `N`, `S` (sites), `omega`, `dN`, `dS`, `NdN`, `SdS` (expected substitutions); `flags` and `reliable_free`. `paml/<orthogroup>/<model>/start<k>/` keeps codeml's inputs and outputs.

## phylo_absrel

Experimental, development code. `_results.tsv`, one row per orthogroup: as for `phylo_busted`, with `n_branches_tested`, `n_branches_holm_005`, `p_value` (the smallest `p_holm`), `p_adj` (BH across orthogroups), `branch_min_p`, `lnl_baseline`, `lnl_adaptive`. `_branches.tsv`, one row per orthogroup and branch: the branch identity columns, `tested_by_hyphy`, `baseline_omega` and `baseline_length` (aBSREL's one-omega-per-branch MG94 fit), `adaptive_length`, `adaptive_dn`, `adaptive_ds`, `rate_classes`, `omega_max` and `omega_max_prop` (the largest omega class and its weight), `lrt`, `p_raw`, `p_holm` (HyPhy's correction across the branches tested in the gene, Holm's step without the running maximum), `p_bh_lineage` and `lineage_family_size` (BH per species-tree branch across orthogroups; one family per branch, nothing corrected across branches: [lineage methods](lineage-selection.md#phylo_absrel)). For every phylogenetic command, `_provenance.json` has `alignment_check` when orthogroups were tested: counts of tested and `dissimilar_window`-flagged orthogroups, flagged ones with an adjusted p at most 0.05, and what the flag does and does not show.

## phylo_relax

Experimental, development code. `_results.tsv`, one row per orthogroup: as for `phylo_busted`, with `k`, `direction` (`relaxation` or `intensification`), `k_at_bound`, `lrt`, `p_value`, `p_adj`, `lnl_alternative`, `lnl_null`, `mg94_omega_test`, `mg94_omega_reference` (one omega per set), and the alternative model's classes (`test_omega_<k>`, `test_prop_<k>`, `ref_omega_<k>`, `ref_prop_<k>`), and `calibration` (`p_value:anticonservative_in_simulation` on tested rows: RELAX's p-values gave too many false positives in simulation; see `phylo_branch_omega` above). `_branches.tsv` adds `set` (`test`, `reference`, or `unclassified`).

## phylo_summarize

Experimental, development code. In `<output_dir>/<run_name>/` (default `summary`):

| File | Contents |
|---|---|
| `_branch_ids.tsv` | One row per species-tree branch: `branch_id`, `tree_label`, `terminal`, `ape_node`, `parent_id`, `descendant_taxa`, `n_descendants`, `species_tree_length` |
| `_species_tree_ids.nwk` | The species tree with internal nodes labelled by `tree_label` |
| `_lineages.tsv` | One row per branch: `n_orthogroups_tested`, `n_estimated`, `n_merged_or_absent` (split into `n_merged`, in a gene where pruning joined the branch to a neighbour, and `n_absent`, where its taxa are missing so the branch does not exist), `n_reliable`, `n_pooled`, `pooled_omega` with `pooled_omega_lo` and `pooled_omega_hi` (bootstrap over orthogroups), `median_omega_reliable`, `q25_omega_reliable`, `q75_omega_reliable`, `frac_reliable_omega_gt1`, `median_ds`, `median_branch_length`, `sum_syn_changes`, `frac_high_ds`, `frac_few_synonymous`, `frac_omega_upper_bound`, `frac_zero_length`, `frac_starts_disagree`; with an aBSREL run, `absrel_n_tested`, `absrel_n_sig_lineage`, `absrel_n_sig_holm`, `absrel_median_baseline_omega` |
| `_lineage_pairs.tsv` | Matched comparisons on shared reliable orthogroups: `branch_a`, `branch_b`, `n_shared`, `pooled_omega_a`, `pooled_omega_b`, `median_log2_ratio`, `n_a_higher`, `n_b_higher`, `sign_test_p`, `sign_test_p_adj` |
| `_orthogroup_branches.tsv` | One row per orthogroup and branch: the free-ratio estimates with `branch_id` (empty for merged branches), flags, and `reliable` |
| `_foreground.tsv` | Two-ratio summary: genes tested, significant with the foreground higher or lower, `p_fg_omega_gt1` count, median log2 foreground/background omega, pooled foreground and background omega |
| `_relax.tsv` | RELAX summary: genes tested, significant relaxation and intensification, median log2 `k` |
| `_provenance.json` | Inputs and their checksums, settings, and the definitions used |
