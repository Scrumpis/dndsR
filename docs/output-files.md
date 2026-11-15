# Output files

All outputs for a comparison go in `<output_dir>/<comparison_name>/`. Contrasts go in `<output_dir>/contrast/`. Tables are tab-separated.

## extract_cds

| File | Contents |
|---|---|
| `<genome>_CDS.fasta` | One CDS per gene (or per transcript with `--group-by tx`), named by ID |
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

`<comparison>_dnds_reliability.tsv` (development version) has one row per pair, in the same order. It flags estimates to treat with care and changes no value: the dN/dS table is written exactly as without it.

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

### McDonald-Kreitman (development version)

Written when the comparison has `snp_file` and `snp_side` columns. Method, inputs, and checks: [mcdonald-kreitman.md](mcdonald-kreitman.md).

`<comparison>_mk.tsv` has one row per 1:1 pair that was counted.

| Column | Meaning |
|---|---|
| `query_id`, `subject_id` | The pair |
| `reference_side` | `query` or `subject`: whose SNPs were counted |
| `pn`, `ps` | Polymorphic nonsynonymous (0-fold) and synonymous (4-fold) sites |
| `dn`, `ds` | Divergent nonsynonymous (0-fold) and synonymous (4-fold) differences, path-averaged and so possibly fractional |
| `dn_int`, `ds_int` | `round()` of the above, the counts Fisher's test uses |
| `snps_midfold` | Polymorphisms at 2- and 3-fold sites, excluded from `pn` and `ps` |
| `div_steps_midfold` | Divergence steps at 2- and 3-fold sites, excluded from `dn` and `ds` |
| `codons_used` | Aligned codons counted for divergence |
| `codons_excluded_polymorphic` | Codons left out because they carry a SNP |
| `alpha` | `1 - (ds * pn) / (dn * ps)` |
| `neutrality_index` | `(pn * ds) / (ps * dn)` |
| `dos` | Direction of selection: `pn / (pn + ps) - dn / (dn + ds)` |
| `fisher_p` | Fisher's exact test on `[pn, ps; dn_int, ds_int]`, NA when a margin is zero |
| `fisher_q` | Benjamini-Hochberg across the finite p-values in the file |

`<comparison>_mk_pooled.tsv` has one row: `comparison`, `reference_side`, `n_pairs`, the summed `pn`, `ps`, `dn`, `ds`, `dn_int`, `ds_int`, `snps_midfold`, `div_steps_midfold`, `codons_used`, `codons_excluded_polymorphic`, and `alpha`, `neutrality_index`, `dos`, `fisher_p` recomputed from the sums.

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
| `start_seeds`, `lrt_min_starts`, `lrt_max_starts`, `lnl_unconstrained_best`, `lnl_null_best`, `lrt_best_fits`, `p_best_fits`, `primary_unconstrained_is_best`, `primary_null_is_best` | Diagnostics across starts, with `--extra-seeds` only. They never change `p_value` or `p_adj` |
| `gt_missing_species_splits`, `gt_extra_splits`, `gt_missing_splits`, `gt_discordant_near_foreground`, `gene_tree_mode` | With gene trees: species-tree branches absent from the gene tree, gene-tree branches absent from the species tree, which are missing, and whether one is or touches a foreground branch |
| `terminal_stops_masked`, `max_gap_fraction`, `max_ambiguous_codon_fraction`, `gappy_column_fraction`, `median_p_distance`, `max_p_distance`, `dissimilar_windows` | Alignment QC |
| `qc_flags` | Flags separated by `;`. See [phylo-methods.md](phylo-methods.md#alignment-checks-reported) |
| `alignment_md5`, `json`, `json_md5` | The input alignment's MD5, and HyPhy's JSON with its MD5 |

Branches columns: `orthogroup`, `branch` (HyPhy's name), `terminal`, `taxa_side` (the smaller side of the branch's split of taxa), `foreground`, `species_tree_branches` (species-tree branches joined into this one, `+`-separated; internal ones are named by the tree's node labels or `mrca(x,y)`), `in_species_tree` (whether the branch exists in the pruned species tree), `topology_source` (`species_tree` or `gene_tree`).
