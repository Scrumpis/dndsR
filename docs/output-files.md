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
