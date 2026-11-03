# dndsR 0.4.0

Nothing here changes an existing default, algorithm, or the content of an existing table.

- `go_enrichment` and `term_enrichment` now warn on every run that they do not adjust for gene length or dS, with the false-term rates measured in simulation. They record this, with their settings, in `<comparison>_go_enrichment_notice.json` or `<comparison>_term_enrichment_notice.json`.
- `term_enrichment` and `ipr_enrichment` now stop if `pos_threshold` is not one number. A text value such as `q0.95` used to be compared with dN/dS as text, so no gene was positive and the run looked like a real result with no terms.
- New `covariate_enrichment()` and `dndsr covariate_enrichment`: GO and InterPro enrichment adjusted for gene length and dS. It is opt-in, and adds output files only. `--pos-threshold` takes a number or a quantile such as `q0.95`, as in `go_enrichment`; the value used is in the `threshold` column. Its simulation evidence does not establish false discovery control in general, and it has not been validated for detecting real selection. Evidence: `docs/validation/enrichment-confounding.md`.
- `calculate_dnds` reports undefined, low-dS, high-dS, and outlier pairs, and says when fewer than 30 genes exceed 1. This is a message only.
- `extract_cds` explains a failure from a trans-spliced gene, such as plastid rps12 in a RefSeq GFF3, and suggests removing organelles.
- `split_comparisons --mode haplotype` warns above 8 labels, which may be scaffold numbers.
- `ipr_enrichment` no longer looks up the InterPro release when entries and tree are local or not used, so it works offline. It runs about 25 times faster, and `append_annotations` about 2.5 times, with byte-identical tables on the public pair.
- New evidence: estimator benchmark, ortholog comparison with OMA, and enrichment calibration by simulation. `go_enrichment`'s excess of false terms is larger than 0.3.0's documents said.

# dndsR 0.3.0

## Rerun earlier analyses

Defects in earlier versions can change results. Rerun analyses made with earlier versions, following `docs/migrating-to-0.3.md`. Between *A. thaliana* and *A. lyrata*, public genomes, every affected gene was wrong, and calls of dN/dS above 1 fell from 1,503 to 206: `docs/validation/public-genome.md`.

- **CDS assembly (all results).** With the default `group_by = "gene"`, multi-exon genes on the minus strand were assembled in the wrong exon order, and genes with several isoforms became chimeras of all isoforms. A designed gene with a true dN/dS near 0.11 was reported as 1.6. See `docs/validation/cds-extraction.md`.
- **Wrong-genome annotations (`append_annotations`).** A gene ID present in both genomes, as with tools that name genes `g1`, `g2`, and so on, took its terms and coordinates from the first genome read.
- **Empty annotations.** RefSeq-style IDs such as `gene-AT1G01010` were not matched, so `append_annotations` and `dnds_ideogram` found no annotations or coordinates. It now also stops if no IDs match.
- **Lost regions.** A headerless regions file lost its first region in `dnds_summary` and `dnds_contrast`.
- **False positives (`term_enrichment`, `ipr_enrichment`).** `min_pos` removed terms before multiple-testing correction, so a term's own outcome decided whether it counted as a test, and `p_adj` was far too small. Correction now spans every term meeting `min_total` (and `max_prop`), and `min_pos` only chooses which rows are written. Expect fewer significant terms. Evidence: `docs/validation/enrichment.md`.
- **Lost annotations (`append_annotations`, `dnds_ideogram`).** A `#` inside a GFF3 attribute, which is legal, cut off the rest of that line, and any terms on it, with no warning. Only a `#` at the start of a line is now a comment. Sequence after `##FASTA` is no longer read as features.

`extract_cds` regenerates old outputs automatically. `calculate_dnds` keeps an existing result but tells you it may be out of date. Add `--overwrite` to recalculate.

## New

- `dndsr example` and `dndsr verify_example` (`dndsr_example()` and `verify_example()` in R) run a designed example and check it against independently established values.
- `extract_cds` writes a QC table, and excludes CDS with in-frame stops or lengths not divisible by 3. `--keep-internal-stops true` keeps them.
- `calculate_dnds` records software versions, parameters, and input hashes, and notices when inputs have changed.
- New documentation: tutorials, output files, interpretation, and validation evidence. A portable Agent Skill is in `skills/dndsr/`.
- `tools/dndsr-install install --bind` makes other folders available in the container.
- Command names accept hyphens: `dndsr extract-cds`.
- `ipr_enrichment` works without internet access, using its bundled InterPro list, and says so.

## Behavior changes

- Failures now stop with an error and a nonzero exit code: a failed parallel worker, a missing input file, an unknown command, or an unknown option. These could previously look like success.
- `--warnings` defaults to `summary`, not `off`.
- Option values are kept as typed: `none` is no longer treated as unset (so `--fdr-method none` works), and names such as `007` and paths containing commas are preserved.
- `filter_expr` that cannot be evaluated is now an error in every command, instead of running unfiltered in some. Rows where it gives NA are removed.
- Each gene's CDS is its longest valid isoform.
- **Default changed: `go_enrichment` now corrects `p_adj` by BH**, once across every tested term in all requested ontologies, for each comparison and side. It was a copy of `raw_p`. New options keep the old or other behavior: `p_adjust = "none"` (`--p-adjust none`), and `adjust_scope = "per_ontology"` (`--adjust-scope per_ontology`). New columns `p_adjust_method`, `adjust_scope`, and `n_tested`. Expect far fewer terms at `p_adj` ≤ 0.05: on a public genome pair, 4 rather than 31. BH does not guarantee nominal error control under all conditions: under a null keeping dN/dS linked to gene length and dS, 9.5% of runs had a false term. Significant terms remain exploratory. Evidence: `docs/validation/enrichment.md`.
- **Default changed: `ipr_enrichment` now corrects once across all eligible InterPro entry types**, including with `method = "parent_child"`, which previously did not support this. `adjust_scope = "per_type"` (`--adjust-scope per_type`) keeps correction within each type. New column `adjust_scope`. On a public genome pair, 15 significant terms rather than 20.
- `go_enrichment` no longer stops with an error when an ontology has no positive genes.
- `go_enrichment`: `enrichment` now has the same definition as the other enrichment commands, the proportion among positives divided by the proportion among other genes.
- `ipr_enrichment` with `method = "parent_child"`: corrected the removal of genes explained by significant child terms, and terms with no parent or child in the tested set are now tested.
- `ipr_enrichment` with stratification: terms of unknown type are no longer tested in every stratum.
- `ipr_enrichment` and `term_enrichment`: new columns `p_adjust_method` and `n_tested`, the number of tests corrected for. It can exceed the number of rows written, so recalculating `p_adj` from the written rows will not match. A type with no positive genes counts toward `n_tested` under `adjust_scope = "global"`.
- `go_enrichment(include_definition = FALSE)` no longer stops with an error.
- `append_annotations`: terms percent-encoded as GFF3 allows, such as `GO%3A0000001`, are now recognized.
- `dnds_ideogram`: dN/dS exactly 1 no longer counts as positive, matching other commands.
- Annotation term patterns no longer match inside other identifiers, and `KEGG:K00001` and `K00001` are the same term.
- `dnds_state_contrast --level gene` is disabled and stops with an error, because it gave no meaningful result: a single table returned an odds ratio of 1 and p = 1 whatever the data.

## Dependencies

- The bundled Arial font is replaced with Liberation Sans Bold (SIL Open Font License). Plots look nearly the same.
- `tools/environment.yml` gives a tested conda installation.
- `orthologr` moved to `Suggests`, and `calculate_dnds` explains how to get it if missing. The orthologr patch is no longer needed.
- `tools/Dockerfile` pins its base image and tools, and installs dndsR into the image. It could not previously be rebuilt, because upstream orthologr had absorbed the patch.
