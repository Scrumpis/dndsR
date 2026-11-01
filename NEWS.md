# dndsR 0.3.0

## Rerun earlier analyses

Defects in earlier versions can change results. Rerun analyses made with earlier versions, following `docs/migrating-to-0.3.md`. On real genomes, every affected gene was wrong, about half of pairs changed, and most dN/dS calls above 1 were artifacts: `docs/validation/real-genomes.md`.

- **CDS assembly (all results).** With the default `group_by = "gene"`, multi-exon genes on the minus strand were assembled in the wrong exon order, and genes with several isoforms became chimeras of all isoforms. A designed gene with a true dN/dS near 0.11 was reported as 1.6. See `docs/validation/cds-extraction.md`.
- **Wrong-genome annotations (`append_annotations`).** A gene ID present in both genomes, as with tools that name genes `g1`, `g2`, and so on, took its terms and coordinates from the first genome read.
- **Empty annotations.** RefSeq-style IDs such as `gene-AT1G01010` were not matched, so `append_annotations` and `dnds_ideogram` found no annotations or coordinates. It now also stops if no IDs match.
- **Lost regions.** A headerless regions file lost its first region in `dnds_summary` and `dnds_contrast`.

`extract_cds` regenerates old outputs automatically. `calculate_dnds` keeps an existing result but tells you it may be out of date. Add `--overwrite` to recalculate.

## New

- `dndsr example` and `dndsr verify_example` (`dndsr_example()` and `verify_example()` in R) run a designed example and check it against independently established values.
- `extract_cds` writes a QC table, and excludes CDS with in-frame stops or lengths not divisible by 3. `--keep-internal-stops true` keeps them.
- `calculate_dnds` records software versions, parameters, and input hashes, and notices when inputs have changed.
- New documentation: tutorials, output files, interpretation, and validation evidence. A portable Agent Skill is in `skills/dndsr/`.
- `tools/dndsr-install install --bind` makes other folders available in the container.
- Command names accept hyphens: `dndsr extract-cds`.

## Behavior changes

- Failures now stop with an error and a nonzero exit code: a failed parallel worker, a missing input file, an unknown command, or an unknown option. These could previously look like success.
- `--warnings` defaults to `summary`, not `off`.
- Option values are kept as typed: `none` is no longer treated as unset (so `--fdr-method none` works), and names such as `007` and paths containing commas are preserved.
- `filter_expr` that cannot be evaluated is now an error in every command, instead of running unfiltered in some. Rows where it gives NA are removed.
- Each gene's CDS is its longest valid isoform.
- `go_enrichment`: `enrichment` now has the same definition as the other enrichment commands, the proportion among positives divided by the proportion among other genes. A new column, `p_adjust_method`, states whether `p_adj` was adjusted.
- `ipr_enrichment` with `method = "parent_child"`: corrected the removal of genes explained by significant child terms, and terms with no parent or child in the tested set are now tested.
- `ipr_enrichment` with stratification: terms of unknown type are no longer tested in every stratum.
- `dnds_ideogram`: dN/dS exactly 1 no longer counts as positive, matching other commands.
- Annotation term patterns no longer match inside other identifiers, and `KEGG:K00001` and `K00001` are the same term.
- `dnds_state_contrast --level gene` is marked experimental.

## Dependencies

- The bundled Arial font is replaced with Liberation Sans Bold (SIL Open Font License). Plots look nearly the same.
- `tools/environment.yml` gives a tested conda installation.
- `orthologr` moved to `Suggests`, and `calculate_dnds` explains how to get it if missing. The orthologr patch is no longer needed.
- `tools/Dockerfile` pins its base image and tools, and installs dndsR into the image. It could not previously be rebuilt, because upstream orthologr had absorbed the patch.
