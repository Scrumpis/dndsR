# Migrating to dndsR 0.3

**v0.1.0 and v0.2.0 were unvalidated pre-releases. Results from them should not be reused for dN/dS, enrichment, or summaries.** Please rerun them. Evidence: [real-genome validation](validation/real-genomes.md).

## Which results are affected

| Bug | Affected if | Effect |
|---|---|---|
| CDS assembled in the wrong exon order | Multi-exon genes on the minus strand: nearly every genome | All such genes are wrong. False dN/dS above 1, high dS, undefined ratios, lost orthologs ([public example](validation/public-genome.md)) |
| Isoforms concatenated into chimeras | Genes with more than one annotated isoform | Sequences about twice the expected length, usually with stops |
| Annotation from the wrong genome | A gene ID also appears in a genome processed earlier in the same run, at one thread | Wrong terms and coordinates |
| Empty annotations | RefSeq-style GFF IDs, such as `gene-` and `rna-` | No annotations or coordinates, without an error |
| First BED region lost | Headerless regions file | One region was missing |
| Enrichment | Results from any of the above, plus the fixes listed in [NEWS](../NEWS.md) | Significant terms can change substantially |

Single-exon genes and plus-strand genes with one isoform were assembled correctly, and their dN/dS do not change, though downstream summaries and enrichment that include other genes do.

## How to rerun

1. Update to 0.3 and confirm it: `dndsr verify_example`.
2. Rerun from the start, with a new output folder if you want to keep the old results to compare:

   ```
   dndsr extract_cds -c comparison_file.txt -o results_v0.3 -t 8
   dndsr calculate_dnds -c comparison_file.txt -o results_v0.3 -t 40
   dndsr append_annotations -c comparison_file.txt -o results_v0.3 -t 8
   ```

   In an existing folder, `extract_cds` regenerates old CDS automatically, but `calculate_dnds` keeps existing results and warns you. Add `--overwrite` there.
3. Rerun every later step, such as enrichment, summaries, contrasts, and ideograms.
4. Revisit conclusions drawn from specific genes or terms. Comparing `*_dnds.tsv` between old and new shows which genes changed.

## Interface changes that may affect scripts

- Mistyped options and commands are now errors. A script with an option that was silently ignored will now stop.
- `--warnings` defaults to `summary`.
- Option values stay as typed: `none` is a real value, and `007` stays text for names and paths.
- `extract_cds` leaves invalid CDS out and writes a QC table. Use `--keep-internal-stops true` to include them.
- A failing `--filter-expr` is now an error.
- `go_enrichment`: `enrichment` has a new definition, and there is a new `p_adjust_method` column.
- `ipr_enrichment` and `term_enrichment`: there are new `p_adjust_method` and `n_tested` columns, and `p_adj` is larger, because it corrects for every eligible term. `min_pos` no longer affects `p_adj`.
- The published image's dndsR is installed differently from 0.3.0 on. See [INSTALL.md](../INSTALL.md).
