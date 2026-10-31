# dndsR commands

The list comes from `dndsr --help`. Descriptions come from the README and have not been verified by running the commands. Use `dndsr <command> --help` for arguments.

| Command | Documented purpose |
|---|---|
| `split_comparisons` | Split genomes by subgenome or haplotype and write a new comparison file |
| `extract_cds` | Extract CDS or proteins from FASTA and GFF3 |
| `calculate_dnds` | Calculate dN/dS with `orthologr` |
| `append_annotations` | Add GFF attributes and coordinates to dN/dS results |
| `ipr_enrichment` | InterPro term enrichment under positive selection |
| `go_enrichment` | Gene Ontology enrichment with topGO |
| `term_enrichment` | Enrichment of other terms, such as KEGG or PANTHER |
| `dnds_ideogram` | Ideogram of binned dN/dS across a genome |
| `dnds_summary` | dN/dS distributions and statistics, optionally by region |
| `dnds_contrast` | Contrasts between comparisons |
| `dnds_state_contrast` | Gene-level contrasts |

## Global options

| Option | Meaning |
|---|---|
| `-c`, `--comparison-file` | Comparison file for batch mode |
| `-t`, `--threads` | Workers |
| `-o`, `--output-dir` | Output directory |
| `-v`, `--verbose` | Verbose logging |
| `-w`, `--warnings` | `off`, `summary`, or `all` |

Option names accept `-` or `_` as separators.

## Comparison file

Whitespace-delimited with one comparison per row, documented as `comparison_name query_fasta query_gff subject_fasta subject_gff`. Paths containing spaces must be double-quoted.

```
Cform_v_Calbum ./path_to/Cformosanum.fasta ./path_to/Cformosanum.gff ./path_to/Calbum.fasta ./path_to/Calbum.gff
```

The default `calculate_dnds` output path, from the code, is `<output_dir>/<comparison_name>/<comparison_name>_dnds.tsv`.
