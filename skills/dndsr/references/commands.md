# dndsR commands

Every command is also an R function with the same name and options, except `example`, which is `dndsr_example()`. `phylo_busted` takes `--alignment-dir`, `--species-tree-file`, `--foreground-file`, and optionally `--gene-tree-dir`, `--run-name`, `--seed`, `--hyphy-cpu`, `--hyphy`, and HyPhy's settings (`--srv`, `--multiple-hits`, `--rates`, `--syn-rates`, `--error-sink`); its outputs go in `<output_dir>/<run_name>/`. Use `dndsr <command> --help` for arguments.

| Command | Documented purpose |
|---|---|
| `split_comparisons` | Split genomes by subgenome or haplotype and write a new comparison file |
| `extract_cds` | Extract CDS or proteins from FASTA and GFF3 |
| `calculate_dnds` | Calculate dN/dS with `orthologr` |
| `append_annotations` | Add GFF attributes and coordinates to dN/dS results |
| `ipr_enrichment` | InterPro term enrichment under positive selection |
| `go_enrichment` | Gene Ontology enrichment with topGO |
| `term_enrichment` | Enrichment of other terms, such as KEGG or PANTHER |
| `covariate_enrichment` | GO and InterPro enrichment adjusted for gene length and dS. `--pos-threshold` takes a number or a quantile such as `q0.95` |
| `dnds_ideogram` | Ideogram of binned dN/dS across a genome |
| `dnds_summary` | dN/dS distributions and statistics, optionally by region |
| `dnds_contrast` | Contrasts between comparisons |
| `dnds_state_contrast` | Contrasts of selection states between groups; `--level gene` is disabled |
| `phylo_busted` | HyPhy BUSTED tests across species on chosen branches, corrected across orthogroups. Development code, not in v0.4.0; needs HyPhy 2.5.103 |
| `phylo_example` | Copy the bundled phylogenetic example to a folder |
| `verify_phylo_example` | Check the phylogenetic example against HyPhy run directly |
| `example` | Copy the bundled example to a folder |
| `verify_example` | Check example results against independent values |

## Global options

| Option | Meaning |
|---|---|
| `-c`, `--comparison-file` | Comparison file for batch mode |
| `-t`, `--threads` | Workers |
| `-o`, `--output-dir` | Output directory |
| `-v`, `--verbose` | Verbose logging |
| `-w`, `--warnings` | `off`, `summary`, or `all` |

Option and command names accept `-` or `_` as separators. Unknown options are errors.

## Comparison file

Whitespace-delimited with one comparison per row: `comparison_name query_fasta query_gff subject_fasta subject_gff`. Paths containing spaces must be double-quoted.

```
sp1_v_sp2 ./path_to/sp1.fasta ./path_to/sp1.gff ./path_to/sp2.fasta ./path_to/sp2.gff
```

The default `calculate_dnds` output path, from the code, is `<output_dir>/<comparison_name>/<comparison_name>_dnds.tsv`.
