# Analyzing your own genomes

This assumes you have completed [getting started](getting-started.md). The steps are the same, but the commands below use a placeholder file, so they are not run automatically.

## 1. Check your inputs

You need a genome FASTA and a GFF3 annotation with `gene`, `mRNA`, and `CDS` features for each genome.

- FASTA sequence names must match the first column of the GFF3.
- Leave out organelles. They use other genetic codes, and a trans-spliced gene such as plastid rps12 can stop `extract_cds`. For a RefSeq GFF3: `awk -F'\t' '$3=="region" && /genome=(chloroplast|mitochondrion)/{drop[$1]=1} !($1 in drop)' in.gff > nuclear.gff`
- Use annotations of comparable quality. Errors in either one can look like selection.
- For enrichment, annotate both genomes with the same term database version. Different versions can give wrong results.
- For polyploids, consider splitting subgenomes first, so that best hits occur between matching subgenomes: `dndsr split_comparisons -c comparison_file.txt`.

## 2. Write a comparison file

Use one line per comparison, separated by spaces or tabs. Use double quotes around paths containing spaces.

```
sp1_v_sp2 data/sp1.fasta data/sp1.gff3 data/sp2.fasta data/sp2.gff3
```

Paths are relative to the folder where you run `dndsr`. If your data is elsewhere, reinstall with `--bind` (see [installation options](../../INSTALL.md)) and use absolute paths, or run from a parent folder that contains your data.

## 3. Run

On a cluster, submit `calculate_dnds` through your scheduler, because it can take many hours on whole genomes.

```
dndsr extract_cds -c comparison_file.txt -t 8
dndsr calculate_dnds -c comparison_file.txt -t 40
dndsr append_annotations -c comparison_file.txt -t 8
```

## 4. Check before you interpret

- **Look at exclusions** in `<comparison>/<genome>_CDS_qc.tsv`. Many exclusions usually mean an annotation or FASTA problem.
- **Count pairs.** A very low number of orthologous pairs compared with your gene count suggests mismatched inputs or very distant genomes.
- **Look at the dS distribution.** `calculate_dnds` reports how many pairs are undefined, below dS 0.02, above dS 0.8, and above 3 times the median dS. Very low dS makes ratios unstable. High dS brings saturation and wrong pairings. See [interpretation](../interpretation.md).
- **Count genes above 1.** `calculate_dnds` says so when there are fewer than 30, because enrichment then has little power.

## 5. Continue

```
dndsr ipr_enrichment -c comparison_file.txt -t 8
dndsr go_enrichment -c comparison_file.txt -t 8
dndsr term_enrichment -c comparison_file.txt -t 8
dndsr covariate_enrichment -c comparison_file.txt -t 8
dndsr dnds_ideogram -c comparison_file.txt
dndsr dnds_summary -c comparison_file.txt --regions-bed regions.bed
```

`ipr_enrichment` and `go_enrichment` correct `p_adj` once across all their tables by default. To correct each table separately, add `--adjust-scope per_type` to `ipr_enrichment` or `--adjust-scope per_ontology` to `go_enrichment`. To get `p_adj` with no correction from `go_enrichment`, add `--p-adjust none`. See [interpretation](../interpretation.md) before reporting terms.

Short genes and genes with small dS pass dN/dS > 1 far more often, and the first three commands can call a term significant because of the kind of gene it holds. `covariate_enrichment` adjusts for length and dS, and is the one to report a term list from: see [its calibration](../validation/enrichment-confounding.md). It writes `<comparison>_<q or s>_<GO or IPR>_covariate_enrichment.tsv`. It reads `GO:` and `IPR` terms from the annotated table, and needs v0.4.0 or later. Its evidence comes from simulation at 4 plant pairs, and it has not been validated for detecting real selection, so its terms remain leads. `go_enrichment` and `term_enrichment` give a calibration warning on every run, and record it in a `_notice.json` file beside their tables.

To contrast comparisons that share a genome, make a contrast file with `name comparison_A side_A comparison_B side_B`, and run:

```
dndsr dnds_contrast -c comparison_file.txt --contrast-file contrast_file.txt
```

`dnds_summary` and `dnds_contrast` read `--regions-bed` as BED by default: 0-based start, end not included. Use `--regions-coord gff1` with `dnds_contrast` for 1-based inclusive coordinates.

## Re-running

Existing results are kept unless you add `--overwrite`. dndsR tells you when an existing result no longer matches its inputs, using content hashes rather than file dates. Results record their software versions and parameters in `<comparison>_dnds_provenance.json`. Keep that file with your results and report it in your methods.
