# Analyzing your own genomes

This assumes you have completed [getting started](getting-started.md). The steps are the same, but the commands below use a placeholder file, so they are not run automatically.

## 1. Check your inputs

You need a genome FASTA and a GFF3 annotation with `gene`, `mRNA`, and `CDS` features for each genome.

- FASTA sequence names must match the first column of the GFF3.
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
- **Look at the dS distribution.** Very high dS suggests saturation, and very low dS makes ratios unstable. See [interpretation](../interpretation.md).

## 5. Continue

```
dndsr ipr_enrichment -c comparison_file.txt -t 8
dndsr go_enrichment -c comparison_file.txt -t 8
dndsr term_enrichment -c comparison_file.txt -t 8
dndsr dnds_ideogram -c comparison_file.txt
dndsr dnds_summary -c comparison_file.txt --regions-bed regions.bed
```

To contrast comparisons that share a genome, make a contrast file with `name comparison_A side_A comparison_B side_B`, and run:

```
dndsr dnds_contrast -c comparison_file.txt --contrast-file contrast_file.txt
```

`dnds_summary` and `dnds_contrast` read `--regions-bed` as BED by default: 0-based start, end not included. Use `--regions-coord gff1` with `dnds_contrast` for 1-based inclusive coordinates.

## Re-running

Existing results are kept unless you add `--overwrite`. dndsR tells you when an existing result no longer matches its inputs, using content hashes rather than file dates. Results record their software versions and parameters in `<comparison>_dnds_provenance.json`. Keep that file with your results and report it in your methods.
