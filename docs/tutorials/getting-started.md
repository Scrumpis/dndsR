# Getting started with dndsR

This tutorial installs dndsR on Linux with Singularity, runs a complete small example, and shows how to read the results. Every command block here is run automatically by `tools/dev/test-tutorial.sh`.

You need bash, git, Singularity, and network access. The image is a 1.3 GB download. If you use Apptainer, replace `singularity` with `apptainer`. For other setups, see [installation options](../../INSTALL.md).

## 1. Get the code

```
git clone --branch dev https://github.com/Scrumpis/dndsR
cd dndsR
```

## 2. Get the container

This digest always identifies the same image.

```
singularity pull dndsr.sif docker://scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831
```

## 3. Install the `dndsr` command

This installs into your home directory and needs no administrator rights. Keep the cloned folder.

```
tools/dndsr-install install --engine singularity --sif dndsr.sif
export PATH="$HOME/.local/bin:$PATH"
```

To make the `PATH` change permanent, add the `export` line to your shell's startup file, such as `~/.bashrc`.

Check that it works. You should see a list of commands.

```
dndsr --help
```

## 4. Get the example

The example is two small synthetic genomes with designed substitutions, so the right answers are known in advance. Its sequences have no biological meaning.

```
dndsr example
cd dndsr_example
cat comparison_file.txt
```

The comparison file has one line per pairwise comparison: a name, then the query FASTA and GFF3, then the subject FASTA and GFF3.

## 5. Extract coding sequences

```
dndsr extract_cds -c comparison_file.txt
```

This creates a folder named after the comparison, `example_AvB`, with one CDS FASTA per genome and a QC table. Each gene is represented by its longest valid isoform.

```
grep -c '>' example_AvB/spA_CDS.fasta
head -1 example_AvB/spA_CDS_qc.tsv
```

`spA` has 35 genes with coding sequences, and 34 are written. The one excluded, `spA_g32`, has a length that is not a multiple of 3. See the `written` column in the QC table.

## 6. Calculate dN/dS

```
dndsr calculate_dnds -c comparison_file.txt -t 2
```

This finds reciprocal best hits with DIAMOND, aligns each pair, and estimates dN and dS. The result is `example_AvB/example_AvB_dnds.tsv`.

```
cut -f1-5 example_AvB/example_AvB_dnds.tsv | head -5
```

| Gene | dN | dS | dN/dS | What it means here |
|---|---:|---:|---:|---|
| g01 | 0 | 0 | NA | Identical sequences, so the ratio is undefined |
| g02 | 0 | 0.14 | 0 | Only synonymous changes |
| g20 | 0.09 | 0.02 | 5.6 | Designed with far more nonsynonymous changes |

Small differences in dN and dS from the second decimal are expected between estimators, but the same version on the same input gives the same values.

## 7. Verify

```
dndsr verify_example
```

This should end with `All 14 example checks passed.` It compares your results with the designed sequences, an independent calculation, and a reference run.

## 8. Add annotations and test for enrichment

```
dndsr append_annotations -c comparison_file.txt
dndsr term_enrichment -c comparison_file.txt
ls example_AvB/*KEGG*
```

`example_AvB_q_KEGG_enrichment.tsv` tests whether a term is over-represented among genes with dN/dS above 1, on the query side. The example's terms are real identifiers assigned arbitrarily. Term `K00001` was assigned to the genes designed with excess nonsynonymous changes, so it should be enriched.

The next example command relaxes a filter that suits real genomes but not this tiny one: by default, terms that annotate more than 20% of genes are skipped as too general.

```
dndsr ipr_enrichment -c comparison_file.txt --max-prop 1
ls example_AvB/*IPR*tsv
```

## 9. Summarize a region

```
cat regions.bed
dndsr dnds_summary -c comparison_file.txt --regions-bed regions.bed
cut -f2-4,10-11 example_AvB/example_AvB_dnds_summary.tsv
```

The `region` rows are genes overlapping the BED interval, and `background` rows are all other genes.

## Next

- Read [how to interpret results](../interpretation.md) before drawing biological conclusions. A dN/dS above 1 is not proof of positive selection.
- Run [your own genomes](your-own-data.md).
- See what each [output file](../output-files.md) contains.
