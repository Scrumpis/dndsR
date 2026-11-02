# Using dndsR from R

Every command-line step is an ordinary R function with the same options, so you can call dndsR from scripts, notebooks, or RStudio. Help is available as usual, for example `?calculate_dnds`.

You need dndsR and its dependencies (see [installation options](../../INSTALL.md)), or an R session inside the dndsR container. DIAMOND must be on the `PATH` for `calculate_dnds`.

## Run the example

```r
library(dndsR)

folder <- dndsr_example(file.path(tempdir(), "dndsr_example"))

# The comparison file names its inputs relative to the working directory.
setwd(folder)
extract_cds(comparison_file = "comparison_file.txt")
calculate_dnds(comparison_file = "comparison_file.txt", comp_cores = 2)
verify_example()
```

The last call ends with `All 14 example checks passed.`

## Work with the results

Results are plain tab-separated tables.

```r
dnds <- read.delim("example_AvB/example_AvB_dnds.tsv")
head(dnds[, c("query_id", "subject_id", "dN", "dS", "dNdS")])

# Genes with dN/dS above 1, keeping dS and dN visible
subset(dnds, dNdS > 1, c(query_id, dN, dS, dNdS))
```

## Annotate and test

```r
append_annotations(comparison_file = "comparison_file.txt")
term_enrichment(comparison_file = "comparison_file.txt", make_plots = FALSE)

read.delim("example_AvB/example_AvB_q_KEGG_enrichment.tsv")
```

`append_annotations()` with `dnds_file`, `query_gff`, and `subject_gff` returns a data frame instead of writing a file, which is convenient in a script.

## Single comparisons without a comparison file

```r
extract_cds(
  comparison_name = "my_pair",
  query_fasta = "spA.fasta", query_gff = "spA.gff3",
  subject_fasta = "spB.fasta", subject_gff = "spB.gff3"
)

calculate_dnds(
  comparison_name = "my_pair",
  query_seq = "my_pair/spA_CDS.fasta",
  subject_seq = "my_pair/spB_CDS.fasta"
)
```

Before interpreting anything, read [how to interpret results](../interpretation.md).
