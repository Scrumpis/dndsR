# Methods: how dndsR computes its results

A plain description of what the software does, with its defaults, to support reproducing and reporting analyses. Validation evidence: [validation/](validation/README.md). Limits: [interpretation](interpretation.md).

## Coding sequences

`extract_cds` builds gene models from a GFF3 file with Bioconductor (`txdbmaker`, `GenomicFeatures`), and extracts each transcript's CDS in transcript order from the genome FASTA (`Rsamtools`, `Biostrings`).

By default, each gene is represented by its longest valid isoform. A valid CDS has a length divisible by 3 and no in-frame stop in the chosen genetic code. Ties are broken by transcript name. Other CDS are left out and recorded in a QC table. With `--group-by tx`, every transcript is kept instead.

## Orthologs and alignment

`calculate_dnds` calls `orthologr::dNdS()` with these defaults:

| Step | Default |
|---|---|
| Ortholog pairs | Reciprocal best hits of translated CDS, using DIAMOND `fast` mode, e-value 1e-5 |
| Protein alignment | Pairwise global alignment (Needleman-Wunsch) |
| Codon alignment | PAL2NAL |
| Estimator | Comeron (1995), through libsequence's `Comeron95` as included in orthologr |

Each pair is treated separately. dN/dS is reported as undefined (NA) where dS is 0.

## Annotation

`append_annotations` adds coordinates and terms from the GFF3 to each pair, matching gene IDs without prefixes such as `gene-` and `rna-`. Terms come from all attributes of the gene and its descendants, other than `ID` and `Parent`.

## Enrichment

Genes are called positive if `dNdS > pos_threshold` (default 1), among rows where dN/dS is defined and below `max_dnds` (default 10). This is a threshold on a point estimate, not a statistical test of selection.

| Command | Test | Multiple testing |
|---|---|---|
| `term_enrichment` | One-sided Fisher's exact test per term, positive against other genes | BH by default |
| `ipr_enrichment` | The same, separately for each InterPro entry type; optionally a parent-child method that removes genes explained by significant child terms | BH by default, within each type |
| `go_enrichment` | topGO `weight01` with Fisher's test, on genes with GO terms in that ontology | None by default |

`enrichment` is the proportion of positive genes with a term divided by the proportion of other genes with it.

## Summaries, contrasts, and ideograms

- `dnds_summary`: mean, median, and interval, for all genes and for those overlapping given regions, with a Wilcoxon rank-sum test of region against background.
- `dnds_contrast`: paired comparison of one genome's genes between two comparisons, using a signed-rank test and McNemar's test on calls above 1, with no correction across contrasts.
- `dnds_state_contrast`, pairwise level: logistic regression of a selection state on a group column. The gene level is disabled.
- `dnds_ideogram`: mean of log(dN/dS + 1) in genome windows, by the window containing a gene's start, with default windows of 300 kb.

## Reporting a dndsR analysis

Report the dndsR version and commit, image digest, and the provenance file `*_dnds_provenance.json`, which records software versions, parameters, and input hashes. Also cite orthologr, DIAMOND, and the estimator.
