# Interpreting dndsR results

## What dN/dS is

dN is the rate of amino-acid-changing substitutions per available site, and dS is the rate of silent substitutions per available site. Their ratio is often read as:

| dN/dS | Usual reading |
|---|---|
| Well below 1 | Most amino-acid changes were removed: purifying selection |
| About 1 | Consistent with neutral evolution, or a mix of selection pressures |
| Above 1 | Consistent with positive selection |

These are readings, not proofs. A whole-gene ratio averages over all sites and the whole time since the two genomes diverged, so most genes under positive selection at a few sites still have ratios below 1, and a ratio above 1 can arise without selection.

## A ratio above 1 is not proof of positive selection

dndsR's `dNdS > 1` is a threshold on a point estimate with no uncertainty attached. In our simulations, 600-codon genes evolving under true neutrality had an estimated ratio above 1 in up to 35% of cases with the default method (see [estimator validation](validation/estimators.md)). In a genome with many thousands of nearly neutral genes, this can be a large number of genes.

Treat dN/dS above 1 as a candidate list. To support a claim about particular genes, use methods that model selection, such as likelihood-ratio tests with more than two sequences, and look at the alignments.

Established methods can also disagree about which side of 1 a gene falls on.

## Before you trust a value

- **Small dS.** With few synonymous changes, the ratio is unstable: one substitution can change it greatly. Look at dN and dS, not just the ratio. dndsR does not filter on dS by default. You can add one, for example `--filter-expr "dS > 0.01"` in the enrichment commands.
- **Large dS.** Above roughly 1 to 2, synonymous sites are close to saturation and dS is unreliable. High dS can also make a ratio look conserved.
- **dS = 0.** The ratio is reported as NA and left out of enrichment tests, even if dN is greater than 0. Such genes are not reported as positive.
- **The method.** The default Comeron method tends to overestimate dS and underestimate dN relative to our simulated truth, which pulls ratios down. For example, a true ratio of 2 gives a median near 1.45 to 1.9. It does not account for codon usage. This makes `> 1` somewhat conservative but does not remove false positives.
- **The site model.** The Comeron method assigns sites to degeneracy classes. A gene designed with only amino-acid changes can therefore still get a small nonzero dS.

## Pipeline assumptions

- **Orthologs.** Pairs are reciprocal best protein hits, a heuristic. Paralogs, gene losses, and whole-genome duplications can pair the wrong genes. A wrong pairing can look like fast evolution.
- **One isoform.** Each gene uses its longest valid isoform by default. That may not be the isoform under selection, and differences between genomes in which isoform is annotated can inflate divergence.
- **Annotation errors.** Frame errors, gene-model errors, and undetected pseudogenes can inflate dN. CDS with in-frame stops or invalid lengths are excluded and recorded in the QC table, but other errors cannot be detected.
- **Alignment.** Each pair is aligned alone by global protein alignment and converted back to codons.
- **Recombination, gene conversion, and heterogeneous base composition** are not modeled.

## Enrichment tests

- Each row is treated as independent. Genes are not independent in practice (tandem duplicates, shared function), so p-values can be too small.
- `go_enrichment` defaults to topGO's `weight01` algorithm with unadjusted p-values, a common topGO practice, because its tests depend on each other. Applying BH to them has no standard FDR meaning.
- `min_pos` filters terms using the positive count, which is part of the test statistic. Terms removed this way shrink the set of tests corrected for multiple testing, which can make adjusted p-values less conservative.
- `max_dnds` (default 10) drops the most extreme ratios, which are usually artifacts but occasionally real, and the enrichment commands do not report how many.
- Compare only genomes annotated with the same term database version.

## Regional summaries and contrasts

- Genes are tested as independent, although neighboring genes share history. Region p-values can be too small.
- `dnds_contrast` compares genes in one genome across two comparisons, so it is confounded by each comparison's divergence time. No correction is applied across contrasts.
- `dnds_state_contrast --level gene` is disabled. See [known issues](validation/known-issues.md).

## What to report

The dndsR version and commit, image digest, command and options (recorded in `*_provenance.json`), the numbers of genes excluded, and how you treated small and large dS.

## Inappropriate uses

- Very closely related genomes or strains: few substitutions, and differences may be polymorphism rather than fixed.
- Very divergent genomes: dS saturates, and reciprocal best hits become unreliable.
- Sequences without reliable gene models.
- Calling individual genes as positively selected from this ratio alone.

## Citing methods

Cite dndsR (see [CITATION.cff](../CITATION.cff)), orthologr (run `citation("orthologr")` in R), DIAMOND, and the estimator: Comeron JM (1995) *J Mol Evol* for the default method.
