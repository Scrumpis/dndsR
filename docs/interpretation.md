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

dndsR's `dNdS > 1` is a threshold on a point estimate with no uncertainty attached. In simulations with realistic codon usage, **about half of truly neutral genes exceed 1** at ordinary divergence, as any roughly unbiased estimate would ([evidence](validation/estimators.md)). In a genome with many thousands of nearly neutral genes, this can be a large number of genes. Established methods can also disagree about which side of 1 a gene falls on.

A bootstrap interval would not fix this. A 95% interval from resampling codons covered the true ratio about 90% of the time for lower ratios at moderate dS, and almost never where it matters most: at a true ratio of 3 and dS 0.4 or above, or a true ratio of 1 or above at dS 0.8 when rates vary among sites. So dndsR does not offer one ([evidence](validation/estimators.md#benchmark-against-independent-estimators)).

Treat dN/dS above 1 as a candidate list. To support a claim about particular genes, use methods that model selection, such as likelihood-ratio tests with more than two sequences, and look at the alignments.

## Before you trust a value

- **Small dS.** With few synonymous changes, the ratio is unstable: one substitution can change it greatly. Look at dN and dS, not just the ratio. dndsR does not filter on dS by default. You can add one, for example `--filter-expr "dS > 0.01"` in the enrichment commands.
- **Large dS.** In simulations, a true ratio of 1 read 0.89 at dS 0.8 and 0.74 at dS 1.5, and a true ratio of 3 read 0.38 of itself at dS 0.8. Genes with a true ratio of 1.5 passed 1 in 62–94% of replicates at dS 0.8, and 9.5% or fewer above 400 codons at dS 1.5. So high divergence can hide a ratio above 1 ([evidence](validation/estimators.md)). Pairs at high dS are also often the wrong genes: over half of A. thaliana and A. lyrata pairs with OMA evidence at dS 0.3–0.8 disagreed with OMA, and 92% at 0.8–1.5 ([evidence](validation/orthology.md)).
- **Too few genes above 1.** Between A. thaliana and increasingly distant relatives, 206 genes passed 1 at median dS 0.14, 17 at 0.27, and 3 at 0.44. Enrichment on `dNdS > 1` loses power with divergence, so a short list says little.
- **dS = 0.** The ratio is reported as NA and left out of enrichment tests, even if dN is greater than 0. Such genes are not reported as positive.
- **The method.** The default Comeron method compresses ratios: in simulations with realistic codon usage, it was accurate near a ratio of 0.5 to 1, read high at very low ratios, and read low at high ones (a true ratio of 2 gave about 1.6 to 1.9), and more so at high divergence. Related methods do the same. It does not account for codon usage.
- **The site model.** The Comeron method assigns sites to degeneracy classes. A gene designed with only amino-acid changes can therefore still get a small nonzero dS.

## Pipeline assumptions

- **Orthologs.** Pairs are reciprocal best protein hits, a heuristic. Paralogs, gene losses, and whole-genome duplications can pair the wrong genes. A wrong pairing can look like fast evolution. Compared with OMA, 3.6% of pairs with evidence disagreed within *Arabidopsis*, and 15% against triplicated *B. rapa*. A pair above 3 times the median dS disagreed 85% and 64% of the time ([evidence](validation/orthology.md)).
- **One isoform.** Each gene uses its longest valid isoform by default. That may not be the isoform under selection, and differences between genomes in which isoform is annotated can inflate divergence.
- **Annotation errors.** Frame errors, gene-model errors, and undetected pseudogenes can inflate dN. CDS with in-frame stops or invalid lengths are excluded and recorded in the QC table, but other errors cannot be detected.
- **Alignment.** Each pair is aligned alone by global protein alignment and converted back to codons.
- **Recombination, gene conversion, and heterogeneous base composition** are not modeled.

## Enrichment tests

**What a result means.** A significant term says that genes with a high dN/dS carry that term more often than other annotated genes in this table. It is not evidence that the term's function evolved under positive selection, and it does not show that any single gene did. Treat it as a lead to test by other means, such as codon-site or branch-site models and independent data.

**Raw and adjusted values.**

- `p_value` and `raw_p` are not corrected for multiple testing. With thousands of terms, many fall below 0.05 by chance.
- In `term_enrichment`, `p_adj` is corrected across `n_tested` terms by `p_adjust_method`.
- In `ipr_enrichment`, by default, `p_adj` is corrected once across all eligible InterPro entry types, so `n_tested` is the same in every table and the tables can be read together. On a real genome pair with no true association, 1–3% of runs had a significant term. `--adjust-scope per_type` in the CLI, or `adjust_scope = "per_type"` in R, corrects within each type: 9–17% of runs then had a significant term in some table, although each table alone was near 5%, so do not pool per-type terms as if corrected together.
- In `go_enrichment`, `p_adj` by default is corrected by BH once across every tested term in all requested ontologies, and `raw_p` is not. `--adjust-scope per_ontology` corrects within each ontology, and `--p-adjust none` copies `raw_p` to `p_adj`. topGO's `weight01` tests depend on each other, so BH has no proven error rate for them. Measured on a real genome pair with no true association, 1% of runs had a significant term under two null models, and **9.5%** under a third that keeps dN/dS linked to gene length and dS. Most of that excess came from one broad term, extracellular region. So a large term among short or low-dS genes can pass this correction with no link to selection: see [enrichment calibration](validation/enrichment.md#go-enrichment).
- **BH correction does not guarantee nominal error control under all conditions: it assumes p-values that are valid under the null, and enrichment p-values do not account for gene length or dS.** On the null above, 9.5% of GO runs had a false term, against a nominal 5%. Significant enrichment terms remain exploratory: they are not independently validated evidence of selection.
- These rates are false positives when *nothing* is associated. BH controls the false discovery rate, the expected share of false terms in a list, not the chance of any false term, and none of these tests measures it when real associations exist.

**What selects the candidates.** Positive genes are not a random sample of genes, so a term can be enriched with no link to selection:

- **Gene length.** Short genes have noisier ratios and pass dN/dS > 1 far more often. Terms common among short genes, such as secreted peptides, are favored.
- **dS.** The ratio is unstable when dS is small, so genes with few synonymous differences pass more easily. At high divergence, dS saturates.
- **Annotation completeness.** Unannotated genes take no part in these tests, and in practice many positive genes are unannotated. Well-studied gene families are favored.
- **Dependence.** Tandem duplicates and genes sharing a term are not independent, although Fisher's test assumes they are.
- **The candidate rule.** Results change with `pos_threshold`, `max_dnds`, and filters. Check whether a conclusion survives other reasonable choices.

`go_enrichment`, `ipr_enrichment`, and `term_enrichment` do not correct for length, dS, or annotation, and passing correction does not remove these effects: multiple-testing correction limits false positives when there is no association, and says nothing about why an association exists. See [enrichment calibration](validation/enrichment.md). **This is larger than the 9.5% above suggests.** When sequences were simulated with omega independent of annotation, and a real estimator produced the labels, 43% of `go_enrichment` runs had a false term between *A. thaliana* and *A. lyrata*, and 70% with Fisher's test on all GO ancestors.

`covariate_enrichment` adjusts for length and dS, keeping false-term rates at or below nominal in all those tests, including new seeds and a held-out pair from another plant family, with power close to Fisher's test at known signals ([evidence](validation/enrichment-confounding.md)). Terms near the threshold change with the annotation release, so check that a term is stable before you rely on it. It does not adjust for annotation completeness or dependence between genes, and it doesn't make a term evidence of selection. Terms remain leads.

**Also:**

- `min_total` and `max_prop` choose what is tested, and `min_pos` chooses what is reported, after correction.
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
