# dndsR 0.4.0 release notes

Draft until tagged. Fill in the image digest once the publish workflow has run.

**This release adds one new command and two warnings. It changes no default, algorithm, or the content of an existing table.**

## New

- **`dndsr covariate_enrichment`**, and `covariate_enrichment()`: GO and InterPro enrichment adjusted for gene length and dS. Under simulated evolution with no link between omega and annotation, Fisher's test reported a false term in up to 70% of runs, and `go_enrichment`'s defaults in up to 43%. This command stayed at or below 5% at 4 public plant pairs, including a pair held out of its selection. Those simulations do not establish false discovery control in general: polyploids and other clades are untested, and it has not been validated for detecting real selection in real genomes. Terms remain leads, not evidence of selection. [Evidence](../validation/enrichment-confounding.md).
- **`go_enrichment` and `term_enrichment` warn on every run** that they do not adjust for length or dS, and record it in a `_notice.json` file beside their tables. Defaults and tables are unchanged.
- `term_enrichment` and `ipr_enrichment` stop if `pos_threshold` is not one number. A value such as `q0.95` used to give a run with no positive genes.
- `calculate_dnds` reports undefined, low-dS, high-dS, and outlier pairs, and says when few genes exceed 1. This is a message only.
- `extract_cds` explains a failure from trans-spliced genes, and suggests removing organelles.
- A [full tutorial on public genomes](../tutorials/public-genomes.md), and a [checklist for publication](../publication-protocol.md).

## Revised evidence

`go_enrichment`'s excess of false terms is larger than v0.3.0's documents say: 43% of simulated runs at one pair, where they say 9.5%. [Estimator](../validation/estimators.md) and [ortholog](../validation/orthology.md) benchmarks are new.

## Faster

`ipr_enrichment` about 25 times and `append_annotations` about 2.5 times, with byte-identical tables on the public pair.

## Install

See [INSTALL.md](../../INSTALL.md). Image digest: fill in at release.
