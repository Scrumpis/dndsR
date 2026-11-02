# Draft release notes, next minor version

Draft only. Nothing is tagged or published. Fill in the version and image digest at release.

**This release adds one new command. It changes no default, algorithm, or existing output value.**

## New

- **`dndsr covariate_enrichment`**, and `covariate_enrichment()`: GO and InterPro enrichment adjusted for gene length and dS. Under simulated evolution with no link between omega and annotation, Fisher's test reported a false term in up to 70% of runs, and `go_enrichment`'s defaults in up to 43%. This command stayed at or below 5% at 4 public plant pairs, including a pair held out of its selection. Terms remain leads, not evidence of selection. [Evidence](../validation/enrichment-confounding.md).
- `calculate_dnds` reports undefined, low-dS, high-dS, and outlier pairs, and says when few genes exceed 1. This is a message only.
- `extract_cds` explains a failure from trans-spliced genes, and suggests removing organelles.
- A [full tutorial on public genomes](../tutorials/public-genomes.md), and a [checklist for publication](../publication-protocol.md).

## Revised evidence

`go_enrichment`'s excess of false terms is larger than v0.3.0's documents say: 43% of simulated runs at one pair, where they say 9.5%. [Estimator](../validation/estimators.md) and [ortholog](../validation/orthology.md) benchmarks are new.

## Faster

`ipr_enrichment` about 25 times and `append_annotations` about 2.5 times, with byte-identical tables on the public pair.

## Install

See [INSTALL.md](../../INSTALL.md). Image digest: fill in at release.
