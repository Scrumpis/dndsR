# Validation criteria

What must hold before a dndsR release, based on what can reasonably be demonstrated. Numbers below are tied to evidence in this folder. They are not targets to tune toward: a change that misses one is investigated and explained, not adjusted until it passes.

## Exact: any failure blocks a release

| Property | Check |
|---|---|
| CDS equal independently derived sequences | gffread on a public genome: 0 differences among compared genes. Designed example: exact |
| Designed zero-change cases | Identical sequences give dN = dS = 0; synonymous-only gives dN = 0 |
| Fisher p-values, BH, and enrichment definitions | Equal to the hypergeometric answer on constructed cases, to 1e-8 |
| Coordinates, regions, contrasts, and summary counts | Equal to base R recalculation |
| Determinism | Same inputs and version give byte-identical dN/dS tables across thread counts |
| Failure behavior | Failed steps exit nonzero; unknown options and missing inputs stop |
| Documented commands | Tutorials and container checks run as written |
| Disabled analyses | Stay disabled, with a test |

## Agreement: a miss is investigated, and blocks a release until explained

| Property | Current evidence |
|---|---|
| Comeron estimates against simulated truth, at realistic codon usage and dS up to about 0.3 | Ratio within about 15% of ω for ω from 0.1 to 2 ([estimators](estimators.md)) |
| Correlation with PAML and KaKs_Calculator on real pairs, on one shared alignment | Spearman on the ratio about 0.97 or higher in the public pair ([public genome](public-genome.md)) |
| Reference run on the bundled example | dN and dS within 0.001 between versions |

## Reported, not pass or fail

These describe behavior users need to know: the share of neutral genes above 1, power by gene length, and sensitivity to isoform choice. They are published in `docs/validation/`, and interpretive text must match them.

## Rules for evidence

- Expected values come from outside dndsR.
- Regression snapshots show stability, not correctness, and are labeled as such.
- Each piece of evidence names its commit, software versions, and inputs by hash.
- Unpublished data and findings stay out of Git. Public evidence uses public genomes or synthetic data.
- Real-genome comparisons rerun when CDS, pairing, annotation, or enrichment code changes.
