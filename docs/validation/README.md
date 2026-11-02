# Validation evidence

What has actually been checked, and how. Commit SHAs are in `git log` and the progress log (`docs/dev/PROGRESS.md`). Provenance relies on content hashes, not file dates.

| Topic | Evidence | Where |
|---|---|---|
| Public genome pair | Old against current code, with gffread, PAML, and KaKs_Calculator | [public-genome.md](public-genome.md) |
| Other real genomes | Method and public results only | [real-genomes.md](real-genomes.md) |
| Validation rules | What must hold before a release | [criteria.md](criteria.md) |
| Installation routes | What ran, and what could not | [installation.md](installation.md) |
| CDS assembly | Output equals designed sequences | [cds-extraction.md](cds-extraction.md), `tests/testthat/test-example-pipeline.R` |
| dN/dS estimation | Compared with PAML, KaKs_Calculator, independent NG86, and simulated truth | [estimators.md](estimators.md) |
| Enrichment calibration | False-positive rates before and after the correction, on real public annotations, under three nulls; GFF recovery | [enrichment.md](enrichment.md) |
| Thread scaling | Time and memory per step, and identical results across thread counts | [performance.md](performance.md) |
| Fisher tests | p-values equal an independent hypergeometric calculation | `test-example-pipeline.R` |
| Region summaries and contrasts | Counts, means, deltas, and Wilcoxon p-values recalculated in base R | `test-example-pipeline.R` |
| Coordinates, filters, parent-child logic | Unit tests with hand-derived results | `tests/testthat/test-statistics.R` |
| GO and InterPro p-values, BH, descendants, ideogram windows | Hypergeometric and hand-calculated answers | `tests/testthat/test-enrichment-analytic.R` |
| CLI behavior and exit codes | Unit tests and real subprocess runs | `tests/testthat/test-cli.R` |
| Installation and tutorial | Commands run from the tutorial in a clean `HOME` | `tools/dev/test-tutorial.sh` |
| Fixed defects | Each defect mapped to the test that catches it | [regression-coverage.md](regression-coverage.md) |
| R tutorial | Every block is run | `tools/dev/test-r-tutorial.R` |
| User self-check | `dndsr verify_example`, and a test that it catches wrong values | `R/example.R` |

Not validated: [known issues](known-issues.md); ideogram appearance beyond rendering checks.
