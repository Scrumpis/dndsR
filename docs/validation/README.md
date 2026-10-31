# Validation evidence

What has actually been checked, and how. Commit SHAs are in `git log` and the progress log (`docs/dev/PROGRESS.md`). Provenance relies on content hashes, not file dates.

| Topic | Evidence | Where |
|---|---|---|
| CDS assembly | Output equals designed sequences | [cds-extraction.md](cds-extraction.md), `tests/testthat/test-example-pipeline.R` |
| dN/dS estimation | Compared with PAML, KaKs_Calculator, independent NG86, and simulated truth | [estimators.md](estimators.md) |
| Fisher tests | p-values equal an independent hypergeometric calculation | `test-example-pipeline.R` |
| Region summaries and contrasts | Counts, means, deltas, and Wilcoxon p-values recalculated in base R | `test-example-pipeline.R` |
| Coordinates, filters, parent-child logic | Unit tests with hand-derived results | `tests/testthat/test-statistics.R` |
| CLI behavior and exit codes | Unit tests and real subprocess runs | `tests/testthat/test-cli.R` |
| Installation and tutorial | Commands run from the tutorial in a clean `HOME` | `tools/dev/test-tutorial.sh` |
| User self-check | `dndsr verify_example`, and a test that it catches wrong values | `R/example.R` |

Not validated: real genomes; [known issues](known-issues.md); `go_enrichment` p-values beyond confirming it runs, and `dnds_ideogram` beyond confirming files are written.
