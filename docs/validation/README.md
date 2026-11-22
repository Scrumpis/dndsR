# Validation evidence

What has actually been checked, and how. Commit SHAs are in `git log` and the progress log (`docs/dev/PROGRESS.md`). Provenance relies on content hashes, not file dates.

| Topic | Evidence | Where |
|---|---|---|
| Public genome pair | Old against current code, with gffread, PAML, and KaKs_Calculator | [public-genome.md](public-genome.md) |
| Other real genomes | Method and public results only | [real-genomes.md](real-genomes.md) |
| Validation rules | What must hold before a release | [criteria.md](criteria.md) |
| Installation routes | What ran, and what could not | [installation.md](installation.md) |
| CDS assembly | Output equals designed sequences | [cds-extraction.md](cds-extraction.md), `tests/testthat/test-example-pipeline.R` |
| dN/dS estimation | Compared with PAML, KaKs_Calculator, independent NG86, and simulated truth, across 540 settings, with failure conditions and bootstrap coverage | [estimators.md](estimators.md) |
| Enrichment and length and dS | Labels from simulated evolution and a real estimator, on three public pairs; ten methods under a rule fixed before results; confirmation on new seeds and a held-out pair | [enrichment-confounding.md](enrichment-confounding.md) |
| Image stacks | The same code on two published images: pairs, values, and why they differ | [stack-comparison.md](stack-comparison.md) |
| Orthologs | Reciprocal best hits against OMA, by dS, and the effect on term lists | [orthology.md](orthology.md) |
| Enrichment calibration | False-positive rates before and after the correction and default changes, on real public annotations, under three nulls; GFF recovery | [enrichment.md](enrichment.md). Earlier evidence: [59faf30](enrichment-59faf30.md), [1a363fb](enrichment-1a363fb.md) |
| Phylogenetic wrapper (`phylo_busted`) | Exact agreement with HyPhy run directly on independently prepared inputs; HyPhy behaviors the design relies on; seed, thread, and version effects. Wrapper correctness only, not calibration | [phylo-busted-wrapper.md](phylo-busted-wrapper.md) |
| Phylogenetic preparation (`phylo_prepare`), end to end | Synthetic genomes with designed failures: alignments equal designed CDS and a true indel alignment; exclusions and their basis; OrthoFinder formats; gene-tree mode and the whole workflow equal HyPhy run directly; hypothesis guards. Engineering correctness only | [phylo-prepare.md](phylo-prepare.md) |
| Phylogenetic calibration pilot (P3) and real-genome benchmark | BUSTED with defaults is conservative under correct nulls; multinucleotide changes and misaligned windows inflate false positives severely; single-start fits miss optima under selection; on five public Brassicaceae genomes, alignment artifacts are concentrated among significant genes, and dndsR equals direct HyPhy. A pilot, not the calibration | [phylo-calibration-pilot.md](phylo-calibration-pilot.md) |
| Phylogenetic calibration, stage 2A | Four BUSTED configurations, triple hits, optimizer starts, and masking of flagged windows on simulated clean, misspecified, and realistically corrupted alignments: all calibrated under correct nulls; BUSTED-E loses much power; double hits fix multinucleotide changes; masking removes corruption-driven false positives without measurable power loss. Intermediate batch, not the full study | [phylo-calibration-stage2a.md](phylo-calibration-stage2a.md) |
| Enrichment error criterion (study N) | Complete-null FWER, mixed-null FDR, and power with labels drawn at known probabilities on a public pair: correction and term dependence controlled; the excess is confounding by length and dS, removed by `covariate_enrichment` | [enrichment-error-calibration.md](enrichment-error-calibration.md) |
| RELAX calibration, stage 1 (study R) | Version, gene-length, heterogeneity, and held-out size, power, and genome-wide BH of two bootstrap-calibrated p-values; the excess is HyPhy's and finite-sample; both corrections pass; confirmation pending | [relax-calibration-stage1.md](relax-calibration-stage1.md) |
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
