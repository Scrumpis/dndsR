# Regression coverage of fixed defects

Each fixed defect, and the test that fails if it returns. Tests are in `tests/testthat/` unless noted. Expected values come from outside dndsR.

| Defect | Test | Independent basis |
|---|---|---|
| Wrong exon order, minus strand | `test-example-pipeline.R`: "multi-exon minus-strand genes keep their exon order", "every written CDS matches the designed sequence exactly" | Sequences designed before the genomes were generated |
| Isoform chimeras | Same, plus "the QC tables record exclusions and isoform choices", and `test-select-isoform.R` | Designed sequences |
| Invalid longest isoform chosen over a valid one | `test-select-isoform.R`: "an invalid longest isoform loses to a shorter valid one", plus gene g36 in the example | Designed gene |
| Internal stops and bad lengths reaching dN/dS | "the QC tables...", "pairs with excluded or missing sequences give no result", `test-cds-quality.R` | Designed genes |
| Stale outputs from older versions reused | "existing outputs without a QC table are not reused", "a provenance record identifies the inputs" | Constructed state |
| Wrong-genome annotations | `test-statistics.R`: "a shared gene ID gets annotations from its own genome" | Constructed GFFs |
| Empty annotations, RefSeq-style IDs | `test-statistics.R`: "RefSeq-style IDs connect sequences to their annotations", "identifier prefixes are removed only at the start" | Constructed GFF from the example |
| Terms lost after encoded commas | `test-statistics.R`: "every term is found in GFF3 lists, including encoded separators" | Attribute text in the real GFF3 format |
| Terms matched inside identifiers | `test-statistics.R`: "feature names are not mistaken for functional terms" | Constructed attributes |
| Headerless BED lost a region | `test-statistics.R`: "headerless BED files keep their first region" | Hand-derived coordinates |
| Failed filters running unfiltered, NA rows | `test-statistics.R`: "filter expressions fail loudly and NA rows are removed" | Constructed table |
| InterPro parent-child logic | `test-enrichment-analytic.R`: "the parent-child method removes genes explained by a significant child"; `test-statistics.R`: propagation test | Hypergeometric answers |
| Unknown-type terms in every stratum | `test-enrichment-analytic.R`: "terms of unknown type are not tested in any stratum" | Constructed case |
| GO counts, descendants, and definition of enrichment | `test-enrichment-analytic.R`: "GO tests count descendants and match hand calculation" | Hypergeometric answers |
| Fisher tests and BH | `test-enrichment-analytic.R`, and "downstream steps run and match independent calculations" | Hypergeometric answers |
| Regions, summaries, contrasts | `test-example-pipeline.R`: "region summaries count the right genes", "contrasts have the right direction and pairing" | Base R recalculation |
| Pairwise state contrast | `test-statistics.R`: "pairwise state contrasts match a directly fitted logistic model" | Hand-calculated odds ratio |
| Gene-level state contrast giving no information | `test-statistics.R`: "gene-level state contrast fails instead of returning meaningless values" | Records why it is disabled |
| Ideogram binning, boundary at 1 | `test-enrichment-analytic.R`: "ideogram windows average the right genes" | Hand calculation |
| Offline InterPro failure | `test-enrichment-analytic.R`: "InterPro enrichment works with no internet access" | Blocked network |
| Worker failures reported as success, wrong exit codes | `test-cli.R`: "failed workers...", "a worker that dies...", "the command line exits nonzero on failure" | Real subprocesses |
| Values altered by the parser, unknown options ignored | `test-cli.R` | Constructed arguments |
| Tutorial facts going stale | `tools/dev/test-tutorial.sh` | Runs the commands and checks stated facts |
| Container problems | `tools/dev/container-checks.sh` in CI | Runs documented workflows in the image |

`dndsr verify_example` also gives users 14 of these checks, and a test confirms it fails on wrong values. `tools/dev/validation/` holds the real-genome comparisons, which are run before a release rather than on every commit, because they need genomes and a cluster.
