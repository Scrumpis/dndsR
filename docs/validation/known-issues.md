# Known issues

Open issues that affect results or use. Fixed defects are in [NEWS](../../NEWS.md).

## Scientific

| Issue | Status | What to do |
|---|---|---|
| `dnds_state_contrast --level gene`: a code audit found that genes from separate comparison files are not aggregated, one regions file is applied to both genomes, and when both sides fall into the two groups, effects can be zero by construction | Experimental, not validated, prints a message. It also stops on genomes whose sequence names do not show a subgenome. The rest is not confirmed by running; redesign is a scientific decision | Do not use for conclusions |
| `dnds_state_contrast` regions use a different coordinate convention from `dnds_summary`: BED starts are not converted | Documented in its help | Check your coordinates |
| Zero-cell state contrasts produce no row, shrinking the multiple-testing family; mixed-model warnings are suppressed | Open | Treat `dnds_state_contrast` results as exploratory |
| Default CDS with a first exon of nonzero GFF phase (partial genes) are not trimmed to frame | Open; they are often excluded by the in-frame stop check | Check annotation quality |
| `split_comparisons --mode haplotype` treats any trailing `_<number>` as a haplotype, so `scaffold_12` becomes haplotype 12 | Open | Check the labels it reports |
| Term filters using `min_pos` select on the test statistic before multiple-testing correction | Open | See [interpretation](../interpretation.md) |
| Fixed 1.96 in normal confidence intervals; dependence between genes; no adjustment across contrasts | Documented | See [interpretation](../interpretation.md) |

Interpretation limits and the estimator's bias are covered in [interpretation](../interpretation.md) and [estimator validation](estimators.md).

## Software

| Issue | Status |
|---|---|
| The published image, `v0.2.0`, predates this code. It works with the current code on the tested example. A new image is needed for each release | Needs maintainer approval to publish |
| The updated `tools/Dockerfile`, which includes dndsR in the image, has not been built: no Docker is available in the development environment. A CI workflow builds it | Runs when `dev` is pushed |
| Docker, Apptainer, macOS, and native installation are untested | Open |
| A `git pull` in the clone reinstalls dndsR automatically on the next run, so the analysis version can change silently. Provenance records the version | Open; pin a tag for production analyses |
| Windows runs sequentially and has not been tested | Open |
| `inst/fonts/ArialBold.ttf`: Arial is proprietary and redistribution rights are unconfirmed | Needs maintainer decision before release or CRAN |
