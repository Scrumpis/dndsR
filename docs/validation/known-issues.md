# Known issues

Open issues that affect results or use. Fixed defects are in [NEWS](../../NEWS.md).

## Scientific

| Issue | Status | What to do |
|---|---|---|
| `dnds_state_contrast --level gene`: confirmed by running that a single table gives no effect by construction (log odds 0, p = 1), and that genes in separate files are not combined. See `docs/dev/state-contrast-design.md` | **Disabled: stops with an error.** Redesign is a scientific decision | Use `--level pairwise` with a `group` column |
| `dnds_state_contrast` regions use a different coordinate convention from `dnds_summary`: BED starts are not converted | Documented in its help | Check your coordinates |
| Zero-cell state contrasts produce no row, shrinking the multiple-testing family; mixed-model warnings are suppressed | Open | Treat `dnds_state_contrast` results as exploratory |
| Default CDS with a first exon of nonzero GFF phase (partial genes) are not trimmed to frame | Open; they are often excluded by the in-frame stop check | Check annotation quality |
| `split_comparisons --mode haplotype` treats any trailing `_<number>` as a haplotype, so `scaffold_12` becomes haplotype 12 | Open | Check the labels it reports |
| Enrichment does not correct for gene length, dS, annotation completeness, or dependence between genes, all of which affect results. See [enrichment calibration](enrichment.md) | Documented; modeling is on the roadmap to v1.0.0 | See [interpretation](../interpretation.md) |
| Most of the default `go_enrichment` list at p < 0.05 cannot be told apart from noise. See [enrichment calibration](enrichment.md) | Documented | Use a stricter threshold, and validate terms independently |
| A dN/dS has no uncertainty attached, and `Comeron95` truncates undefined terms without flagging them | Open; see `docs/dev/ROADMAP-1.0.md` | See [interpretation](../interpretation.md) |
| Fixed 1.96 in normal confidence intervals; dependence between genes; no adjustment across contrasts | Documented | See [interpretation](../interpretation.md) |

Interpretation limits and the estimator's bias are covered in [interpretation](../interpretation.md) and [estimator validation](estimators.md).

## Software

| Issue | Status |
|---|---|
| The published image, `v0.2.0`, predates this code. It works with the current code on the tested example. A new image is needed for each release | Needs maintainer approval to publish |
| What was run on which route is in [installation](installation.md) | See there |
| A `git pull` in the clone reinstalls dndsR automatically on the next run, so the analysis version can change silently. Provenance records the version | Open; pin a tag for production analyses |
| `--threads` gives `append_annotations` no speedup, and `extract_cds` and the enrichment commands use at most 2 workers per comparison. See [thread scaling](performance.md) | Open |
| `calculate_dnds` needs about 1 GB per core. Called from an R session that already ran the workflow, 32 cores was killed at 96 GB with an unclear error | Open; use fewer cores or a fresh session |
| Windows runs sequentially and has not been tested | Open |
