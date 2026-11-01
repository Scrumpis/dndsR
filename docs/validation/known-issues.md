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
| `term_enrichment` and `ipr_enrichment` filter by `min_pos` before multiple-testing correction. Measured on a real genome pair with no true association: 58% and 95% of runs report a significant term, where 5% is nominal. See [enrichment calibration](enrichment.md) | **Open, needs a decision** | Treat `p_adj` as a ranking, or recalculate BH across all terms of adequate size |
| The GFF reader truncates attributes at a `#` inside a line, which is legal GFF3. In *A. thaliana* RefSeq, 71 pairs silently lose terms | Open | Check for `#` in attributes before annotating |
| Most of the default `go_enrichment` list at p < 0.05 cannot be told apart from noise. See [enrichment calibration](enrichment.md) | Documented | Use a stricter threshold, and validate terms independently |
| A dN/dS has no uncertainty attached, and `Comeron95` truncates undefined terms without flagging them | Open; see `docs/dev/ROADMAP-1.0.md` | See [interpretation](../interpretation.md) |
| Fixed 1.96 in normal confidence intervals; dependence between genes; no adjustment across contrasts | Documented | See [interpretation](../interpretation.md) |

Interpretation limits and the estimator's bias are covered in [interpretation](../interpretation.md) and [estimator validation](estimators.md).

## Software

| Issue | Status |
|---|---|
| The published image, `v0.2.0`, predates this code. It works with the current code on the tested example. A new image is needed for each release | Needs maintainer approval to publish |
| `go_enrichment(include_definition = FALSE)` always stops with an error | Open |
| What was run on which route is in [installation](installation.md) | See there |
| A `git pull` in the clone reinstalls dndsR automatically on the next run, so the analysis version can change silently. Provenance records the version | Open; pin a tag for production analyses |
| Windows runs sequentially and has not been tested | Open |
