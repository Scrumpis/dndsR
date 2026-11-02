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
| Default `go_enrichment` exceeds 5% under the null that keeps dN/dS linked to gene length and dS: 9.5% (7.0–12.8) of runs have a false term, on a real genome pair, and 1% under the other two nulls. One broad term, extracellular region, accounts for most of it, and 3.0% without it. See [enrichment calibration](enrichment.md#go-enrichment) | Open: modeling these confounders is on the roadmap to v1.0.0 | Do not read a significant large GO term as evidence about its function. Check its genes' lengths and dS |
| `adjust_scope = "per_type"` in `ipr_enrichment` and `"per_ontology"` in `go_enrichment` correct each table separately, with higher rates across tables, as does `term_enrichment` across its term types | Documented | Do not pool these tables as corrected together. The defaults of `ipr_enrichment` and `go_enrichment` span all their tables |
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
| `ipr_enrichment` looks up the current InterPro release even with `entries_source = "local"`, and with no network stops with a message saying local files would avoid it. No wrong results. Use `auto_detect_release = FALSE` | Open |
| Windows runs sequentially and has not been tested | Open |
