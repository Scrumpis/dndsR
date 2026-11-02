# Verification record

Only checks actually run are listed. Update these when the relevant software changes. Identify evidence by commit and content hash, never by file date.

## Environment

Linux amd64 HPC node, singularity-ce 4.1.2, no Docker. Image `sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831`: R 4.5.2, DIAMOND 2.1.14, OrthoFinder 2.5.4, orthologr 0.4.2, pwalign 1.6.0, topGO 2.62.0, IHW 1.38.0, qvalue 2.42.0, RIdeogram 0.2.2. dndsR itself is not in this image.

## Tested

| Check | Evidence | Result |
|---|---|---|
| Every command in the getting-started tutorial, in order, in a clean `HOME` | `tools/dev/test-tutorial.sh` | Pass |
| `dndsr verify_example` | Tutorial and tests | 14 of 14 |
| Package tests, run in the image | `tools/dev/run-tests.sh` | Pass, none skipped |
| `R CMD check --no-manual` | Run in the image with `lme4` missing | No errors or warnings |
| `extract_cds`, `calculate_dnds`, `append_annotations`, `term_enrichment`, `ipr_enrichment`, `dnds_summary`, `dnds_contrast`, `split_comparisons`, pairwise `dnds_state_contrast` | Results compared with independent values in tests | Pass |
| `go_enrichment`, `ipr_enrichment` (including parent-child), `dnds_ideogram` windows | Constructed cases with hypergeometric or hand-calculated answers | Pass |
| Rendered plots | Run on the example with the bundled Liberation Sans font | See `docs/validation/` |
| Public real-genome pair, *A. thaliana* and *A. lyrata* | `docs/validation/public-genome.md` | CDS match gffread for all compared genes; dN/dS agree with PAML and KaKs_Calculator; byte-identical across thread counts |
| Enrichment on real public annotations | `docs/validation/enrichment.md`, three null models | Empirical false-positive rates recorded; `p_adj` agrees with an independent Python recalculation |
| `dndsr ipr_enrichment`, default and `--adjust-scope per_type`; with `--fdr-method none`; with `--method parent_child`, both scopes | Run on the real public table, 8 threads; recalculated independently except parent-child | One `n_tested` across all tables by default, per type with the option, and `p_adj` = `p_value` with `none`; exit 0 |
| `dndsr go_enrichment`, default, `--adjust-scope per_ontology`, `--p-adjust none`, `--algorithm classic`, `--include_definition false` | Run on the real public table, 8 threads; BH recalculated independently, as are `classic` raw p-values | One `n_tested` across ontologies by default; `p_adj` = `raw_p` with `none`; exit 0 |
| GFF lines with `#` in attributes | `docs/validation/enrichment.md`, against source annotation files | 73 pairs recovered; no term lost; no other field changed |
| `dndsr covariate_enrichment`, GO and InterPro, at 1 and 2 threads | `tests/testthat/test-covariate-enrichment.R` with `pbinom`, `phyper`, and hand arithmetic; `covariate-verify.py` recalculates counts, Fisher's p-values, and BH without dndsR, on AtAl and SlSt tables; calibration at `docs/validation/enrichment-confounding.md` | Pass, no differences, identical across threads. Splines are not independently refit |
| `dndsr covariate_enrichment --dnds-annot-file ... --sides query --filter-expr "dS <= 3 * median(dS, na.rm = TRUE)"` | Run on the AtAl table | Runs, and calibrated on 2 pairs: worst null cell 5.0%. `docs/validation/enrichment-confounding.md` |
| Every command in `docs/tutorials/public-genomes.md`, in order, in a clean `HOME`, from commit `eb92852` | `tools/dev/test-public-tutorial.sb` | Pass: 30 min, 26 GB, results in the tutorial |
| `tools/check-first-exon-phase.py`, `tools/annotate-from-uniprot.py` | The same tutorial run | Pass |
| `extract_cds` stopping on a trans-spliced gene, with a hint to remove organelles | `tests/testthat/test-statistics.R`, and tomato's plastid *rps12* | Pass |
| Apptainer, conda, candidate image | `docs/validation/installation.md` | Pass |
| Skill discovery | Claude Code 2.1.316 and Codex CLI 0.165.0 in a clone | Found |

## Not verified

Docker; macOS; Windows; R from CRAN and Bioconductor on a working system; GitHub Actions; gene-level `dnds_state_contrast`; an agent completing the full workflow, unless recorded in `docs/dev/PROGRESS.md`.
