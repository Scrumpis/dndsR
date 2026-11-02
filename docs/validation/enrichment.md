# Enrichment calibration

Do the enrichment tests give false positives at the rate they claim? This page tests them on a real, publicly annotated genome pair. Everything here can be rebuilt with the scripts named below.

**Result.** Before v0.3.0, the defaults of `term_enrichment` and `ipr_enrichment` were badly miscalibrated: with dN/dS assigned to genes at random, where any significant term is false, they reported one in 59–74% and 95–100% of runs. The cause, a filter applied before multiple-testing correction, is corrected. `ipr_enrichment` and `go_enrichment` now also apply BH across one family by default. `term_enrichment` and `ipr_enrichment` then report a significant term in 1–3% of null runs. `go_enrichment` does in 1.0% and 0.8% under `permute` and `rotate`, **but in 9.5% (7.0–12.8) under `strata`, above nominal.** Most of that excess is one term, extracellular region. See [GO enrichment](#go-enrichment) for what this means. BH correction does not guarantee nominal error control under all conditions: it assumes p-values that are valid under the null, and enrichment p-values do not account for gene length or dS. Significant enrichment terms remain exploratory, not independently validated evidence of selection.

These results concern false positives under no association. They do not show that the remaining terms reflect selection: see [what selects the candidates](#what-selects-the-candidates).

## Data

| Item | Source | SHA-256 (first 16) |
|---|---|---|
| Genomes | *A. thaliana* (query, RefSeq TAIR10.1) and *A. lyrata* (subject), as in [public-genome.md](public-genome.md) | See that page |
| GO terms | GO Consortium TAIR GAF, without `NOT` or `ND` annotations | `582ea7a8cbbae942` |
| InterPro terms | UniProt reference proteome UP000006548 | `159a133c5b9fb88f` |
| Per-gene table, before | `AtAl_dnds_annot.tsv`, from `annotate-gff.py`, then the standard workflow | `603fdcd5d8f67ab7` |
| Per-gene table, after | The same workflow rerun with the corrected code | `56321306d1b52a23` |
| InterPro tree, parent-child only | `ParentChildTreeFile.txt`, current release when downloaded | `3a40c7856308a0a8` |

Code: original is commit `59faf30`; the first correction is `1a363fb`; final defaults are the commit recorded in [RELEASE-READINESS.md](../dev/RELEASE-READINESS.md). Image `dndsr_v0.2.0_sha256-726657e9`. Seeds are 20261101 to 20261116 for real-function runs (`--adjust-scope per_type`: the first 8), 20262101 to 20262110 for parent-child (`per_type`: the first 5), 20263101 to 20263110 for `1a363fb`'s global runs, 20261101 for 5,000 draws, and 1001 to 1020 for GO. Earlier evidence is preserved in [enrichment-59faf30.md](enrichment-59faf30.md) and [enrichment-1a363fb.md](enrichment-1a363fb.md), with their data files beside the current ones.

The workflow ran end to end. There are 22,146 pairs, of which 22,091 pass `max_dnds`, and 206 have dN/dS > 1. See [GFF recovery](#gff-recovery) for how appended terms were checked.

## Null models

The dN/dS column is reassigned, so no term can be truly associated with it:

| Null | What it keeps | Why |
|---|---|---|
| `permute` | Nothing | The usual label permutation. It ignores all structure |
| `strata` | Gene length, dS, and annotation richness (permuted within strata) | Holds the confounders found below constant |
| `rotate` | Genomic order (circular shift) | Keeps tandem clusters and local dependence together |

## Term and InterPro enrichment

**The defect.** `min_pos` removed terms with few positive genes before BH. That count is the test's outcome, so only terms that could be significant remained, and BH corrected for about 100 tests instead of about 7,200. Among all terms, 2% have p < 0.05, because Fisher's test is conservative on small counts. Among those kept, 26–34% do. In parent-child mode, the same filter also removed terms before the child-to-parent pass.

**The correction.** BH now covers every term meeting `min_total` and `max_prop`, which do not depend on the outcome. `min_pos` only chooses which rows are written. `n_tested` gives the size of the family. Terms that have no positive gene still count in the family, including an InterPro type with no positive gene.

**The final default.** `ipr_enrichment` now corrects once across all eligible InterPro entry types (`adjust_scope = "global"`). Before, each type had its own family, and a run had up to seven chances, so a list drawn from several tables was not corrected across them. `adjust_scope = "per_type"` (`--adjust-scope per_type`) keeps the old behavior. This is an intentional change to a default. In parent-child mode, the child-to-parent decisions are still made within each entry type, and the final `p_adj` is corrected across every type with parent-child edges.

**False positives.** The real functions, except that `term_enrichment` is given `terms = "ipr", exclude_terms = NULL`, because by default it leaves InterPro and GO to their own commands. On null assignments ([data](data/enrichment_real_null.tsv)). The figure is the share of runs with any significant term, with its 95% interval.

| Function | Null | Runs: original / final | Original (`59faf30`) | Final default | Mean false terms, original → final |
|---|---|---|---|---|---|
| `term_enrichment` | `permute` | 104 / 208 | 58.7% | **1.0%** (0.3–3.4) | 4.96 → 0.01 |
| | `strata` | 104 / 208 | 94.2% | **2.9%** (1.3–6.1) | 8.50 → 0.03 |
| | `rotate` | 104 / 208 | 74.0% | **1.9%** (0.8–4.8) | 6.18 → 0.04 |
| `ipr_enrichment` | `permute` | 104 / 208 | 95.2% | **1.0%** (0.3–3.4) | 6.79 → 0.01 |
| | `strata` | 104 / 208 | 100% | **2.9%** (1.3–6.1) | 9.16 → 0.03 |
| | `rotate` | 104 / 208 | 98.1% | **1.9%** (0.8–4.8) | 7.31 → 0.04 |
| `ipr_enrichment`, parent-child | `permute` | 50 / 100 | 72.0% | **0.0%** (0.0–3.7) | 3.76 → 0.00 |
| | `strata` | 47 / 90 | 80.9% | **2.2%** (0.6–7.7) | 3.74 → 0.03 |
| | `rotate` | 50 / 100 | 78.0% | **0.0%** (0.0–3.7) | 3.64 → 0.00 |

The selectable per-type correction, from [`1a363fb`'s runs](data/enrichment_real_null_1a363fb.tsv), which final code reproduces exactly on the seeds rerun (8 and 5 tasks): `ipr_enrichment` 8.7% (5.5–13.3), 16.8% (12.4–22.5), and 8.7% (5.5–13.3); parent-child 1.0% (0.2–5.4), 8.0% (4.1–15.0), and 0.0% (0.0–3.7), in the same order of nulls. Those are rates across all tables. Each table alone was at most 5.3%. An independent set of 100 runs per null with global correction, on other seeds at `1a363fb`, gave 1.0%, 6.0%, and 0.0%.

Over the final default's 208 runs per null, no single `ipr_enrichment` table had a false term in more than 1.9% of runs. `term_enrichment` and `ipr_enrichment` have the same rates, on different families (7,208 and 5,202 terms), because the same draws produced nearly the same false terms: 5 of 624 runs differ in their count, and none changes whether the run has one.

**Which error rate this is.** Under these nulls no association exists, so every significant term is false. For one run, any significant term makes the false discovery proportion 1, and none makes it 0. So the rate above is at once the family-wise error rate and the false discovery rate under the complete null, and BH at 0.05 aims to keep it at or below 5%. *It does not measure false discovery when real associations exist*, because none of these nulls contains one. BH does not control the family-wise rate then, only the false discovery rate, and only under assumptions on dependence that these terms need not meet. Results are also conditional on one genome pair.

The two parent-child `strata` figures exclude 10 runs of 100 that were stopped by a validation accident, not by dndsR: a concurrent job replaced the output directory, and each task's run failed with `cannot open the connection`. They are recorded and left out, and the original's 3 of 50 are left out in the same way, for its own reason, no parent-child edges.
**Corroboration.** A reimplementation, verified identical to the old `term_enrichment` on the 67 observed terms (largest p-value difference under 1e-10), allowed 5,000 draws per null on the corrected table ([data](data/enrichment_null_calibration.tsv)). "BH on all terms with ≥ 2 genes" is what `term_enrichment` now does:

| Null | Old default: filter, then BH | BH on all terms with ≥ 2 genes | ≥ 5 genes | ≥ 10 genes |
|---|---|---|---|---|
| `permute` | 60.9% | 0.7% | 1.3% | 1.2% |
| `strata` | 90.6% | 3.4% | 5.6% | 5.5% |
| `rotate` | 66.3% | 0.9% | 2.2% | 1.5% |

**Independent check.** `enrichment-verify.py` recomputes Fisher's test and BH in plain Python, with none of dndsR's code, from the annotated table. On the real data it agrees for `ipr_enrichment` with global, per-type, and no correction on which terms are reported, on `n_tested`, and on `p_value` and `p_adj` (largest relative difference 4.5e-11). It agreed for `term_enrichment` at `1a363fb` (5.7e-11), whose code has not changed. Parent-child's conditional values are covered by hand-calculated values in `tests/testthat/test-enrichment-family.R`, not by this script. The package's tests also compare against `phyper` and a BH calculated by hand, and show that shuffling the input's row order changes nothing.

### Real data, before and after

Significant terms at `p_adj` ≤ 0.05 ([data](data/enrichment_before_after.tsv)). "Per type" is `1a363fb`'s default. No term is significant after that was not significant before. The same terms are significant on the old table too, so the recovered annotations did not change the first correction's list.

| Function, table | Reported terms | Original | Per type (`1a363fb`) | **Final default** | Final family size |
|---|---|---|---|---|---|
| `term_enrichment`, all InterPro | 67 | 48 | 18 | **18** | 7,208 |
| `ipr_enrichment`, Domain | 21 | 16 | 4 | **4** | 5,202 |
| `ipr_enrichment`, Family | 10 | 10 | 7 | **5** | 5,202 |
| `ipr_enrichment`, Homologous superfamily | 15 | 9 | 6 | **4** | 5,202 |
| `ipr_enrichment`, Repeat | 4 | 3 | 3 | **2** | 5,202 |
| `ipr_enrichment`, all four | 50 | 38 | 20 | **15** | 5,202 |
| `ipr_enrichment`, parent-child, Domain | 21 (19 before) | 14 | 3 | **3** | 3,880 |
| `ipr_enrichment`, parent-child, Family | 10 (none before) | not tested | 7 | **5** | 3,880 |

Before the correction, parent-child skipped Family, because the filter left it no parent-child edges.

The 18 `term_enrichment` survivors, with Westfall-Young adjusted p-values from each null as a stricter check that respects dependence between terms:

| Term | Name | Positive of total | `p_adj` before → after | WY `permute` / `strata` / `rotate` |
|---|---|---|---|---|
| IPR011043 | Galactose oxidase/kelch, beta-propeller | 13 of 160 | 3.26e-09 → 3.5e-07 | 0.0002 / 0.0002 / 0.0002 |
| IPR001810 | F-box domain | 16 of 484 | 4.44e-06 → 0.000476 | 0.0006 / 0.0006 / 0.001 |
| IPR036047 | F-box-like domain superfamily | 15 of 478 | 1.45e-05 → 0.00155 | 0.0006 / 0.006 / 0.0012 |
| IPR050354 | (not in bundled list) | 6 of 51 | 1.76e-05 → 0.00189 | 0.0018 / 0.0068 / 0.0016 |
| IPR057499 | (not in bundled list) | 6 of 53 | 1.77e-05 → 0.0019 | 0.0018 / 0.0104 / 0.0018 |
| IPR017451 | F-box associated interaction domain | 9 of 167 | 1.9e-05 → 0.00203 | 0.0018 / 0.0108 / 0.002 |
| IPR050796 | (not in bundled list) | 7 of 95 | 3.06e-05 → 0.00329 | 0.0032 / 0.018 / 0.0036 |
| IPR034086 | Pectinesterase inhibitor, plant | 4 of 18 | 4.35e-05 → 0.00467 | 0.005 / 0.0244 / 0.0058 |
| IPR011498 | Kelch repeat type 2 | 4 of 21 | 7.44e-05 → 0.008 | 0.0102 / 0.0448 / 0.014 |
| IPR052421 | (not in bundled list) | 4 of 24 | 0.000107 → 0.0115 | 0.015 / 0.0624 / 0.0206 |
| IPR010851 | Defensin-like protein | 5 of 49 | 0.000107 → 0.0115 | 0.015 / 0.0652 / 0.0206 |
| IPR006527 | F-box associated domain, type 1 | 6 of 87 | 0.000135 → 0.0145 | 0.017 / 0.0748 / 0.0218 |
| IPR022618 | Defensin-like protein 20-27 | 3 of 10 | 0.000169 → 0.0181 | 0.025 / 0.0824 / 0.0282 |
| IPR008801 | Rapid ALkalinization Factor | 4 of 30 | 0.000209 → 0.0225 | 0.0302 / 0.189 / 0.0368 |
| IPR016140 | Bifunctional inhibitor/plant lipid transfer protein/seed storage helical domain | 6 of 99 | 0.000211 → 0.0227 | 0.0358 / 0.196 / 0.0436 |
| IPR036312 | Bifunctional inhibitor/plant lipid transfer protein/seed storage helical domain superfamily | 6 of 99 | 0.000211 → 0.0227 | 0.0358 / 0.196 / 0.0436 |
| IPR015915 | Kelch-type beta propeller | 6 of 110 | 0.000358 → 0.0384 | 0.135 / 0.234 / 0.161 |
| IPR005508 | Protein of unknown function DUF313 | 3 of 14 | 0.000362 → 0.0389 | 0.152 / 0.246 / 0.18 |

Nine are also below 0.05 by Westfall-Young under all three nulls, and these are the most robust: Kelch and F-box entries, a pectinesterase inhibitor family, and three InterPro entries not in the bundled names list. Seven more pass under `permute` and `rotate` but not `strata`, so the gene properties held constant by `strata` could explain them. Tandem clusters do not explain the original list: collapsing them kept 45 of the 48. **Whether these terms make biological sense has not been compared with the literature.** Several names suggest small secreted peptides, such as defensin-like and rapid alkalinization factor, which would fit the [length bias](#what-selects-the-candidates), but the gene lengths behind each term have not been tested.

## GO enrichment

**The default.** `go_enrichment` now applies BH by default, in `p_adj`, and `raw_p` stays as topGO gives it. Before, `p_adj` by default was a copy of `raw_p`. One family covers every tested term in all requested ontologies, separately for each comparison and side, so each table's `n_tested` is the same and the tables can be pooled. `adjust_scope = "per_ontology"` gives each ontology its own family. `p_adjust = "none"` copies `raw_p` to `p_adj`, labelled `none` in `p_adjust_method` and `adjust_scope`. An ontology with no positive gene writes no table, and its terms stay in the family.

**Why.** Thousands of terms are tested, and there was no compelling reason not to correct. topGO's authors do not advocate against adjustment. They note that adjusted values can be misleading and very conservative because topology-aware tests are conditional, and add that this "is not necessary a bad thing" (topGO 2.62.0 vignette). BH's guarantee is not established for `weight01`'s conditional, dependent p-values, so the case rests on the measurements below: BH is applied to whatever p-values the chosen `algorithm` gives, and the same applies to `elim`, `weight`, and `classic`. `go_enrichment` has no parent-child mode. `go_enrichment` never filtered by `min_pos`.

**Independent check.** `enrichment-verify.py go` recalculates BH from `raw_p` for global, per-ontology, and unadjusted tables, largest relative difference 5.2e-15. For `classic`, whose tests are not conditional, it recomputes each term's genes and raw p-value from the annotated table and GO.db's ancestor table, without topGO. All terms and counts agree, and p-values agree to 8.2e-11. `weight01`'s raw p-values have no such independent recalculation. Tests in `tests/testthat/test-enrichment-family.R` compare BH with hand values for all four algorithms.

**False positives.** The real `go_enrichment` on null assignments, 400 runs per null ([data](data/enrichment_go_null.tsv), [by ontology](data/enrichment_go_null_by_ontology.tsv)). As above, the rate is both the family-wise rate and the false discovery rate under the complete null.

| Algorithm, runs | Null | Mean terms at raw p < 0.05 | **Default, global** | Mean false terms | `per_ontology` | Default, without GO:0005576 |
|---|---|---|---|---|---|---|
| `weight01`, 400 | `permute` | 21.0 | **1.0%** (0.4–2.5) | 0.010 | 2.8% (1.5–4.9) | 1.0% (0.4–2.5) |
| | `rotate` | 21.1 | **0.8%** (0.3–2.2) | 0.007 | 2.2% (1.2–4.2) | 0.8% (0.3–2.2) |
| | `strata` | 11.6 | **9.5%** (7.0–12.8) | 0.100 | 20.5% (16.8–24.7) | 3.0% (1.7–5.2) |
| `classic`, 200 | `permute` | 55.7 | 1.0% (0.3–3.6) | 0.010 | 5.0% (2.7–9.0) | 1.0% (0.3–3.6) |
| | `rotate` | 59.5 | 0.0% (0.0–1.9) | 0.000 | 3.5% (1.7–7.0) | 0.0% (0.0–1.9) |
| | `strata` | 25.9 | 10.5% (7.0–15.5) | 0.175 | 20.5% (15.5–26.6) | 4.0% (2.0–7.7) |

**The remaining problem.** Under `strata`, the default's interval is above 5%. This is not a defect in the correction, which reaches about 1% where nothing is held linked to dN/dS. BH assumes p-values are valid under the null, and this term's are not. GO:0005576, extracellular region, is significant in 28 of the 38 `strata` runs that have a false term, and has raw p < 0.05 in 311 of 400 (78%). `strata` shuffles dN/dS only among genes of similar length, dS, and GO richness, so genes in the term keep more than their share of positive labels. Fisher's test takes those labels to be exchangeable across all genes, so its p-value is too small for this term whichever correction follows. It is a confounded term, in a default that does not model confounding. The final rate depends on a term's size too: this term is 14% of the universe, so a small shift is highly significant.

How much of the real signal do length and dS reproduce? [`enrichment-go-confounding.py`](data/enrichment_go_confounding.tsv) shuffles within finer strata and counts positive genes from each term, with `classic` p-values, 300 draws. At random these terms would hold 14.8, 0.3, 0.1, and 0.2 positive genes:

| Term (members; median length, others ~400 aa) | Observed positive | 5 × 5 bins | 10 × 10 bins | Draws as extreme as observed, 5 × 5 / 10 × 10 |
|---|---|---|---|---|
| GO:0005576 extracellular region (2,104; 346 aa) | 44 | 24.5 | 31.3 | 0 / 0 |
| GO:0007267 cell-cell signaling (77; 101 aa) | 7 | 1.5 | 2.3 | 0 / 0 |
| GO:0005179 hormone activity (46; 97 aa) | 4 | 1.0 | 1.4 | 1 / 1 |
| GO:0019722 calcium-mediated signaling (66; 152 aa) | 4 | 1.1 | 1.1 | 1 / 1 |

Holding length and dS within 10 bins each brings extracellular region to 31 of its 44 positive genes, where chance alone would give 15. So most of its excess can arise from those properties, though not all of it. The figures are descriptive: finer bins leave fewer genes to exchange, about 74 per stratum at 10 bins and 5 at 40 (in the data file), so at the extreme the shuffle returns the observed labels, and they are not error rates. Because these are `classic` counts, they differ from the `weight01` p-values below. No correction for this is made in v0.3.0, no term is removed, and no threshold is changed: the maintainer accepted this as a documented limitation, because a correction needs confounder modeling and biological validation, which are on the [v1.0.0 roadmap](../dev/ROADMAP-1.0.md).

**Real data.** Significant terms at `p_adj` ≤ 0.05 among 2,978 tested. `p_adj` below 0.05 by the old default meant raw p below 0.05.

| Ontology (terms) | Old default, `p_adj` = `raw_p` | **Final default** | `per_ontology` | `classic`, global |
|---|---|---|---|---|
| BP (1,963) | 19 | **2** | 1 | 2 |
| MF (684) | 9 | **1** | 1 | 6 |
| CC (331) | 3 | **1** | 1 | 1 |
| All | 31 | **4** | 3 | 9 |

| Term | Ontology | Positive of annotated | `raw_p` | `p_adj` | Smallest `raw_p` of any term in that ontology, over null runs: `permute` / `rotate` / `strata` |
|---|---|---|---|---|---|
| GO:0005576 extracellular region | CC | 44 of 2,104 | 2.87e-10 | 8.56e-07 | 3.65e-5 / 2.16e-5 / 9.08e-9 |
| GO:0007267 cell-cell signaling | BP | 7 of 77 | 9.80e-08 | 1.46e-04 | 5.73e-6 / 1.06e-5 / 2.4e-7 |
| GO:0005179 hormone activity | MF | 4 of 46 | 7.26e-06 | 7.20e-03 | **5.64e-6** / 6.35e-5 / **4.12e-7** |
| GO:0019722 calcium-mediated signaling | BP | 4 of 66 | 6.57e-05 | 4.89e-02 | **5.73e-6** / **1.06e-5** / **2.4e-7** |

Bold marks a null run in which some term went beyond this term's observed value. Extracellular region should not be read as evidence about extracellular proteins, as the table of finer strata shows. Hormone activity and calcium-mediated signaling are not reliably beyond chance. All are leads, as in [interpretation](../interpretation.md).

## GFF recovery

**The defect.** `append_annotations` and `dnds_ideogram` read GFF files with `#` as a comment character (`data.table::fread(comment = "#")`, which matches `comment.char` by partial name, or `read.table`), and that cuts a line at any `#`. GFF3 comments start at the beginning of a line, and `#` can appear in attribute values. 79 lines in the *A. thaliana* GFF, 78 of them genes, contain one; the *A. lyrata* GFF has none. Attributes after it were lost.

**The correction.** A shared reader treats a line as a comment only if `#` comes first, stops at `##FASTA`, and rejects lines that do not have nine columns. Percent-encoded characters defined by GFF3 are now decoded when terms are extracted, with `%25` last, so that `%2523` becomes `%23`.

**Verification** by `gff-recovery.py`, comparing every field of every row before and after, then checking terms against source files that dndsR did not produce:

| Check | Result |
|---|---|
| Rows, columns, and order | Unchanged, 22,146 rows |
| Fields changed | `q_attributes`, `q_go`, and `q_ipr` only, on 73 rows |
| Pairs that lost terms | 73: GO on 71, InterPro on 71, and both on 69. Earlier reports of "71 pairs" counted one column at a time |
| Terms recovered | 682 GO, 168 InterPro |
| Those 73 pairs against GO source | All identical, 0 GO terms not in the source |
| Those 73 pairs against InterPro source | All contain the source's terms. 41 further IDs are absent from it, and all 41 are on that gene's own GFF lines, in RefSeq `Note` attributes |
| Terms lost by any pair | 0 |
| All pairs in the GO source | 20,746 identical to it, up from 20,675. The other 1,400 pairs are absent from it and have no GO terms |

Unrelated outputs are byte-identical before and after the default changes: `AtAl_dnds.tsv`, `AtAl_dnds_annot.tsv`, `AtAl_dnds_summary.tsv`, both ideogram SVGs, and both CDS FASTA files. On the final code, the annotated table and every `ipr_enrichment` and `go_enrichment` table are byte-identical at 1 and 8 threads.

## What selects the candidates

The positive genes are not a random sample ([data](data/enrichment_confounders.tsv)):

- **Length.** Positive rate falls from 3.05% in the shortest fifth of genes to 0.09% in the longest. Positive genes have a median of 125.5 amino acids, and others 364.
- **dS.** It falls from 3.10% in the lowest fifth to 0.16% in the highest: dN/dS is a ratio, and is unstable when dS is small.
- **Annotation.** 0.66% of annotated genes are positive, and 3.98% of unannotated genes. So 72 of 206 positive genes (35%) cannot contribute to any enrichment test.

Under `strata`, which holds these constant within coarse strata, the corrected tests stay near 5% for InterPro terms, and not for one broad GO term ([above](#go-enrichment)). That does not rule out confounding in the real labels: a significant term can still reflect short genes, so report term results with their gene lengths in mind.

The result also depends on the candidate rule ([data](data/enrichment_candidate_rules.tsv)):

| Rule | Significant terms (old procedure) | Shared with dN/dS > 1 | Jaccard |
|---|---|---|---|
| dN/dS > 1 (default) | 48 | 48 | 1 |
| dN/dS > 0.5 | 139 | 37 | 0.25 |
| dN/dS > 0.75 | 86 | 42 | 0.46 |
| dN/dS > 1.5 | 3 | 3 | 0.06 |
| dN/dS > 1 and ≥ 150 amino acids | 40 | 35 | 0.66 |
| dN/dS > 1 and dS ≥ 0.05 | 45 | 43 | 0.86 |

Length matters more than a dS floor here. These figures use the old procedure, filter then BH, as in the first column of the calibration table, so read them only for how much the list moves between rules.

## Defects found

| Defect | Effect | Status |
|---|---|---|
| `min_pos` filtered before BH in `term_enrichment` and `ipr_enrichment` | Severe, above | Corrected at `1a363fb` |
| The GFF reader treated `#` inside attributes as a comment | 79 lines truncated, and 73 pairs silently lost terms | Corrected at `1a363fb` |
| `go_enrichment(include_definition = FALSE)` always stopped with an error | Explicit failure, no wrong results. Now removes `term_def`; every other field on the real data is identical | Corrected at `1a363fb` |
| `go_enrichment` stopped when an ontology had no positive genes ("allGenes must be a factor with 2 levels") | Explicit failure, no wrong results. It now writes no table for that ontology, and its terms stay in the family | Corrected with the final defaults |
| `p_adj` in `go_enrichment` was uncorrected by default; `ipr_enrichment` was corrected only within each type | Not defects, but defaults that overstated significance across tables. Maintainer chose to change both | Changed with the final defaults |
| `ipr_enrichment` looks up the current InterPro release even with `entries_source = "local"`, and stops if it cannot connect, with a message that local files would avoid it | No wrong results. Set `auto_detect_release = FALSE` | Open: [known issues](known-issues.md) |

## Limits

One genome pair, at one divergence. The correction controls false positives under no association; it does nothing for confounding when an association is real, and it does not stop a term confounded with length or dS from passing under `strata`, as GO:0005576 shows. `permute` and `rotate` do not hold confounders constant, and `strata` uses coarse strata and leaves dN/dS linked to length and dS within them. Nothing here measures false discovery with real associations, or power. No result has been compared with independent biological evidence. Modeling confounders and biological benchmarking stay on the [v1.0.0 roadmap](../dev/ROADMAP-1.0.md). The rank test in `enrichment-calibration.R` was run but is not interpreted here. `enrichment-verify.py` does not cover parent-child, whose expected values come from the package's constructed tests. It does not recompute `weight01`'s raw p-values, only `classic`'s, and BH on `weight01`'s. `term_enrichment` writes one table for each term type, each with its own BH family, as `ipr_enrichment`'s `per_type` does, so do not pool terms across its tables as corrected together. Old-code runs were fewer than new-code runs.

## Rebuild

Run in the image, with `dndsR` installed, after building the table with `annotate-gff.py` and the standard workflow. `enrichment-calibration.R` checks its reimplementation against a `term_enrichment` table from the old code, so give it one from `59faf30`.

```
Rscript tools/dev/validation/enrichment-real-null.R <dnds_annot.tsv> <runs per null> <seed> <outdir> [InterPro tree]
IPR_ADJUST_SCOPE=per_type Rscript tools/dev/validation/enrichment-real-null.R <dnds_annot.tsv> <runs per null> <seed> <outdir> [InterPro tree]
python3 tools/dev/validation/enrichment-null-summary.py <label>=<outdir> ...
Rscript tools/dev/validation/enrichment-calibration.R <dnds_annot.tsv> <old IPR enrichment.tsv> <term column> 5000 20261101 <outdir>
GO_ALGORITHM=classic Rscript tools/dev/validation/enrichment-go-null.R <dnds_annot.tsv> <task> <runs per null> <outdir>   # omit GO_ALGORITHM for the default
python3 tools/dev/validation/enrichment-go-null-summary.py GO:0005576 <rates.tsv> <by_ontology.tsv> <label>=<outdir> ...
Rscript tools/dev/validation/enrichment-go-ancestors.R <ancestors.tsv>
python3 tools/dev/validation/enrichment-verify.py go <dnds_annot.tsv> <ancestors.tsv> <dir> <prefix> <algorithm> [node_min]
python3 tools/dev/validation/enrichment-go-confounding.py <dnds_annot.tsv> <ancestors.tsv> <term> 1.68e-5 300 20261102
python3 tools/dev/validation/enrichment-verify.py term <dnds_annot.tsv> <table> [min_total] [min_pos]
python3 tools/dev/validation/enrichment-verify.py ipr <dnds_annot.tsv> <entry list> <dir> <prefix> [min_total] [min_pos] [max_prop]
python3 tools/dev/validation/gff-recovery.py <old dnds_annot.tsv> <new dnds_annot.tsv> <truth.tsv> <query GFF>
python3 tools/dev/validation/enrichment-before-after.py <old dir> <new dir> <observed_terms.tsv> <entry list>
```
