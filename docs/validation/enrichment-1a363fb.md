# Enrichment calibration

> **Superseded.** This page records the evidence at `1a363fb`, where `ipr_enrichment` defaulted to per-type correction and `go_enrichment` to no correction. Current results: [enrichment.md](enrichment.md).

Do the enrichment tests give false positives at the rate they claim? This page tests them on a real, publicly annotated genome pair. Everything here can be rebuilt with the scripts named below.

**Result.** Before v0.3.0, the defaults of `term_enrichment` and `ipr_enrichment` were badly miscalibrated: with dN/dS assigned to genes at random, where any significant term is false, they reported one in 59–74% and 95–100% of runs. The cause was a filter applied before multiple-testing correction, and it is corrected in this release. `term_enrichment` now reports one in 1–3% of runs. Each `ipr_enrichment` table is at or near 5%, but its default corrects each InterPro entry type separately, so 9–17% of runs have a significant term in *some* table; `adjust_scope = "global"` brings that to 0–6%. `go_enrichment` was not affected by this defect, and its default list at p < 0.05 still cannot be told apart from noise by its size.

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

Code: before is commit `59faf30`, and after is `1a363fb`. Image `dndsr_v0.2.0_sha256-726657e9`. Seeds are 20261101 to 20261116 for real-function runs, 20262101 to 20262110 for parent-child, 20263101 to 20263110 for global correction, and 20261101 for 5,000 draws.

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

**The correction.** BH now covers every term meeting `min_total` and `max_prop`, which do not depend on the outcome. `min_pos` only chooses which rows are written. `n_tested` gives the size of the family. The interface is unchanged. Terms that have no positive gene still count in the family, including an InterPro type with no positive gene under `adjust_scope = "global"`.

**False positives.** The real functions at their defaults, except that `term_enrichment` is given `terms = "ipr", exclude_terms = NULL`, because by default it leaves InterPro and GO to their own commands. On null assignments, before (`59faf30`) and after (`1a363fb`) ([data](data/enrichment_real_null_1a363fb.tsv)). The figure is the share of runs with any significant term, with its 95% interval after correction. With no association, about 5% is nominal for a single BH family.

| Function | Null | Runs before / after | Before | After | Mean false terms, before → after |
|---|---|---|---|---|---|
| `term_enrichment` | `permute` | 104 / 208 | 58.7% | **1.0%** (0.3–3.4) | 4.96 → 0.01 |
| | `strata` | 104 / 208 | 94.2% | **2.9%** (1.3–6.1) | 8.50 → 0.03 |
| | `rotate` | 104 / 208 | 74.0% | **1.9%** (0.8–4.8) | 6.18 → 0.04 |
| `ipr_enrichment`, default | `permute` | 104 / 208 | 95.2% | **8.7%** (5.5–13.3) | 6.79 → 0.12 |
| | `strata` | 104 / 208 | 100% | **16.8%** (12.4–22.5) | 9.16 → 0.20 |
| | `rotate` | 104 / 208 | 98.1% | **8.7%** (5.5–13.3) | 7.31 → 0.12 |
| `ipr_enrichment`, `adjust_scope = "global"` | `permute` | not run / 100 | | **1.0%** (0.2–5.4) | 0.03 after |
| | `strata` | not run / 100 | | **6.0%** (2.8–12.5) | 0.08 after |
| | `rotate` | not run / 100 | | **0.0%** (0.0–3.7) | 0.00 after |
| `ipr_enrichment`, parent-child | `permute` | 50 / 100 | 72.0% | **1.0%** (0.2–5.4) | 3.76 → 0.01 |
| | `strata` | 47 / 100 | 80.9% | **8.0%** (4.1–15.0) | 3.74 → 0.11 |
| | `rotate` | 50 / 100 | 78.0% | **0.0%** (0.0–3.7) | 3.64 → 0.00 |

**Why default `ipr_enrichment` stays above 5%.** Its default, `adjust_scope = "per_type"`, gives each of up to seven InterPro entry types its own BH family, so a run has seven chances. Each table alone is near its nominal rate: over 208 runs per null, the highest rate is 5.3% (`strata`, Conserved site, 11 of 208), and all others are at most 3.8%. `p_adj` therefore means "corrected within this entry type". A list drawn from all the tables together has not been corrected across them. For one family across all types, use `adjust_scope = "global"` (`--adjust-scope global`), which brings the rate across all tables to 0–6%, and gives 5,202 tests on this data. The price on the real data is fewer terms: 15 significant, against 20 per type. Parent-child has two tables with edges in this data (Domain and Family) and does not support `global`. Under `strata`, its Domain and Family tables reached 3% and 5%.

**Corroboration.** A reimplementation, verified identical to the old `term_enrichment` on the 67 observed terms (largest p-value difference under 1e-10), allowed 5,000 draws per null on the corrected table ([data](data/enrichment_null_calibration.tsv)). "BH on all terms with ≥ 2 genes" is what `term_enrichment` now does:

| Null | Old default: filter, then BH | BH on all terms with ≥ 2 genes | ≥ 5 genes | ≥ 10 genes |
|---|---|---|---|---|
| `permute` | 60.9% | 0.7% | 1.3% | 1.2% |
| `strata` | 90.6% | 3.4% | 5.6% | 5.5% |
| `rotate` | 66.3% | 0.9% | 2.2% | 1.5% |

**Independent check.** `enrichment-verify.py` recomputes Fisher's test and BH in plain Python, with none of dndsR's code, from the annotated table. For `term_enrichment` and each of the four `ipr_enrichment` tables on the real data, it agrees on which terms are reported, on `n_tested`, and on `p_value` and `p_adj` (largest relative difference 5.7e-11). The package's tests also compare against `phyper` and a BH calculated by hand, and show that shuffling the input's row order changes nothing.

### Real data, before and after

Significant terms at `p_adj` ≤ 0.05 ([data](data/enrichment_before_after_1a363fb.tsv)). "Corrected" is the new code on the new table. No term is significant after that was not significant before. The same terms are significant on the old table too, so the recovered annotations did not change this list.

| Function, table | Reported terms | Before | Corrected | Family size |
|---|---|---|---|---|
| `term_enrichment`, all InterPro | 67 | 48 | **18** | 7,208 |
| `ipr_enrichment`, Domain | 21 | 16 | **4** | 1,972 |
| `ipr_enrichment`, Family | 10 | 10 | **7** | 1,908 |
| `ipr_enrichment`, Homologous superfamily | 15 | 9 | **6** | 913 |
| `ipr_enrichment`, Repeat | 4 | 3 | **3** | 45 |
| `ipr_enrichment`, `adjust_scope = "global"`, all four tables | 50 | not run | **15** | 5,202 |
| `ipr_enrichment`, parent-child, Domain | 21 (19 before) | 14 | **3** | 1,972 |
| `ipr_enrichment`, parent-child, Family | 10 (none before) | not tested | **7** | 1,908 |

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

The real `go_enrichment` at its defaults, on the corrected table, 200 runs per null ([data](data/enrichment_go_null_1a363fb.tsv)). It reports unadjusted `weight01` p-values, and never filtered by `min_pos`, so this procedure did not change. The recovered annotations add 7 tested BP terms and 1 MF term. Observed counts did not change, and the terms at p < 0.001 are the same.

| Ontology (terms) | Threshold | Observed | `permute` mean | `rotate` mean | `strata` mean |
|---|---|---|---|---|---|
| BP (1,963) | p < 0.05 | 19 | 12.4 | 12.2 | 5.2 |
| | p < 0.001 | 6 | 0.32 | 0.20 | 0.14 |
| MF (684) | p < 0.05 | 9 | 5.7 | 5.8 | 4.0 |
| | p < 0.001 | 4 | 0.17 | 0.13 | 0.17 |
| CC (331) | p < 0.05 | 3 | 3.0 | 3.2 | 2.6 |
| | p < 0.001 | 1 | 0.07 | 0.07 | 0.31 |

Observed counts at p < 0.05 are mostly inside the null's range, so most of that list could be noise. At p < 0.001, BP and MF exceed the nulls by more than ten-fold, and CC exceeds `strata` only three-fold. `weight01` p-values are conservative (well under 5% of terms fall below 0.05), so they are not ordinary p-values either way.

The strongest single term is not reliably beyond chance:

| Ontology | Strongest term | Observed p | Smallest null p: `permute` / `rotate` / `strata` |
|---|---|---|---|
| BP | cell-cell signaling (GO:0007267), 7 of 77 | 9.8e-8 | 5.7e-6 / 1.2e-5 / **5.6e-8** |
| MF | hormone activity (GO:0005179), 4 of 46 | 7.3e-6 | **5.6e-6** / 3.0e-5 / **4.3e-6** |
| CC | extracellular region (GO:0005576), 44 of 2,104 | 2.9e-10 | 3.7e-5 / 8.8e-6 / 5.4e-7 |

Bold marks a null run that went beyond the observed value. These minima also move with small changes: on the table before the recovered annotations, the smallest BP `strata` p-value was above the observed one, and the smallest CC `strata` p-value was 1.9e-9. Extracellular region's strength under `strata` shows that it can arise from gene properties alone, because `strata` keeps dN/dS tied to length and dS, so positive genes stay short. It should not be read as evidence about extracellular proteins. BH on CC under `strata` reaches 18% of runs with a significant term, again because of that link.

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

Unrelated outputs from the two candidates are byte-identical: `AtAl_dnds.tsv`, `AtAl_dnds_summary.tsv`, and both ideogram SVGs. On the corrected code, the annotated table and every `term_enrichment`, `ipr_enrichment`, and `go_enrichment` table are byte-identical at 1 and 8 threads.

## What selects the candidates

The positive genes are not a random sample ([data](data/enrichment_confounders.tsv)):

- **Length.** Positive rate falls from 3.05% in the shortest fifth of genes to 0.09% in the longest. Positive genes have a median of 125.5 amino acids, and others 364.
- **dS.** It falls from 3.10% in the lowest fifth to 0.16% in the highest: dN/dS is a ratio, and is unstable when dS is small.
- **Annotation.** 0.66% of annotated genes are positive, and 3.98% of unannotated genes. So 72 of 206 positive genes (35%) cannot contribute to any enrichment test.

Under `strata`, which holds these constant, the corrected Fisher tests stay near 5%. That does not rule out confounding in the real labels: a significant term can still reflect short genes, so report term results with their gene lengths in mind.

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
| `ipr_enrichment` looks up the current InterPro release even with `entries_source = "local"`, and stops if it cannot connect, with a message that local files would avoid it | No wrong results. Set `auto_detect_release = FALSE` | Open: [known issues](known-issues.md) |

## Limits

One genome pair, at one divergence. The correction controls false positives under no association; it does nothing for confounding when an association is real. `permute` and `rotate` do not hold confounders constant, and `strata` uses coarse strata and leaves dN/dS linked to length and dS within them. Nothing here measures power. No result has been compared with independent biological evidence. Modeling confounders and biological benchmarking stay on the [v1.0.0 roadmap](../dev/ROADMAP-1.0.md). The rank test in `enrichment-calibration.R` was run but is not interpreted here. `enrichment-verify.py` does not cover parent-child, whose expected values come from the package's constructed tests. Old-code runs were fewer than new-code runs, and parent-child old-code strata has 3 runs of 50 stopped by the error for no parent-child edges, left out and counted.

## Rebuild

Run in the image, with `dndsR` installed, after building the table with `annotate-gff.py` and the standard workflow. `enrichment-calibration.R` checks its reimplementation against a `term_enrichment` table from the old code, so give it one from `59faf30`.

```
Rscript tools/dev/validation/enrichment-real-null.R <dnds_annot.tsv> <runs per null> <seed> <outdir> [InterPro tree]
IPR_ADJUST_SCOPE=global Rscript tools/dev/validation/enrichment-real-null.R <dnds_annot.tsv> <runs per null> <seed> <outdir>
python3 tools/dev/validation/enrichment-null-summary.py <label>=<outdir> ...
Rscript tools/dev/validation/enrichment-calibration.R <dnds_annot.tsv> <old IPR enrichment.tsv> <term column> 5000 20261101 <outdir>
Rscript tools/dev/validation/enrichment-go-null.R <dnds_annot.tsv> <task> <runs per null> <outdir>
python3 tools/dev/validation/enrichment-verify.py term <dnds_annot.tsv> <table> [min_total] [min_pos]
python3 tools/dev/validation/enrichment-verify.py ipr <dnds_annot.tsv> <entry list> <dir> <prefix> [min_total] [min_pos] [max_prop]
python3 tools/dev/validation/gff-recovery.py <old dnds_annot.tsv> <new dnds_annot.tsv> <truth.tsv> <query GFF>
python3 tools/dev/validation/enrichment-before-after.py <old dir> <new dir> <observed_terms.tsv> <entry list>
```
