# Enrichment calibration

Do the enrichment tests give false positives at the rate they claim? This page tests them on a real, publicly annotated genome pair. Everything here can be rebuilt with the scripts named below.

**Result: at their defaults, `term_enrichment` and `ipr_enrichment` are not calibrated.** With dN/dS assigned to genes at random, where any significant term is false, they report one in 58% and 95% of runs. The cause is a filter applied before multiple-testing correction. Correcting first brings the error rate close to its nominal level. `go_enrichment` is not affected by this, but its default list at p < 0.05 cannot be told apart from noise by its size.

## Data

| Item | Source | SHA-256 (first 16) |
|---|---|---|
| Genomes | *A. thaliana* (query, RefSeq TAIR10.1) and *A. lyrata* (subject), as in [public-genome.md](public-genome.md) | See that page |
| GO terms | GO Consortium TAIR GAF, without `NOT` or `ND` annotations | `582ea7a8cbbae942` |
| InterPro terms | UniProt reference proteome UP000006548 | `159a133c5b9fb88f` |
| Per-gene table analyzed | `AtAl_dnds_annot.tsv`, from `annotate-gff.py`, then the standard workflow | `603fdcd5d8f67ab7` |

Code: commit `59faf30` (`R/` has not changed since), image `dndsr_v0.2.0_sha256-726657e9`, seed 20261101.

The workflow ran end to end. There are 22,146 pairs, of which 22,091 pass `max_dnds`, and 206 have dN/dS > 1. Terms appended by `append_annotations` were checked against the source files: GO terms are identical for 22,075 pairs. Of the other 71, all lost terms to the [defect below](#defects-found); none gained any.

## Null models

The dN/dS column is reassigned, so no term can be truly associated with it:

| Null | What it keeps | Why |
|---|---|---|
| `permute` | Nothing | The usual label permutation. It ignores all structure |
| `strata` | Gene length, dS, and annotation richness (permuted within strata) | Holds the confounders found below constant |
| `rotate` | Genomic order (circular shift) | Keeps tandem clusters and local dependence together |

## Term and InterPro enrichment

The real functions, on `permute`, 200 runs each ([data](data/enrichment_real_null.tsv)):

| Function | Runs with any significant term | Mean significant terms | Mean terms tested |
|---|---|---|---|
| `term_enrichment` | 57.5% | 4.25 | 102 |
| `ipr_enrichment` | 94.5% | 6.16 | 76 |

With no association, about 5% of runs should have any. A reimplementation, verified identical to `term_enrichment` on the 67 observed terms (largest p-value difference 4.4e-16), allowed 5,000 draws per null ([data](data/enrichment_null_calibration.tsv)):

| Null | Default: filter, then BH | BH on all terms with ≥ 2 genes | ≥ 5 genes | ≥ 10 genes |
|---|---|---|---|---|
| `permute` | 61.0% | 0.7% | 1.3% | 1.2% |
| `strata` | 90.5% | 3.6% | 5.5% | 5.5% |
| `rotate` | 66.0% | 1.0% | 1.8% | 1.6% |

**Cause.** `min_pos` removes terms with few positive genes before BH. That count is the test's outcome, so only terms that could be significant remain, and BH corrects for about 100 tests instead of about 7,200. Among all terms, 2% have p < 0.05, because Fisher's test is conservative on small counts. Among those kept, 26–34% do. This was documented as "can make adjusted p-values less conservative"; measurement shows that understated it.

**Smallest correction.** Apply BH across all terms meeting `min_total`, which does not depend on the outcome, and use `min_pos` only to choose which rows are reported. This keeps the interface. It changes reported `p_adj`, so it needs the maintainer's approval.

Effect on the real data:

| Procedure | Significant InterPro terms |
|---|---|
| Default | 48 of 67 tested |
| BH on all terms, ≥ 2 / ≥ 5 / ≥ 10 genes | 18 / 19 / 24 |
| Westfall-Young, by `permute` / `strata` / `rotate` | 16 / 9 / 16 |

The strongest terms survive every procedure, such as IPR011043 (galactose oxidase/kelch β-propeller, 13 of 160 genes, BH on all terms 3.5e-07) and IPR001810 (F-box, 16 of 484). Collapsing tandem clusters kept 45 of the 48 under the default filter, so tandem duplication is not what inflates this list. Whether these terms make biological sense has not yet been compared with the literature.

## GO enrichment

The real `go_enrichment` at its defaults, 200 runs per null. It reports unadjusted `weight01` p-values, and does not filter by `min_pos`.

| Ontology (terms) | Threshold | Observed | `permute` mean | `rotate` mean | `strata` mean |
|---|---|---|---|---|---|
| BP (1,956) | p < 0.05 | 19 | 12.3 | 11.9 | 5.2 |
| | p < 0.001 | 6 | 0.31 | 0.26 | 0.20 |
| MF (683) | p < 0.05 | 9 | 5.7 | 5.7 | 4.2 |
| | p < 0.001 | 4 | 0.17 | 0.15 | 0.20 |
| CC (331) | p < 0.05 | 3 | 3.0 | 3.0 | 2.5 |
| | p < 0.001 | 1 | 0.07 | 0.04 | 0.34 |

Observed counts at p < 0.05 are mostly inside the null's range (BP 5–95%: 6–21 under `permute`), so most of that list could be noise. At p < 0.001, BP and MF exceed the nulls by more than ten-fold; CC does not exceed `strata` by much. `weight01` p-values are conservative (well under 5% of terms fall below 0.05), so they are not ordinary p-values either way. The strongest p-value is not always beyond chance. In BP, cell-cell signaling (GO:0007267, 7 of 77, p = 1e-7) is below the smallest p-value of all 600 null runs (5.9e-6). In MF, hormone activity (p = 7.4e-6) is not: `permute` reached 5.8e-6. In CC, extracellular region (44 of 2,100, p = 3e-10) is below every run, but `strata` reached 1.9e-9, far beyond `permute` (3.6e-5). `strata` keeps dN/dS tied to length and dS, so positive genes stay short. A null that reaches this term's strength shows it can arise from those properties alone, with no link to function. It should not be read as evidence about extracellular proteins. BH on CC under `strata` also reaches 13% of runs with a significant term.

## What selects the candidates

The positive genes are not a random sample ([data](data/enrichment_confounders.tsv)):

- **Length.** Positive rate falls from 3.05% in the shortest fifth of genes to 0.09% in the longest. Positive genes have a median of 125.5 amino acids, and others 364.
- **dS.** It falls from 3.10% in the lowest fifth to 0.16% in the highest: dN/dS is a ratio, and is unstable when dS is small.
- **Annotation.** 0.66% of annotated genes are positive, and 3.98% of unannotated genes. So 72 of 206 positive genes (35%) cannot contribute to any enrichment test.

Under `strata`, which holds these constant, the corrected Fisher tests stay near 5%. That does not rule out confounding in the real labels: a significant term can still reflect short genes, so report term results with their gene lengths in mind.

The result also depends on the candidate rule ([data](data/enrichment_candidate_rules.tsv)):

| Rule | Significant terms (default procedure) | Shared with dN/dS > 1 | Jaccard |
|---|---|---|---|
| dN/dS > 1 (default) | 48 | 48 | 1 |
| dN/dS > 0.5 | 139 | 37 | 0.25 |
| dN/dS > 0.75 | 86 | 42 | 0.46 |
| dN/dS > 1.5 | 3 | 3 | 0.06 |
| dN/dS > 1 and ≥ 150 amino acids | 40 | 35 | 0.66 |
| dN/dS > 1 and dS ≥ 0.05 | 45 | 43 | 0.86 |

Length matters more than a dS floor here. These figures use the default procedure, so they are not corrected for the defect above.

## Defects found

| Defect | Effect | Status |
|---|---|---|
| `min_pos` filters before BH in `term_enrichment` and `ipr_enrichment` | Severe, above | Needs a decision |
| The GFF reader treats `#` inside attributes as a comment. Legal GFF3 | 79 lines in the *A. thaliana* GFF are truncated, and 71 pairs silently lose their terms | Open: treat only a `#` at line start as a comment |
| `go_enrichment(include_definition = FALSE)` always stops with an error | Explicit failure, no wrong results | Open |

## Limits

One genome pair, at one divergence. `permute` and `rotate` do not hold confounders constant, and `strata` uses coarse strata. Parent-child mode in `ipr_enrichment` was not assessed separately. The rank test in `enrichment-calibration.R` was run but is not interpreted here. Calibration under the global null says nothing about power.

## Rebuild

Run in the image, with `dndsR` installed, after building the table with `annotate-gff.py` and the standard workflow:

```
Rscript tools/dev/validation/enrichment-calibration.R <dnds_annot.tsv> <real IPR enrichment.tsv> <term column> 5000 20261101 <outdir>
Rscript tools/dev/validation/enrichment-go-null.R <dnds_annot.tsv> <task> <runs per null> <outdir>
Rscript tools/dev/validation/enrichment-real-null.R <dnds_annot.tsv> <runs> <seed> <outdir>
```
