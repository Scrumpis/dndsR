# Image stacks compared

**Do not change image partway through an analysis.** The same dndsR code on the v0.2.0 and v0.3.0 images chose a slightly different set of pairs, because tiny E-values are read differently. Values for shared pairs were identical.

| | v0.2.0 image, `sha256:726657e9` | v0.3.0 image, `sha256:feb5e4b1` |
|---|---|---|
| R | 4.5.2 | 4.6.2 |
| orthologr | 0.4.2 | 0.4.3, pinned to `06c3098` |
| Biostrings | 2.78.0 | 2.80.2 |
| topGO | 2.62.0 | 2.64.0 |
| DIAMOND | 2.1.14 | 2.1.14 |

## Measured

*A. thaliana* against *A. lyrata*, RefSeq, from the same inputs. The v0.2.0 image ran the v0.3.0 code, and the v0.3.0 image ran the v0.4.0 candidate. The code that calculates these did not change between them.

| Output | Result |
|---|---|
| CDS and CDS QC tables | Byte-identical |
| Pairs | 22,146 and 22,149, with 22,140 shared: 6 only on the old stack, 9 only on the new |
| dN, dS, and dN/dS of the 22,140 shared pairs | All identical |
| Pairs above dN/dS 1, shared pairs | 206 and 206, none crossing |
| `evalue` | Differs in 598 rows, all below about 1e-305 |
| `covariate_enrichment` InterPro | 61 written rows on both, 16 significant on both; largest change in a p-value 0.00013 |
| `covariate_enrichment` GO | 0 significant on both. 210 and 211 rows: GO.db and pairs differ, and non-significant p-values move, by up to 0.24 |
| `go_enrichment` and `ipr_enrichment` tables | Differ, as a result of the above and of GO.db |

## Cause

DIAMOND reports E-values below 1e-305. On the old stack, `7.79e-315` appears in results as `7.79e-305`, and `1.44e-318` as `1.44e-305`. `readr` does exactly this, in both images, while `data.table::fread` and base R read these correctly. Best hits are chosen by the lowest E-value, then the longest alignment, so among very similar long genes the old stack could rank hits wrongly.

Each of the 6 pairs lost has a replacement with a higher bit score, for example AT2G36530: 835 at 93.7% identity on the old stack, and 863 at 98.0% on the new. So the newer stack gives better pairs, in a small share of genes. That `readr` was the reader on the old stack, and not on the new, is an inference from these values: it was not traced through orthologr's dependencies.

**Consequences:** the 22,146 and related counts in earlier evidence come from the old stack. Calibration evidence is unaffected in its conclusion, as below. Ties remain possible at E-value 0, where the longest alignment then decides. Other install routes hold other versions, so their results may differ in the same way. Record the image digest, and check the versions in `<comparison>_dnds_provenance.json`.

Numbers for the old stack: `final/new/AtAl`, tables from code at `eed1096`. New stack: candidate `5adbc5b`.

## Thread counts

On the v0.3.0 image with the candidate, 8 and 32 `calculate_dnds` threads gave byte-identical tables and FASTA, 19 files in all. The only difference was the timestamp in `AtAl_interpro_provenance.json`.
