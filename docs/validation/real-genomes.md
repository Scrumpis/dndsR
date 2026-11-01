# Real-genome validation: old code versus current code

## Summary

On real genomes, the coding-sequence bug in dndsR versions up to `c8b534f` (including v0.1.0 and v0.2.0) is not a rare edge case:

- **Every multi-exon gene on the minus strand was assembled wrongly, and so was every gene with more than one isoform.** That was 40% of genes in *C. formosanum* and *C. quinoa* (about 12,000 per genome), and 60% in *A. thaliana*.
- **About half of the comparable gene pairs changed their CDS** (49% and 47%).
- **Most genes reported with dN/dS above 1 were artifacts.** Over the whole table, calls fell from 1,241 to 436 in *C. formosanum* and from 1,049 to 212 in *C. quinoa*.
- **About 15% more ortholog pairs were found** after the fix, because garbled sequences had lost their matches.
- **Significant enrichment terms changed almost entirely** (Jaccard similarity of significant sets: 0 to 0.09).
- The current code's CDS match independent `gffread` output for **all** tested genes, and its dN/dS agree closely with PAML and KaKs_Calculator (Spearman 0.95 to 0.98 for the ratio).

Earlier analyses are affected. An existing *C. formosanum* result table in the maintainer's project directory has the same dN/dS columns as a fresh run of the old code, so this is the behavior those analyses had.

## Method

| Item | Detail |
|---|---|
| Old code | `c8b534f` |
| Current code | `f6cde64` for CDS and dN/dS, whose code is unchanged since. `fedd59e` for annotation and enrichment |
| Image, both versions | digest in `skills/dndsr/SKILL.md`; R 4.5.2, orthologr 0.4.2, DIAMOND 2.1.14 |
| Comparison 1 | *C. formosanum* subgenomes B and C (`CheFo_v3`) |
| Comparison 2 | *C. quinoa* subgenomes B and A (`CheQu_v2`, a public genome) |
| CDS only | *A. thaliana* TAIR10.1 RefSeq, which has isoforms and RefSeq-style IDs |
| Independent CDS | gffread 0.12.9 |
| Independent estimators | PAML yn00 4.10.10, KaKs_Calculator 2.0.1 |
| Input hashes (SHA-256, first 16) | `CheFo_v3_B.gff` `f4de9c4277d9814c`, `CheFo_v3_C.gff` `69ebf95a4999344e`, `CheQu_v2_A.gff` `447b4bc7c598e486`, `CheQu_v2_B.gff` `17c02c1b03402371`, formosanum B and C FASTA `8c42700b3618b582` and `3d2e2b0609750f1a`, quinoa A and B FASTA `abb700f8fdd63abf` and `deea035314ea91ae` |

Both versions used identical inputs and options through Slurm, with data mounted read-only. Scripts: `tools/dev/validation/run-version.sb`, `compare-versions.R`, `real-estimators.R`. Aggregate tables: [real-genomes/](real-genomes/). Gene-level results are kept outside Git.

## CDS extraction

Against gffread's sequence for the transcript dndsR chose:

| Comparison | Gene class | Genes | Current matches | Old matches |
|---|---|---:|---:|---:|
| *C. formosanum* (both genomes) | Multi-exon, minus strand | 23,989 | 23,989 | **0** |
| | Multi-exon, plus strand | 23,939 | 23,939 | 23,939 |
| | Single exon | 11,403 | 11,403 | 11,403 |
| *C. quinoa* (both genomes) | Multi-exon, minus strand | 23,397 | 23,397 | **0** |
| | Multi-exon, plus strand | 23,496 | 23,496 | 23,496 |
| | Single exon | 11,779 | 11,779 | 11,779 |
| *A. thaliana* | Multi-exon, minus strand | 5,655 | all written | **0** |
| | Multi-isoform | 10,695 | all written | **0** |
| | Plus strand, or single exon | 11,094 | all written | 11,094 |

In *A. thaliana*, 27,417 of 27,444 genes were written, and all 27,417 match gffread. The 27 not written are 22 with a CDS length not divisible by 3 and 5 with an in-frame stop. Old multi-isoform sequences had a median of twice the chosen isoform's length, and 10,678 of 10,695 contained internal stops, as did 5,629 of the 5,655 minus-strand genes. The Chenopodium annotations have one transcript per gene, so the isoform bug could not affect them.

## dN/dS

Pairs with the same partner and a defined ratio in both versions:

| | *C. formosanum* | *C. quinoa* |
|---|---:|---:|
| Pairs compared | 12,426 | 12,284 |
| CDS changed in either gene | 6,041 | 5,812 |
| Median absolute change where CDS changed | 0.19 | 0.20 |
| Largest absolute change | 21.4 | 20.7 |
| Pairs whose ratio changed over twofold | 2,705 | 2,733 |
| Ratio above 1: old, then current | 1,081, then 306 | 941, then 148 |
| Lost the above-1 call | 850 | 824 |
| Gained the above-1 call | 75 | 31 |
| **Pairs with unchanged CDS: value changed** | **0 of 6,385** | **0 of 6,472** |

The last row is important. Where the sequences did not change, every value is identical, so these changes come from the sequence fix and not from numerical variation.

Across whole tables:

| | Old | Current | Old | Current |
|---|---:|---:|---:|---:|
| | *C. formosanum* | | *C. quinoa* | |
| Pairs | 16,053 | 18,458 | 15,081 | 16,880 |
| Ratio undefined (NA) | 1,489 | 363 | 1,092 | 262 |
| dS above 2 | 347 | 52 | 310 | 92 |
| Ratio 10 or more | 9 | 0 | 4 | 0 |
| Ratio above 1 | 1,241 | 436 | 1,049 | 212 |
| Median ratio | 0.358 | 0.252 | 0.337 | 0.239 |

No infinite values occurred in either version. Mis-assembled genes appear as high dS and undefined ratios as well as false positives: frame loss makes sequences look saturated. Note that *C. formosanum* and *C. quinoa* subgenome comparisons are homeologs within a species.

## Pairing

| | *C. formosanum* | *C. quinoa* |
|---|---:|---:|
| Same partner in both | 13,709 | 13,129 |
| Different partner | 917 | 706 |
| Only in old | 1,427 | 1,246 |
| Only in current | 3,832 | 3,045 |

## Independent estimators on real pairs

800 pairs per comparison, with genes above 1 oversampled, on one independently built alignment: global protein alignment back-translated to codons, not from orthologr. Compared with dndsR's reported values:

| Method | Spearman, ratio | Agree on side of 1 | Median ratio against dndsR |
|---|---:|---:|---|
| Comeron, same alignment | 0.98, 0.98 | 96%, 98% | 0.36 against 0.34, 0.36 against 0.33 |
| PAML YN00 | 0.96, 0.95 | 91%, 92% | 0.34 against 0.34, 0.35 against 0.33 |
| PAML NG86 | 0.97, 0.96 | 89%, 91% | 0.33 against 0.34, 0.34 against 0.33 |
| KaKs_Calculator MLWL | 0.97, 0.96 | 90%, 92% | 0.37 against 0.34, 0.36 against 0.33 |
| KaKs_Calculator MYN | 0.95, 0.95 | 90%, 92% | 0.34 against 0.34, 0.34 against 0.33 |

Values are *C. formosanum*, then *C. quinoa*. Roughly 8 to 11% of this enriched sample falls on opposite sides of 1 between methods, which again shows why a threshold of 1 is a weak criterion for individual genes.

## Enrichment

Each version ran its own enrichment on its own results, with defaults, so these differences include the statistical fixes as well as the sequence fix. Significant is `p_adj` of 0.05 or less.

| Comparison, side | Set | Significant: old | Significant: current | Shared |
|---|---|---:|---:|---:|
| *C. formosanum*, query | InterPro | 4 | 7 | 0 |
| *C. formosanum*, subject | InterPro | 4 | 8 | 1 |
| *C. formosanum*, query | GO | 20 | 13 | 0 |
| *C. formosanum*, subject | GO | 22 | 16 | 0 |
| *C. quinoa*, query | InterPro | 16 | 6 | 0 |
| *C. quinoa*, subject | InterPro | 16 | 4 | 0 |
| *C. quinoa*, query | GO | 14 | 8 | 0 |
| *C. quinoa*, subject | GO | 21 | 9 | 0 |

These tables use current code at `fedd59e`, after the fix below.

Old enrichment conclusions should not be carried forward. The default `ipr_enrichment` chooses an InterPro release using a network lookup and records it in `*_interpro_provenance.json`.

## Defect found by this work

**RefSeq-style annotations got no annotations at all.** Bioconductor removes `gene-` and `rna-` prefixes, so `extract_cds` writes `AT1G01010`, but `append_annotations` looked for `gene-AT1G01010` and left every annotation empty, without any error. `dnds_ideogram` used the same lookup. Fixed, with regression tests, and `append_annotations` now stops if no IDs match and reports partial mismatches.

**A regression from this release, fixed.** A stricter term pattern lost every GO term after the first in lists encoded as `GO:1%2CGO:2`, the GFF3 standard. It was found in this comparison: 5,391 of 15,992 rows in one comparison lost terms. GFF3 escapes are now decoded first, with a regression test. After the fix, GO terms match the old code wherever the annotation leak does not apply: see [annotation-leak.md](real-genomes/annotation-leak.md).

## Duplicate identifiers across genomes

The tested comparison pairs share no IDs (0 shared, in query and subject of each). But gene IDs across two of the maintainer's assemblies, `CheAl_v01.0` and `CheAl_ToL_v01.0`, overlap by 80,767 of about 85,000, and those were run in one batch. See [real-genomes/annotation-leak.md](real-genomes/annotation-leak.md).

## What this does not show

- Whether the old bugs' effects are similar in other genomes. Strand and exon balance suggest they will be, but this was not measured elsewhere.
- That current dN/dS values are right in an absolute sense. Agreement between estimators on a shared alignment supports the pipeline, not the biology.
- Gene-model quality. Real annotations contain errors that no check here can see.
