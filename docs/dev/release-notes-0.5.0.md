# dndsR 0.5.0 release notes

Draft until tagged. Fill in the image digest once the publish workflow has run.

**Pairwise defaults are unchanged.** With default settings, every CDS, dN/dS table, and annotated table from whole genomes (the public *A. thaliana*-*A. lyrata* pair and two other pairs) is byte-identical to v0.4.0. What changes is listed under Fixed.

## New

- **Other dN/dS estimators**: `--dnds-method NG`, `LWL`, `LPB`, `MLWL`, `YN`, or `MYN`, through KaKs_Calculator 2.0.1 on the same alignments; `Comeron` stays the default. Identical to KaKs_Calculator run directly, with known differences from PAML recorded. [Evidence](../validation/estimators.md).
- **Reliability flags**: `<comparison>_dnds_reliability.tsv` beside each dN/dS table flags pairs with few codons, few synonymous differences, extreme dS, or dN/dS above 10. It changes no value. [Criteria](../output-files.md).
- **Selection tests across species, experimental**: `phylo_prepare` (genomes or CDS and orthogroups to codon alignments and labelled trees) and `phylo_busted` (HyPhy BUSTED on chosen branches, corrected across orthogroups). The wrapper reproduces HyPhy run directly. A simulation study on one five-taxon tree measured false positives and power, and recommends `--multiple-hits Double` and three starts (`--extra-seeds 2,3`), reporting the new `p_best_fits_adj`. That calibration does not cover other trees or real data, and alignment errors remain the main source of false positives. [Methods and limits](../phylo-methods.md).
- **McDonald-Kreitman test, experimental**, contributed and maintained by a collaborator: add `snp_file` and `snp_side` to a comparison row. Tested on designed examples only. [Method](../mcdonald-kreitman.md).
- **`--unique-pairs`** for multi-species comparison files: every step skips a repeated genome pair, unless the row carries its own McDonald-Kreitman inputs.
- The published image includes HyPhy 2.5.103 and KaKs_Calculator 2.0.1.

## Fixed

- **Gene and transcript names now always come from the GFF3 `ID`.** A GFF3 whose gene symbols were all unique was named by symbol, and `append_annotations` then left those genes unannotated (201 of 502 pairs in an *A. thaliana* subset); Ensembl-style IDs (`gene:...`) were never found. Whole RefSeq genomes keep their names; one QC column changes for transcripts without a `transcript_id`.
- **Header rows in comparison files work** in `calculate_dnds`. Before, a header was read as a comparison.
- **Sequence names DIAMOND cannot use** (whitespace, or an NCBI tag such as `gnl|` or `ref|`) now stop with a clear message. Before, runs failed with misleading errors or, for `gi|` names, returned an empty table.
- `calculate_dnds` refuses `--dnds-method` values it does not support, instead of passing them on unchecked.

## Testing

A small public real-genome fixture now runs the whole pipeline in about 2 minutes in CI on every change to the code, byte for byte against a baseline; whole genomes are rerun for releases and core changes. [Testing strategy](TESTING.md).

## Install

See [INSTALL.md](../../INSTALL.md). Image digest: fill in at release.
