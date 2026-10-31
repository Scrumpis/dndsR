# CDS extraction: exon order and isoform handling

## Problem

`extract_cds()` with its default `group_by = "gene"` called `GenomicFeatures::cdsBy(txdb, by = "gene")`. That returns the CDS parts of every isoform of a gene together, in ascending genomic order. Two consequences:

- **Multi-exon genes on the minus strand were assembled in the wrong exon order**, usually producing out-of-frame sequence.
- **Multi-isoform genes were concatenated into chimeras** containing all isoforms.

This code dates from September 2025 and is in v0.1.0, v0.2.0, `main`, and `dev` up to `c8b534f`. Results from those versions are affected, in any genome with such genes. Sequence pairs stay comparable because both species are garbled the same way, so the problem is not obvious from the output.

## Fix

Coding sequences are now grouped by transcript, in transcript order. With `group_by = "gene"`, each gene is represented by its longest isoform with a valid CDS. A valid CDS has a length divisible by 3 and no in-frame stop. Invalid CDS are excluded and recorded in `<genome>_CDS_qc.tsv`.

## Evidence

Designed example: `inst/extdata/example/` (generator `tools/dev/validation/make-example.R`). Its true CDS are written before the genomes are laid out. Reference values come from `tools/dev/validation/ng86.py`, an independent Nei-Gojobori calculation.

Using `calculate_dnds` defaults (Comeron method), original code `c8b534f` versus the fix:

| Gene | Structure | dN/dS before | dN/dS after | NG86 on true CDS | Query protein length before and after |
|---|---|---:|---:|---:|---|
| g28 | 3 exons, minus strand, split mid-codon | **1.615** | 0.1024 | 0.1123 | 302, 302 |
| g29 | Exon-skipping isoform | 0.09304 | 0.09985 | 0.1119 | 504, 302 |
| g30 | Extended exon isoform | 0.09414 | 0.09453 | 0.0971 | 619, 317 |
| g35 | Shortened exon isoform | 0.09218 | 0.08993 | 0.0987 | 589, 302 |
| g36 | Longest isoform has in-frame stop | 0.09784 | 0.09823 | 0.1027 | 619, 302 |
| g31 | In-frame stop in subject | 0.1092 | excluded | not computable | 302, none |

The other 27 pairs have identical values before and after. In this one example, g28's true dN/dS of about 0.11 became an apparent signal of positive selection. The size of the effect in real genomes has not been measured.

After the fix, every written CDS matches its designed sequence, except g30, which by design uses its longer valid isoform (`tests/testthat/test-example-pipeline.R`). Results are byte-identical for `-t 1` and `-t 4`.

Environment: image digest in `skills/dndsr/SKILL.md`, GenomicFeatures 1.62.0, Biostrings 2.78.0, orthologr 0.4.2, DIAMOND 2.1.14.

## Action for users

Rerun `extract_cds` and `calculate_dnds` for any analysis made with an earlier version. `extract_cds` regenerates old CDS automatically. `calculate_dnds` keeps existing results and says so unless `--overwrite` is given.

## Remaining limits

A CDS whose first part has a nonzero GFF phase (a partial gene) is not trimmed to its reading frame. If its length is a multiple of 3 it may be read out of frame, and an in-frame stop would then usually exclude it. Selecting the longest isoform is a convention, and it may not be the isoform under selection.
