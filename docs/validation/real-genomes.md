# Real-genome validation

## Result on a public genome

Old code (`c8b534f`, the same CDS logic as v0.1.0 and v0.2.0) against current code on *Arabidopsis thaliana* TAIR10.1 RefSeq, with gffread 0.12.9 as an independent source of CDS:

| Gene class | Genes | Current matches gffread | Old matches gffread |
|---|---:|---:|---:|
| Multi-exon, minus strand | 5,655 | all written | **0** |
| More than one isoform | 10,695 | all written | **0** |
| Plus strand with one isoform, or single exon | 11,094 | all written | 11,094 |

27,417 of 27,444 genes were written, and all 27,417 match gffread. The 27 not written are 22 with a CDS length not divisible by 3, and 5 with an in-frame stop. In the old code, 10,678 of the 10,695 multi-isoform sequences and 5,629 of the 5,655 minus-strand sequences contained internal stops, and multi-isoform sequences had a median of twice the chosen isoform's length.

## Further checks on unpublished data

The maintainer also compared old and current code end to end on unpublished genome pairs: pairing, dN, dS, dN/dS, annotation, and enrichment. Those data and results are not public, and will be reported with that work. In general terms, they support the following, which should be checked for any re-analysis:

- The structural errors above occurred at the same rate, for every affected gene.
- Pairs whose sequences did not change kept identical values, so changes come from the sequence fix, not numerical variation.
- Mis-assembled genes can appear as false dN/dS above 1, high dS, and undefined ratios, and can lose their orthologs.
- Downstream enrichment can change substantially, so old enrichment conclusions should not be carried forward.
- dndsR's values agreed closely with PAML and KaKs_Calculator on a separately built alignment.

## Method, for reproducing on your own data

Run both versions on identical inputs with data mounted read-only, using `tools/dev/validation/run-version.sb`, then `compare-versions.R` (gene classes, gffread comparison, pairing, changes in dN and dS, problem values, shared IDs) and `real-estimators.R` (PAML yn00, KaKs_Calculator). Keep results and unpublished data outside Git.

## Defects found this way

- **RefSeq-style annotations got no annotations.** Bioconductor removes `gene-` and `rna-` prefixes, so `extract_cds` writes `AT1G01010`, but `append_annotations` looked for `gene-AT1G01010` and left every annotation empty without an error. `dnds_ideogram` used the same lookup. Fixed with regression tests, and `append_annotations` now stops if no IDs match.
- **A regression in this release, fixed before release.** A stricter term pattern lost every GO term after the first in lists encoded as `GO:1%2CGO:2`, the GFF3 standard. GFF3 escapes are now decoded first, with a regression test.
- **Annotation leak.** At one thread, a gene ID shared by two assemblies in a batch took coordinates and terms from the first one read. Results depended on thread count. This was confirmed on real data, and current code gives the same result at 1 and 8 threads.

## What this does not show

- That current dN/dS values are right in an absolute sense. Agreement between estimators supports the pipeline, not the biology.
- Gene-model quality. Real annotations contain errors that no check here can see.
- Results for genomes other than those tested, though the bugs are structural rather than species-specific.
