# Validation on a public genome pair

*Arabidopsis thaliana* against *A. lyrata*, using public RefSeq assemblies, comparing the code used for earlier pre-releases with the current code. Anyone can reproduce this. Tables: [data/](data/).

| Item | Detail |
|---|---|
| Query | *A. thaliana* TAIR10.1 RefSeq, GCF_000001735.4, with chromosomes renamed `Chr1` and so on. SHA-256 prefixes `6d1af248f5361fad` (FASTA), `3243386483ba94ca` (GFF) |
| Subject | *A. lyrata* v.1.0 RefSeq, GCF_000004255.2, NCBI checksums verified. `5f6b7b6bab223eeb` (FASTA), `d8e60f85ef791b5b` (GFF) |
| Old code | `c8b534f`, with the same CDS logic as v0.1.0 and v0.2.0 |
| Current code | `59faf30` |
| Environment | Published image, R 4.5.2, orthologr 0.4.2, DIAMOND 2.1.14 |
| Independent checks | gffread 0.12.9, PAML yn00 4.10.10, KaKs_Calculator 2.0.1 |

Scripts: `tools/dev/validation/run-version.sb`, `compare-versions.R`, `real-estimators.R`, `isoform-effects.R`.

## CDS against gffread

For the transcript dndsR chose. "Compared" leaves out 83 *A. lyrata* transcripts (mostly single-exon) that have no name in common with gffread's output, so they could not be matched.

| Genome | Gene class | Genes | Written | Current: match, differ | Old: match, differ |
|---|---|---:|---:|---:|---:|
| *A. thaliana* | Multi-exon, minus strand | 5,655 | 5,645 | 5,645, **0** | **0, 5,655** |
| | Multi-exon, plus strand | 5,659 | 5,644 | 5,644, **0** | 5,659, 0 |
| | Multi-isoform | 10,695 | 10,695 | 10,695, **0** | **0, 10,695** |
| | Single exon | 5,435 | 5,433 | 5,433, **0** | 5,435, 0 |
| *A. lyrata* | Multi-exon, minus strand | 9,135 | 8,848 | 8,839, **0** | **0, 9,125** |
| | Multi-exon, plus strand | 8,923 | 8,667 | 8,663, **0** | 8,918, 1 |
| | Multi-isoform | 5,868 | 5,868 | 5,868, **0** | **0, 5,868** |
| | Single exon | 5,891 | 5,820 | 5,750, **0** | 5,819, 2 |

Current code matches gffread on every gene compared. CDS left out were invalid: 580 *A. lyrata* CDS had in-frame stops, which is a property of that annotation.

## dN/dS, old against current

| | Old | Current |
|---|---:|---:|
| Ortholog pairs | 20,104 | 22,146 |
| Ratio undefined | 671 | 55 |
| dN/dS above 1 | **1,503** | **206** |
| dN/dS at least 10 | 8 | 0 |
| dS above 2 | 214 | 30 |
| Median dS | 0.149 | 0.141 |
| Median dN/dS | 0.319 | 0.187 |
| Infinite ratios | 0 | 0 |

Among 17,300 pairs with the same partner and a defined ratio in both versions:

- **10,740 (62%) changed CDS.** Their median absolute change in the ratio is 0.17, up to 20.3.
- **6,560 did not change CDS, and every one of their values is identical,** so the changes come from the sequence fix, not numerical variation.
- **4,906 changed more than twofold.**
- **1,206 lost a ratio above 1, and 36 gained one:** old 1,331, current 161.

Pairing: 17,859 pairs kept the same partner, 1,166 changed partner, 1,079 were only in the old results, and 3,121 only in the current ones.

**This means roughly nine in ten old calls above 1 in this comparison were artifacts of mis-assembled sequence.**

## Independent estimators on current results

800 pairs, with those above 1 oversampled, on one separately built alignment: global protein alignment, back-translated to codons.

| Method | Spearman with dndsR, ratio | dN | dS | Same side of 1 | Median ratio: method, dndsR |
|---|---:|---:|---:|---:|---:|
| Comeron, same alignment | 0.99 | 0.99 | 0.95 | 97% | 0.267, 0.254 |
| PAML YN00 | 0.98 | 0.99 | 0.90 | 96% | 0.263, 0.254 |
| PAML NG86 | 0.98 | 0.99 | 0.93 | 94% | 0.247, 0.255 |
| KaKs_Calculator MLWL | 0.99 | 0.99 | 0.92 | 91% | 0.264, 0.260 |
| KaKs_Calculator MYN | 0.97 | 0.99 | 0.89 | 95% | 0.254, 0.260 |

These tools did not share code with dndsR, and agree closely. On this enriched sample, 3 to 9% of pairs fall on opposite sides of 1 between methods, which shows why 1 is a weak cutoff for a single gene.

## Determinism

The dN/dS table is byte-identical between `-t 32` and `-t 8` (`ad656dea0495cf1e`). `extract_cds` output for *A. thaliana* is unchanged by the transcript-name fix (`5366bf48601924f6`).

## Defect found by this work

A public RefSeq genome where only some transcripts have `transcript_id` made the first version of the new `extract_cds` stop. It failed clearly and did not give wrong results. Fixed, with tests.

## Shared identifiers and annotation

No gene IDs are shared between these two genomes. RefSeq-style IDs matched, so annotations and coordinates are present in the current table. Enrichment on this pair, with public GO and InterPro terms added, is in [enrichment.md](enrichment.md).

## Isoform choice

How much does choosing the longest isoform matter? Of 22,146 pairs, 10,576 (48%) have at least one multi-isoform gene. For a random sample of 1,000 of them, every valid isoform combination (up to 16 per pair) was aligned separately with the Comeron method. 996 had at least two combinations with a defined ratio.

| Spread in dN/dS across isoform combinations | Share or value |
|---|---:|
| Median range | 0.018 |
| 90th percentile of range | 0.18 |
| Range above 0.1 | 22% |
| Maximum at least twice the minimum | 14% |
| Crosses 1 | 0.3% (3 pairs) |

Isoform choice rarely matters, but it matters noticeably for about a fifth of multi-isoform pairs, or roughly one in ten of all pairs. It seldom moves a gene across 1 here, partly because most genes here are far from 1. The longest isoform is typically close to the middle of the combinations: the median absolute difference from the median combination is 0.0025. This supports the longest-isoform default as reasonable, not as uniquely correct.

Limits: this uses one genome pair with one alignment method, and most annotated isoforms in these two annotations come from automated pipelines.
