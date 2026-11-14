# McDonald-Kreitman test

The McDonald-Kreitman (MK) test compares polymorphism within a species with divergence between species. Under neutrality, the ratio of nonsynonymous to synonymous changes should be the same in both, so a 2×2 table of polymorphic and fixed differences, synonymous and nonsynonymous, can be tested with Fisher's exact test.

dndsR runs the simplest form: counts only, no site-frequency weighting, on 1:1 orthologs.

## Turning it on

Add two optional columns to the comparison file:

| Column | Meaning |
|---|---|
| `snp_file` | Path to a `.vcf`, `.vcf.gz`, `.bed`, or `.bed.gz` of SNP positions |
| `snp_side` | `query` or `subject`: which species' genome the SNP coordinates refer to |

```
comparison_name query_fasta query_gff subject_fasta subject_gff snp_file snp_side
spA_vs_spB      spA.fasta  spA.gff3  spB.fasta      spB.gff3      spA.vcf  query
```

Both are needed together; `snp_side` must be `query` or `subject`. Without them nothing changes. Headerless comparison files may add them as columns 6 and 7.

Requirements:

- The CDS FASTAs written by `extract_cds` (normally already produced before `calculate_dnds`).
- The GFF3 and genome FASTA of the SNP side, to map SNP coordinates onto codons. They must already be in the comparison file.
- MK runs while dN/dS is computed, because it reuses the codon alignments in the temporary directory. A cached dN/dS result (inputs unchanged) skips MK; rerun with `--overwrite` to recompute both.

## SNPs used

Only bi-allelic SNPs are used, in both the polymorphism counts and the divergence exclusion. In a VCF that means REF is a single A/C/G/T base and there is exactly one ALT allele, also a single A/C/G/T base, different from REF. Multi-allelic records, indels, symbolic alleles, and REF = ALT records are dropped. In a BED file, only 1-bp intervals are kept, and coordinates are read as 0-based half-open and converted to 1-based.

Genotype columns are ignored: only the count of SNP positions matters.

A polymorphic site counts once however many samples carry it. Sites in a VCF are classified from their alleles, on the CDS strand; BED sites have no alleles, so they are classified by codon degeneracy: a position whose three alternative bases are all synonymous is synonymous (4-fold), one where no alternative is synonymous is nonsynonymous (0-fold), and 2- and 3-fold positions are skipped. Sites whose reference codon is a stop are excluded.

## Divergence counted

Dn and Ds come from the same PAL2NAL codon alignments `orthologr::dNdS()` produced, one aligned codon at a time:

- Codons with a gap, a non-ACGT base, or a stop codon on either side are skipped.
- A codon that carries a SNP in the SNP-side species is excluded from the divergence counts, so a site is not counted as both polymorphic and divergent.
- A codon pair that differs at one position counts 1 synonymous or 1 nonsynonymous difference. A pair differing at two or three positions is averaged over every shortest mutational path, with equal weight: each path steps through the intermediate codons, and a path that passes through a stop codon is discarded.

Pairs whose query or subject gene appears in more than one pair are dropped, so only 1:1 orthologs are counted.

## Statistics

With `Pn`, `Ps` (polymorphic synonymous and nonsynonymous sites) and `Dn`, `Ds` (divergent):

| Statistic | Formula |
|---|---|
| `alpha` | `1 - (Ds * Pn) / (Dn * Ps)` |
| `neutrality_index` | `(Pn * Ds) / (Ps * Dn)` |
| `dos` (direction of selection) | `Pn / (Pn + Ps) - Dn / (Dn + Ds)` |

Fisher's exact test uses the 2×2 table `[Pn, Ps; round(Dn), round(Ds)]`: exact tests need integer counts, so the path-averaged Dn and Ds are rounded (`R`'s round-half-to-even). `round()` values are written as `dn_int` and `ds_int`. The p-value is NA when a margin is zero.

`fisher_q` is Benjamini-Hochberg across the finite p-values in the file. The pooled row is a separate single test on the summed counts, not a combination of per-gene p-values.

## Checks on the SNP file

Because a SNP file with the wrong coordinates would silently produce nonsense, the run checks it:

1. Contigs. No contig in the SNP file shared with the SNP-side genome is an error, with a hint if the contigs match the other side instead. Contigs absent from the genome are dropped with a message.
2. Contig lengths, when the VCF header has `##contig` lines with lengths: a disagreement with the FASTA index stops the run. BED input has no such header, so this check is skipped.
3. REF alleles: each SNP that maps into a CDS is compared with the genome base (complemented on the minus strand). If at least 10 sites were checked and more than 10% disagree, the run stops; above 1% it warns.

These checks cannot detect a wrong assembly that shares contig names and lengths. Confirm the SNP file is from the genome you think it is.

## Output

Two files beside the dN/dS table:

- `<comparison>_mk.tsv`, one row per pair
- `<comparison>_mk_pooled.tsv`, one row for the whole comparison

Columns: [output-files.md](output-files.md#mcdonald-kreitman).

## Limits

- One population sample, unweighted: no correction for sample size, allele frequencies, or slightly deleterious polymorphisms (no DoS-style filtering).
- Only biallelic SNPs are used; multiallelic and indel records are ignored.
- Codon-level only: SNPs outside CDS regions are counted nowhere.
- 1:1 orthologs only. Genes in many-to-many pairs are dropped.
- With few substitutions the test has little power; a non-significant result is not evidence of neutrality.
