# McDonald-Kreitman test

The McDonald-Kreitman (MK) test compares polymorphism within a species with divergence between species. Under neutrality, the ratio of nonsynonymous to synonymous changes should be the same in both, so a 2×2 table of polymorphic and fixed differences, synonymous and nonsynonymous, can be tested with Fisher's exact test.

dndsR runs the simplest form: counts only, no site-frequency weighting, on 1:1 orthologs. Sites are classified by codon degeneracy: 4-fold positions count as synonymous, 0-fold as nonsynonymous, and 2- and 3-fold positions are left out of both classes.

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

## Filter your SNP calls first

dndsR does no SNP filtering. It does not look at `FILTER`, quality, depth, allele count, or genotypes, so whatever records you give it are treated as real polymorphisms. Filter before you run; `bcftools` already does this better than a package would:

```
# 1. Only if AC/AN are missing, add them first.
bcftools +fill-tags calls.vcf.gz -- -t AC,AN -Oz -o calls_tagged.vcf.gz

# 2. PASS calls, SNPs only, bi-allelic, ALT seen at least twice (no singletons), not fixed in the sample.
bcftools view -f PASS -v snps -m2 -M2 -i 'AC>=2 && AC<AN' calls_tagged.vcf.gz -Oz -o snps.vcf.gz
```

`AC<AN` drops sites where every sampled allele is ALT: those are fixed in the sample, not polymorphisms. `AC>=2` also drops singletons (ALT seen once). That is a frequency filter, and a choice: it removes some sequencing errors and some rare, often slightly deleterious, variants, which changes `Pn` and `Ps`. To keep singletons, use `-i 'AC>=1 && AC<AN'`. Report which you used.

With a single diploid sample, `AN` is 2, so `AC>=2 && AC<AN` keeps nothing. Use `-i 'AC>=1 && AC<AN'` (or, equivalently there, `-i 'GT="het"'`), which keeps its heterozygous sites. One individual gives few polymorphisms, and the MK test has little power with them.

One low-quality call can move a gene from `Pn` to `Ps`. The same applies to BED input: filter first, then convert.

## SNPs used

Only bi-allelic SNPs are used, in both the polymorphism counts and the divergence exclusion. In a VCF that means REF is a single A/C/G/T base and there is exactly one ALT allele, also a single A/C/G/T base, different from REF. Multi-allelic records, indels, symbolic alleles, and REF = ALT records are dropped. In a BED file, only 1-bp intervals are kept, and coordinates are read as 0-based half-open and converted to 1-based.

Genotype columns are ignored: only the count of SNP positions matters.

A polymorphic site counts once however many samples carry it. Every site, VCF or BED, is classified the same way, by the degeneracy of the changed position in the reference codon: if all three alternative bases are synonymous the position is 4-fold and the site is a synonymous polymorphism (`Ps`); if none is, it is 0-fold and the site is nonsynonymous (`Pn`). 2- and 3-fold positions depend on the allele, so they are excluded from both classes and counted in `snps_midfold`. Sites whose reference codon is a stop are excluded.

**Known limitation.** The class is read from the codon the count starts from, which for divergence is the `snp_side` species' codon. For 38 of 526 single-difference codon pairs (7.2%) the other codon gives a different answer — always 0-fold against mid-fold, never 4-fold; `AGA`→`GGA` at position 1 is one. On 1,472 public *A. thaliana*–*A. lyrata* genes that changes pooled `Dn` by 0.2% and `Ds` not at all, but `Dn` of 66% of genes, by up to 8 steps in one. Per-gene `Dn`, and so `alpha` and NI, can therefore move slightly when `snp_side` is swapped. Report `snp_side` with per-gene results.

## Divergence counted

Dn and Ds come from the same PAL2NAL codon alignments `orthologr::dNdS()` produced, one aligned codon at a time:

- Codons with a gap, a non-ACGT base, or a stop codon on either side are skipped.
- A codon that carries a SNP in the SNP-side species is excluded from the divergence counts, so a site is not counted as both polymorphic and divergent.
- Each difference is classified by the degeneracy of the position it changes, in the codon it starts from: a 4-fold position counts as a synonymous difference (`Ds`), a 0-fold position as nonsynonymous (`Dn`), and a 2- or 3-fold position counts in neither class (`div_steps_midfold`). The limitation above applies.
- A codon pair that differs at two or three positions is averaged over every shortest mutational path, with equal weight: each path steps through the intermediate codons, and a path that passes through a stop codon is discarded.

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
- No SNP filtering of any kind: `FILTER`, quality, depth, allele count, and genotypes are ignored. Filter the VCF with `bcftools` first, as above.
- Only biallelic SNPs are used; multiallelic and indel records are ignored.
- 2- and 3-fold sites are excluded from `Pn`/`Ps` and `Dn`/`Ds`. This makes the classes unambiguous, at the cost of power: on 1,472 public genes, 31.5% of divergence steps are mid-fold, so more genes have `Ps = 0`, `Ds = 0`, and an undefined `alpha` or neutrality index.
- Per-gene `Dn` depends slightly on `snp_side`, as described under "SNPs used".
- Codon-level only: SNPs outside CDS regions are counted nowhere.
- 1:1 orthologs only. Genes in many-to-many pairs are dropped.
- With few substitutions the test has little power; a non-significant result is not evidence of neutrality.
