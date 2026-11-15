# Analyzing your own genomes

This assumes you have completed [getting started](getting-started.md). The steps are the same, but the commands below use a placeholder file, so they are not run automatically.

## 1. Check your inputs

You need a genome FASTA and a GFF3 annotation with `gene`, `mRNA`, and `CDS` features for each genome.

- FASTA sequence names must match the first column of the GFF3.
- Leave out organelles. They use other genetic codes, and a trans-spliced gene such as plastid rps12 can stop `extract_cds`. For a RefSeq GFF3: `awk -F'\t' '$3=="region" && /genome=(chloroplast|mitochondrion)/{drop[$1]=1} !($1 in drop)' in.gff > nuclear.gff`
- Use annotations of comparable quality. Errors in either one can look like selection.
- For enrichment, annotate both genomes with the same term database version. Different versions can give wrong results.
- For polyploids, consider splitting subgenomes first, so that best hits occur between matching subgenomes: `dndsr split_comparisons -c comparison_file.txt`.

## 2. Write a comparison file

Use one line per comparison, separated by spaces or tabs. Use double quotes around paths containing spaces.

```
sp1_v_sp2 data/sp1.fasta data/sp1.gff3 data/sp2.fasta data/sp2.gff3
```

Paths are relative to the folder where you run `dndsr`. If your data is elsewhere, reinstall with `--bind` (see [installation options](../../INSTALL.md)) and use absolute paths, or run from a parent folder that contains your data.

## 3. Run

On a cluster, submit `calculate_dnds` through your scheduler, because it can take many hours on whole genomes.

```
dndsr extract_cds -c comparison_file.txt -t 8
dndsr calculate_dnds -c comparison_file.txt -t 40
dndsr append_annotations -c comparison_file.txt -t 8
```

## 4. Check before you interpret

- **Look at exclusions** in `<comparison>/<genome>_CDS_qc.tsv`. Many exclusions usually mean an annotation or FASTA problem.
- **Count pairs.** A very low number of orthologous pairs compared with your gene count suggests mismatched inputs or very distant genomes.
- **Look at the dS distribution.** `calculate_dnds` reports how many pairs are undefined, below dS 0.02, above dS 0.8, and above 3 times the median dS. Very low dS makes ratios unstable. High dS brings saturation and wrong pairings. See [interpretation](../interpretation.md).
- **Count genes above 1.** `calculate_dnds` says so when there are fewer than 30, because enrichment then has little power.

## 5. Continue

```
dndsr ipr_enrichment -c comparison_file.txt -t 8
dndsr go_enrichment -c comparison_file.txt -t 8
dndsr term_enrichment -c comparison_file.txt -t 8
dndsr covariate_enrichment -c comparison_file.txt -t 8
dndsr dnds_ideogram -c comparison_file.txt
dndsr dnds_summary -c comparison_file.txt --regions-bed regions.bed
```

`ipr_enrichment` and `go_enrichment` correct `p_adj` once across all their tables by default. To correct each table separately, add `--adjust-scope per_type` to `ipr_enrichment` or `--adjust-scope per_ontology` to `go_enrichment`. To get `p_adj` with no correction from `go_enrichment`, add `--p-adjust none`. See [interpretation](../interpretation.md) before reporting terms.

Short genes and genes with small dS pass dN/dS > 1 far more often, and the first three commands can call a term significant because of the kind of gene it holds. `covariate_enrichment` adjusts for length and dS, and is the one to report a term list from: see [its calibration](../validation/enrichment-confounding.md). It writes `<comparison>_<q or s>_<GO or IPR>_covariate_enrichment.tsv`. It reads `GO:` and `IPR` terms from the annotated table, and needs v0.4.0 or later. Its evidence comes from simulation at 4 plant pairs, and it has not been validated for detecting real selection, so its terms remain leads. `go_enrichment` and `term_enrichment` give a calibration warning on every run, and record it in a `_notice.json` file beside their tables.

To contrast comparisons that share a genome, make a contrast file with `name comparison_A side_A comparison_B side_B`, and run:

```
dndsr dnds_contrast -c comparison_file.txt --contrast-file contrast_file.txt
```

`dnds_summary` and `dnds_contrast` read `--regions-bed` as BED by default: 0-based start, end not included. Use `--regions-coord gff1` with `dnds_contrast` for 1-based inclusive coordinates.

## 6. McDonald-Kreitman test (development version)

To ask whether selection shapes a gene's polymorphism as well as its divergence, add a within-species SNP set. It must use the coordinates of one of the two genomes in the comparison: the genome FASTA and GFF3 for that side, already in the comparison file, are used to place each SNP in its codon.

**Filter the SNP calls first.** dndsR applies no filtering: it ignores `FILTER`, quality, depth, allele count, and genotypes, so a fixed site or a low-quality call is counted as a polymorphism. Use `bcftools`:

```
bcftools +fill-tags calls.vcf.gz -- -t AC,AN -Oz -o calls_tagged.vcf.gz   # only if AC/AN are missing
bcftools view -f PASS -v snps -m2 -M2 -i 'AC>=2 && AC<AN' calls_tagged.vcf.gz -Oz -o snps.vcf.gz
```

`AC>=2` suits a population VCF; one diploid sample has at most `AC = 2`, so it keeps only heterozygous sites. Filter BED input the same way before converting it.

Add two columns to the comparison file: `snp_file`, a VCF (`.vcf`, `.vcf.gz`) or BED (`.bed`, `.bed.gz`), and `snp_side`, `query` or `subject`. In a headerless file they are columns 6 and 7.

```
sp1_v_sp2 data/sp1.fasta data/sp1.gff3 data/sp2.fasta data/sp2.gff3 data/sp1_snps.vcf query
```

Then run the same commands as before; `calculate_dnds` writes two extra tables.

```
dndsr calculate_dnds -c comparison_file.txt -t 40
```

- `<comparison>_mk.tsv`: one row per 1:1 pair, with `pn`, `ps`, `dn`, `ds`, `alpha`, `neutrality_index`, `dos`, `fisher_p`, `fisher_q`.
- `<comparison>_mk_pooled.tsv`: the same statistics from the summed counts, one row for the whole comparison.

Only bi-allelic SNPs are used: multiallelic records, indels, and genotype columns are ignored. Each site is classified by codon degeneracy, the same way for VCF and BED: 4-fold positions count as synonymous, 0-fold as nonsynonymous, and 2- and 3-fold positions are excluded from both classes (they appear as `snps_midfold` and `div_steps_midfold`). Codons carrying a SNP are left out of the divergence counts. dN/dS cached from unchanged inputs skips the test; add `--overwrite` to recompute both.

Before trusting a table, check the run's messages. dndsR stops when no contig of the SNP file matches the chosen genome, when a VCF `##contig` length disagrees with the FASTA, or when more than 10% of REF alleles disagree with the genome. Those checks catch a SNP file from the wrong assembly, but not one that shares contig names and lengths, so confirm the source. Method, counting rules, and limits: [McDonald-Kreitman](../mcdonald-kreitman.md).

## Re-running

Existing results are kept unless you add `--overwrite`. dndsR tells you when an existing result no longer matches its inputs, using content hashes rather than file dates. Results record their software versions and parameters in `<comparison>_dnds_provenance.json`. Keep that file with your results and report it in your methods.
