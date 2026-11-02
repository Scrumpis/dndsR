# Running dndsR for a publication

A checklist for an analysis you intend to publish. It fixes choices before you see results. Each step points to its evidence, and thresholds below are starting points to state in your methods, not guarantees.

## 1. Fix the software

Use one tagged release and one image digest for every comparison, and do not `git pull` during the analysis, because a pull reinstalls dndsR on the next run.

```
git clone --branch <tag> https://github.com/Scrumpis/dndsR
cd dndsR
singularity pull dndsr.sif docker://scrumpis/dndsr@sha256:<digest>
tools/dndsr-install install --engine singularity --sif dndsr.sif
dndsr --help
```

Run `dndsr example` and `dndsr verify_example` first. Record the tag, commit, and digest. `<comparison>_dnds_provenance.json` records versions, parameters, and input hashes. Keep it.

## 2. Check the inputs

For each genome:

- [ ] Hash the FASTA and GFF3, and record where they came from.
- [ ] Remove organelles: they use other genetic codes, and can stop `extract_cds`. See [your own data](tutorials/your-own-data.md).
- [ ] Count partial genes, which dndsR does not trim to frame: `python3 tools/check-first-exon-phase.py <gff3>`. Above a few in a thousand, look at them.
- [ ] Use annotations of comparable quality on each side, and the same term database version.
- [ ] For polyploids, decide whether to split subgenomes before you see results. Pairing against a duplicated genome disagreed with OMA far more often ([evidence](validation/orthology.md)).

Decide these before running, and write them down:

- [ ] Which comparisons, and which is query.
- [ ] Which terms, with the database version.
- [ ] The candidate rule. `covariate_enrichment` defaults to `--pos-threshold 1`, `--max-dnds 10`, `--min-total 10`, `--min-pos 2`. Do not change these after seeing term lists. If you want to explore them, report all of what you tried.

## 3. Run

```
dndsr extract_cds -c comparison_file.txt -t 8
dndsr calculate_dnds -c comparison_file.txt -t 32
dndsr append_annotations -c comparison_file.txt -t 8
dndsr covariate_enrichment -c comparison_file.txt -t 8
dndsr go_enrichment -c comparison_file.txt -t 8
dndsr ipr_enrichment -c comparison_file.txt -t 8
dndsr dnds_summary -c comparison_file.txt
```

Allow about 1 GB of memory per `calculate_dnds` thread. Start from a fresh session, and put each comparison through a scheduler.

## 4. Check each comparison

Fill this in for each one, and report it.

| Check | Where | What to look for |
|---|---|---|
| CDS written and excluded, per genome | `<comparison>/<genome>_CDS_qc.tsv`, or the `extract_cds` log | Many exclusions suggest an annotation or FASTA problem |
| Pairs, against genes in each genome | `<comparison>_dnds.tsv` | Few pairs suggest mismatched inputs or high divergence |
| Pairs undefined; below dS 0.02; above dS 0.8; above 3 x median dS; above dN/dS 1 | The `calculate_dnds` report | See below |
| Median dS and its range | Same | Compare with what is known of the divergence |
| Length and dS of genes above 1 | `<comparison>_dnds.tsv` | Compare with other genes. Short genes are expected to dominate |

How to read that report, from [estimator](validation/estimators.md) and [orthology](validation/orthology.md) evidence on public plant genomes:

| Finding | Meaning |
|---|---|
| Fewer than 30 genes above dN/dS 1 | Little power for enrichment. Reporting no terms says little |
| Median dS above about 0.4 | Few true ratios above 1 will pass 1, and high ratios read far too low. Using `dNdS > 1` is probably the wrong question here |
| Pairs above dS 0.8 | Ratios read low, and many pairs may not be orthologs |
| Pairs above 3 x median dS | A wrong pairing is likely at low divergence: 85% of such pairs disagreed with OMA, in one comparison |
| Many pairs below dS 0.02 | Ratios are unstable. Short genes at very low dS are often undefined |

## 5. Choose the term list

Report `covariate_enrichment`'s `p_adj` as the primary list. It is the one test here with simulation evidence of holding its false-term rate at both high and low numbers of positive genes, at `pos_threshold = 1` in 4 plant pairs, and it accounts for length and dS ([evidence](validation/enrichment-confounding.md)). Report `go_enrichment` and `ipr_enrichment` alongside, to show what the adjustment removes. Do not choose among them by which result you prefer.

Check each term you intend to discuss:

- [ ] `pos_count`, `total_count`, `expected`, and `obs_over_expected`, in the table. A term resting on 3 genes is fragile.
- [ ] **Sensitivity to doubtful pairs.** A few terms depended on pairs OMA does not support. Repeat without pairs above 3 x median dS, in a separate folder, and report which terms change:

```
dndsr covariate_enrichment --dnds-annot-file sensitivity/<comparison>_dnds_annot.tsv --sides query --filter-expr "dS <= 3 * median(dS, na.rm = TRUE)"
```


  The repeat is [calibrated](validation/enrichment-confounding.md#the-ds-filter-in-the-checklist), so its list is valid too, and neither run is more correct. Report the terms in both as stable. A term in only one, especially one resting on 2 or 3 positive genes near `p_adj` 0.05, is fragile. Say which you report as primary before you see either: here, the unfiltered run.
- [ ] **Do not rely on one family's genes being independent.** Tandem families share history. Look at whether a term's positive genes cluster on the genome: `q_gff_seqname` and `q_gff_start` in the annotated table.
- [ ] **For a polyploid or another clade, check calibration on your own table.** The evidence does not cover them. See the end of this page. Homeologs carry much the same annotation and share each pair's dN/dS, so do not count agreement between sides or between comparisons of one genome as replication.
- [ ] **Test by other means.** Use site or branch-site models with more than two sequences, on real alignments, before you claim selection on any gene or family.

## 6. What the results can and cannot support

**Can be said:** "genes with dN/dS above 1 hold more genes with term X than their length and dS predict, at this false discovery rate"; "these pairs have this dN/dS, by this estimator".

**Cannot be said from this alone:** that a gene or a family evolved under positive selection; that a function is the target of selection; that a term absent from the list is not under selection, especially at low numbers of positive genes; that comparisons of different divergence differ biologically. `dnds_contrast` is confounded by divergence time.

**Known limits to state:** one-to-one pairs from reciprocal best hits; a whole-gene, whole-branch average; a Comeron estimate that reads low at high ratios and high dS; no adjustment for annotation completeness, or for dependence among genes; calibration measured on public plant genomes, and not on yours.

## 7. Report

The version, commit, and image digest; input hashes; commands; the table from step 4; thresholds, and when they were chosen; the term lists with `n_tested` and all of what was tried. Cite dndsR ([CITATION.cff](../CITATION.cff)), orthologr, DIAMOND, and Comeron (1995), and GOseq only if you compare with it.

To check calibration for your own data, `tools/dev/validation/sim-genome.R` and `covariate-confirm.R` run the simulation on any annotated table. Consider this at a high or low number of positive genes, or for a clade far from the ones tested. See [the evidence](validation/enrichment-confounding.md#rebuild). Keep those outputs out of Git until published.
