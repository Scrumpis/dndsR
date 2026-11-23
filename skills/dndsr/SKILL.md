---
name: dndsr
description: Installs and runs dndsR, an R package and CLI for pairwise dN/dS analysis of orthologous genes from genome FASTA and GFF3 files, with annotation, term enrichment, ideograms, and regional or contrast summaries; from v0.5.0, experimental HyPhy BUSTED selection tests across species and an experimental McDonald-Kreitman test from a VCF or BED of SNP positions. Use when a user wants to calculate, summarize, or interpret dN/dS ratios, BUSTED tests, or McDonald-Kreitman tables with dndsR, or set it up with Singularity, Apptainer, Docker, or R.
license: MIT
compatibility: Needs bash and git, plus Singularity, Apptainer, or Docker, or an R installation with dndsR's dependencies. Only Linux amd64 with Singularity has been tested. Cloning and pulling the image need network access unless both are already present.
metadata:
  software-version: "0.5.0"
  skill-status: "core workflow tested on a designed example and public genomes; see references/verification.md"
---

# dndsR

dndsR estimates and analyzes dN/dS between orthologous genes in pairwise genome comparisons. It is both an R package and a command-line tool.

What has been tested is recorded in [references/verification.md](references/verification.md). Tell the user when you rely on something that has not been tested. Open problems are in [references/gaps.md](references/gaps.md).

## 1. Check the environment

```
command -v apptainer singularity docker Rscript
command -v dndsr && dndsr --help
```

If `dndsr --help` lists commands, go to step 3. Without containers, use the tested conda route in [INSTALL.md](../../INSTALL.md). If R has dndsR (`Rscript -e 'packageVersion("dndsR")'`), every command below is also a function: see [Using dndsR from R](../../docs/tutorials/using-r.md).

## 2. Install

Tell the user first: this downloads 1.3 GB, and writes to `~/.config/dndsr`, `~/.dndsr/rlib`, and `~/.local/bin`. The image holds dependencies; dndsR is installed from the clone, which must stay in place. Replace `singularity` with `apptainer` if that is what is available (untested).

```
git clone --branch v0.5.0 https://github.com/Scrumpis/dndsR
cd dndsR
singularity pull dndsr.sif docker://scrumpis/dndsr@sha256:f1b42d9d338b406ec23daca1bb146ed0329db29d21a24a6871158996aa0572da
tools/dndsr-install install --engine singularity --sif dndsr.sif
export PATH="$HOME/.local/bin:$PATH"
```

Until the v0.5.0 image digest replaces it, the digest above is the v0.4.0 image, which lacks HyPhy and KaKs_Calculator: tell the user so before `phylo_busted` or another `--dnds-method`. Pin digests rather than tags, because tags can move. To see published tags: `curl -s "https://hub.docker.com/v2/repositories/scrumpis/dndsr/tags?page_size=20"`. The Docker route is in [INSTALL.md](../../INSTALL.md) and is untested.

## 3. Verify with the example

Always do this after installing. It takes about a minute.

```
dndsr example
cd dndsr_example
dndsr extract_cds -c comparison_file.txt
dndsr calculate_dnds -c comparison_file.txt -t 2
dndsr verify_example
```

Success ends with `All 14 example checks passed.` Any `[FAIL]` line, or a nonzero exit code, means the installation does not behave as expected: report it rather than continuing. This confirms the installation works. It does not validate results on other data.

## 4. Run a real analysis

Input is a comparison file, one line per comparison: `name query.fasta query.gff3 subject.fasta subject.gff3`. Paths are relative to the current folder, which is visible in the container. For data elsewhere, reinstall with `--bind <folder>` and use absolute paths.

```
dndsr extract_cds -c comparison_file.txt -t 8
dndsr calculate_dnds -c comparison_file.txt -t 40
dndsr append_annotations -c comparison_file.txt -t 8
dndsr term_enrichment -c comparison_file.txt -t 8
dndsr ipr_enrichment -c comparison_file.txt -t 8
dndsr go_enrichment -c comparison_file.txt -t 8
dndsr covariate_enrichment -c comparison_file.txt -t 8
dndsr dnds_ideogram -c comparison_file.txt
dndsr dnds_summary -c comparison_file.txt --regions-bed regions.bed
```

If the comparison file lists the same two genomes twice, in either order, dndsR v0.5.0 and later can skip the repeats: pass `--unique-pairs` to every step that reads the file, or later steps stop on the skipped row's missing output. A repeated row with its own `snp_file` is kept. Details: [docs/tutorials/your-own-data.md](../../docs/tutorials/your-own-data.md).

`calculate_dnds` can take many hours on whole genomes, so on a cluster use the scheduler. Commands and options: [references/commands.md](references/commands.md). Output columns: [docs/output-files.md](../../docs/output-files.md).

Check `<comparison>/<genome>_CDS_qc.tsv` for exclusions, and the number of pairs in `<comparison>_dnds.tsv`, before interpreting results. `calculate_dnds` reports undefined, low-dS, high-dS, and outlier pairs: pass that on to the user.

From v0.5.0, dndsR also runs an experimental McDonald-Kreitman test when a comparison row adds `snp_file` (a VCF or BED of SNP positions for one side's genome) and `snp_side` (`query` or `subject`). It writes `<comparison>_mk.tsv` and `<comparison>_mk_pooled.tsv` beside the dN/dS table. It is opt-in, so the default outputs do not change. Only bi-allelic SNPs are used, and it needs that side's GFF3 and FASTA in the comparison, so it cannot run in `seq_list_file` mode. dN/dS cached from unchanged inputs skips it: rerun with `--overwrite`. It does no SNP filtering: filter the calls with `bcftools` before the run (PASS, bi-allelic, minimum allele count, and `AC < AN`), or a fixed site will be counted as a polymorphism.

`covariate_enrichment` needs v0.4.0 or later, so check that `dndsr --help` lists it. `go_enrichment` and `term_enrichment` give a calibration warning on every run, and record it in `<comparison>_go_enrichment_notice.json` or `<comparison>_term_enrichment_notice.json`: pass it on to the user. Their defaults and tables are unchanged.

A complete run on public genomes, with measured time and memory: [docs/tutorials/public-genomes.md](../../docs/tutorials/public-genomes.md). For a publication, follow [docs/publication-protocol.md](../../docs/publication-protocol.md).

## 4b. Selection tests across species (experimental, v0.5.0 and later)

`phylo_prepare` and `phylo_busted` go from genomes and orthogroups to HyPhy BUSTED tests: does any site show dN/dS above 1 on chosen foreground branches? They are new in v0.5.0. Use them only if `dndsr --help` lists `phylo_prepare`, and tell the user they are experimental, and calibrated only in simulation on one five-taxon tree. They need MAFFT (in the dndsR image) and HyPhy 2.5.103, which published images up to v0.4.0 lack; the tested route installs HyPhy with micromamba and binds it into the container. Steps: [docs/tutorials/phylogenetic-selection.md](../../docs/tutorials/phylogenetic-selection.md).

```
dndsr phylo_example
cd dndsr_phylo_example
dndsr phylo_prepare --genome-file genomes.tsv --orthogroup-file orthogroup_table.tsv --species-tree-file species_tree.nwk --foreground-file foreground.tsv -t 4
dndsr phylo_busted --prepared-dir prepared --hyphy "$HOME/.dndsr/hyphy/bin/hyphy" -t 4
dndsr verify_phylo_example
```

Success ends with `All 11 phylogenetic example checks passed.` On real data the user supplies a genome table (`taxon`, `fasta`, `gff`) or CDS table (`taxon`, `cds`), orthogroups (`orthogroup`, `taxon`, `gene_id`, or OrthoFinder's `Orthogroups.tsv` or `N0.tsv` with gene names matching the CDS), a rooted species tree, and a foreground file (`terminal`, `stem`, or `clade`), decided before seeing results. Do not build alignments or HyPhy inputs by hand. When reporting:

- A significant BUSTED result means some site on some foreground branch had ω above 1 under the model. It is not the same as, and does not confirm, a pairwise dN/dS above 1. Never mix the two in one table or claim.
- Report `p_best_fits_adj` with three starts, or `p_adj` with one, and `family_size`: only orthogroups with status `tested` are corrected together. Report counts of `excluded` (by `reason` and `exclusion_basis`) and `failed` orthogroups, which limit what the results can say. `failed` means a step errored: investigate, do not ignore.
- `foreground_branch_zero_length` (policy) means HyPhy deleted a requested branch as zero-length; it is not a failure.
- `foreground_unidentifiable` (mathematical) means the requested branches cannot be separated from other history with the taxa present; `foreground_partly_unidentifiable` and `foreground_incomplete` are conservative policy. dndsR never relabels a branch.
- Gene-tree discordance is reported (`gt_*` columns, `gene_tree_discordant_near_foreground`), not corrected.
- Results depend on `--seed`, `--hyphy-cpu`, and slightly on the processor (with or without AVX-512); the first two are recorded in `busted/busted_provenance.json`, so record the processor too. `--extra-seeds` shows the spread; it never changes `p_value`. Genes near the threshold can change call with the seed.
- The defaults are HyPhy's. Simulation recommends `--multiple-hits Double` and `--extra-seeds 2,3`, reading `lrt_best_fits`, `p_best_fits`, and its correction `p_best_fits_adj`; `p_adj` corrects only `p_value`, from the first start. Tell the user, and do not change settings without their choice. `--error-sink Yes` (BUSTED-E) is an experimental conservative check with a large power cost. Alignment errors cause most false positives: point out `dissimilar_window` flags. Masking them is a sensitivity analysis, not a dndsR feature. If a user asks what masking costs: in the mixtures, the lower power after correction came from BH once invalid discoveries were removed, not from masked signal; masking does cost real power when selected sites form one block (about 20 points in simulation). Evidence and limits: [docs/phylo-methods.md](../../docs/phylo-methods.md#calibration-in-simulation).
- Methods and limits: [docs/phylo-methods.md](../../docs/phylo-methods.md).

## 4c. Lineage-level selection (experimental, development code)

To compare selection among lineages (for example against LTR activity or genome-size change), four commands run on the same `phylo_prepare` folder. Use them only if `dndsr --help` lists `phylo_branch_omega`; tell the user they are experimental and calibrated only in simulation. They need HyPhy 2.5.103 and PAML 4.10.10 (`--hyphy`, `--codeml`, or on `PATH`). Steps, with the example: [docs/tutorials/lineage-selection.md](../../docs/tutorials/lineage-selection.md); methods: [docs/lineage-selection.md](../../docs/lineage-selection.md).

```
dndsr phylo_branch_omega --prepared-dir prepared -t 4
dndsr phylo_absrel --prepared-dir prepared -t 4
dndsr phylo_relax --prepared-dir prepared -t 4 --models Minimal
dndsr phylo_summarize --branch-omega-dir branch_omega --absrel-dir absrel --relax-dir relax
```

- They answer different questions: branch-model omega is gene-wide dN/dS per branch (not a test of positive selection; omega above 1 on a branch is not proof); aBSREL tests episodic positive selection per branch; RELAX tests relaxation or intensification, not positive selection. Do not merge their results into one claim.
- `pooled_omega`'s interval undercovers with few genes (about 91% at 25), and under site-to-site omega variation the pooled omega was 4-11% low, more on long branches (exploratory, `docs/dev/validation-audit.md`): tell the user, and suggest branch length as a covariate in any trait model.
- For lineages, report `pooled_omega` with its interval, `n_pooled`, and the flag fractions from `summary_lineages.tsv`, and aBSREL counts by `p_bh_lineage`, stating that it corrects within one lineage only (counts summed over lineages are exploratory) and assumes correct alignments; for gene-level aBSREL results report `p_adj`, not `p_value`; never the mean of per-gene omega. Summaries are descriptive: a lineage difference is not evidence of adaptation or of a trait association. For trait hypotheses, point the user to a phylogenetic comparative model (`nlme::gls` with `ape::corBrownian`, `phylolm`, `caper`), using `summary_branch_ids.tsv` and `summary_species_tree_ids.nwk`.
- Calibration (simulation only, `docs/validation/phylo-lineage-calibration.md`): aBSREL's tests held their level on error-free alignments; RELAX's `p_value` was anticonservative on 16 taxa (8-9% false positives at 0.05, reproduced); `p_fg_omega_gt1` at omega 1 (5.5-6.8%) and `p_free_vs_m0` with missing taxa (8%, then 4%) are not established as calibrated. Results carry these in the `calibration` column and a run warning: tell the user, and do not present a single gene's RELAX, free-ratio, or omega-above-1 result as established. A single free-ratio omega is noisy even when `reliable_free`.
- Input errors (simulation, `docs/lineage-selection.md#input-errors`): after realignment of sequences with indels, aBSREL's `p_raw`, `p_value`, `p_adj`, and `p_bh_lineage` all made excess false discoveries, and `dissimilar_window` missed those errors; the branch models were unaffected and BUSTED nearly so (8 and 16 taxa). A corrupted segment in one sequence breaks aBSREL, BUSTED (foreground sequence), and the branch models, and is flagged. When the run warns `Alignment check`, tell the user how many flagged genes are significant and to review them; never say unflagged genes are correctly aligned. Alignments a user supplies are accepted if in frame; in-frame is not evidence of correct homology, so ask which aligner made them and report it. Masking is an optional sensitivity analysis, not protection against false positives.
- Branches merged by pruning, including the two at the root, are left out of lineage summaries: advise an outgroup, and plainly labelled internal nodes in the species tree.
- Only single-copy orthogroups are tested (`multi_copy` exclusions), so gene families that duplicate or are lost are excluded, more so in lineages with more duplication. When the question is about genome evolution, say so, and report the number excluded as `multi_copy` beside any lineage result.

## 5. Interpret honestly

Read [docs/interpretation.md](../../docs/interpretation.md) before explaining results. The essentials:

- **dN/dS above 1 is a candidate, not proof of positive selection.** dndsR applies a threshold to a point estimate. In simulations with realistic codon usage, about half of truly neutral genes exceed 1.
- Look at dN and dS as well as the ratio. Small dS makes ratios unstable, and dS = 0 gives NA. Above dS 0.8, estimates of high ratios fall far short and ortholog pairs are often wrong: [docs/validation/estimators.md](../../docs/validation/estimators.md), [docs/validation/orthology.md](../../docs/validation/orthology.md).
- From v0.5.0, `calculate_dnds` also writes `<comparison>_dnds_reliability.tsv`: flags per pair (few codons, few synonymous differences, dS far below the observed synonymous differences, dN/dS above 10, and low, high, or outlying dS). Mention flagged genes when you discuss them. The flags change no value, and an unflagged estimate is not proven reliable: [docs/output-files.md](../../docs/output-files.md).
- McDonald-Kreitman tables are counts, not dN/dS estimates. `alpha` and the neutrality index compare the polymorphism and divergence ratios; interpret a significant Fisher `p`/`fisher_q` as evidence against neutrality at the gene, and read `dos` alongside because slightly deleterious polymorphisms can bias MK. Only 1:1 pairs are counted, biallelic SNPs only, allele frequencies are not used, and codons carrying a SNP are excluded from divergence. Each site is classed by its allele, as in the original test: a change is synonymous when the amino acid is unchanged, nonsynonymous when it changes, so 2- and 3-fold positions count too. A BED file must carry the ALT allele in column 4; without it the file is refused. The package does no SNP filtering, so say what filter was applied before the run. Report `n_pairs` and the pooled row, and say whether `snp_side` was `query` or `subject`. For several species, give one row per unordered pair (`k` species need `k(k-1)/2` rows) and run with `--unique-pairs` if the file also lists the reversed pairs, so no genome pair is computed twice. It has been tested on designed examples only, not on real population data: [references/verification.md](references/verification.md).
- Other methods (`--dnds-method NG`, `YN`, and others; v0.5.0 and later, needs KaKs_Calculator) use the same orthologs and alignments. Methods differ by a few percent per gene, and some differ systematically from PAML's versions: [docs/validation/estimators.md](../../docs/validation/estimators.md). Do not mix methods within one analysis, and say which one was used.
- The default Comeron method compresses ratios: very low ones read high, and high ones read low, more so at high divergence. Evidence: [docs/validation/estimators.md](../../docs/validation/estimators.md). Agreement with PAML and KaKs_Calculator on real genes: [docs/validation/public-genome.md](../../docs/validation/public-genome.md).
- Reciprocal best hits can mispair paralogs, especially in polyploids. Annotation errors can look like selection.
- `go_enrichment`, `ipr_enrichment`, and `term_enrichment` treat genes as independent, and do not correct for gene length, dS, or annotation completeness, which all affect which genes pass dN/dS > 1. Under simulated evolution where annotation has no link to omega, they can report false terms in up to about half of runs: [docs/validation/enrichment-confounding.md](../../docs/validation/enrichment-confounding.md). That is the chance of any false term, and with every term null also the false discovery rate. The cause is confounding by length and dS, not the BH correction, which held when the call did not depend on them ([study N](../../docs/validation/enrichment-error-calibration.md)); say this when a user asks whether the correction failed. For a term list to report, prefer `covariate_enrichment`, which adjusts for length and dS, and say which you used. It does not adjust for annotation completeness or dependence between genes. Its calibration comes from simulation at 4 plant pairs: it does not establish false discovery control in general, and it has not been validated for detecting real selection in real genomes. Polyploids are untested: run the simulation on the user's own table. Terms near `p_adj` 0.05 change with the annotation release or a filter, so report whether each term you discuss is stable.
- Few genes above 1 means little power. In the public comparisons, genes above 1 fell from 206 at median dS 0.14 to 3 at median dS 0.44. A short list at high divergence is not evidence of absence.
- Multiple-testing correction limits false positives under no association. It does not remove gene length, dS, or annotation bias. A significant term is a lead, not evidence that a function evolved under selection. Report it that way.
- Report `p_adj`, not `p_value` or `raw_p`, which are uncorrected. `p_adj` is corrected across `n_tested` terms, and `min_pos` changes only which rows are written. By default, `ipr_enrichment` and `go_enrichment` use one family across all their tables. `--adjust-scope per_type` (`ipr_enrichment`) and `--adjust-scope per_ontology` (`go_enrichment`) correct each table separately, so do not treat terms pooled across those tables as corrected together; `--p-adjust none` (`go_enrichment`) and `--fdr-method none` (`ipr_enrichment`) leave `p_adj` uncorrected. Check `adjust_scope` and `p_adjust_method` in the table.
- BH does not guarantee nominal error control under all conditions, and significant terms remain exploratory, not validated evidence of selection. A corrected GO term can still come from confounding. Results from before this default change had uncorrected GO `p_adj`, and per-type InterPro `p_adj`.
- `dnds_state_contrast --level gene` is disabled, because its results were not meaningful.
- Results from dndsR before 0.3.0 should not be reused. They can contain false dN/dS above 1. See [docs/migrating-to-0.3.md](../../docs/migrating-to-0.3.md).

Report the version, commit, image digest, command, and the number of excluded genes. Parameters and input hashes are in `<comparison>_dnds_provenance.json`.

## 6. Troubleshooting

- `command not found: dndsr`: add `~/.local/bin` to `PATH`, then run `hash -r`.
- `Unknown command` or `Unknown option`: dndsR stops on typos. Run `dndsr --help` or `dndsr <command> --help`.
- `Sequence files not found`: run `extract_cds` first, with the same `-o`.
- `Some sequence names ... cannot be used` (v0.5.0 and later): GFF3 IDs, which become sequence names, contain whitespace or start with an NCBI tag such as `gnl|` or `ref|`, which DIAMOND rewrites. Rename them in the GFF3, with their `Parent` references, or use `--group-by gene` if only transcript IDs are affected. Evidence: [docs/validation/known-issues.md](../../docs/validation/known-issues.md).
- `jobs failed`: the message gives each error, and logs are in the comparison folder. Retry with `-t 1`, and check memory if a worker died.
- Enrichment writes no table: a small dataset may lose terms to `--max-prop` (default 0.2), `--min-pos`, or `--min-total`. Read the log before relaxing them.
- Singularity cannot start, for example without `/dev/fuse`: try Apptainer or the conda route. Report what you changed.
- `dndsr doctor` and `tools/dndsr-install doctor` check the installation.

## Citation

No DOI exists yet. Use [CITATION.cff](../../CITATION.cff) with the exact commit, and also cite orthologr, DIAMOND, and Comeron (1995).
