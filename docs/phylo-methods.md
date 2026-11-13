# Phylogenetic selection tests: methods

`phylo_prepare` and `phylo_busted` test orthogroups of several species for episodic positive selection on chosen branches, using HyPhy's BUSTED. This is separate from dndsR's pairwise dN/dS: different inputs, a different question, and its own outputs. Tutorial: [tutorials/phylogenetic-selection.md](tutorials/phylogenetic-selection.md).

**Status: experimental and unreleased** (development version 0.4.0.9000; not in v0.4.0). dndsR prepares inputs as designed and reproduces HyPhy run directly, exactly, on synthetic examples ([wrapper](validation/phylo-busted-wrapper.md), [preparation](validation/phylo-prepare.md)). That is engineering correctness. The pipeline has **not** been calibrated as a selection test: its false-positive rate and power under realistic conditions are not measured. Treat results as you would direct HyPhy results, with the caveats below.

## What BUSTED tests

BUSTED (Murrell et al. 2015, *Mol Biol Evol* 32:1365) fits a codon model in which each site has one of three dN/dS (ω) classes on the foreground branches, and asks whether the class allowed above 1 is needed. The null fixes it at 1. A significant result means: **at least one site, on at least one foreground branch, evolved with ω above 1.** It does not say which sites or branches, and it does not mean the gene as a whole evolved under positive selection.

By default, synonymous rates vary among sites (BUSTED-S, Wisotsky et al. 2020, *Mol Biol Evol* 37:2430). dndsR keeps every HyPhy default, including no modelling of multinucleotide substitutions (`--multiple-hits None`). Alternatives are available (`--multiple-hits Double` or `Double+Triple`; `--error-sink Yes` for BUSTED-E, which HyPhy and dndsR treat as experimental). Which should be recommended is left to calibration.

| | Pairwise dN/dS (`calculate_dnds`) | BUSTED (`phylo_busted`) |
|---|---|---|
| Input | Two genomes | One codon alignment per orthogroup, a species tree, a foreground |
| Output | A point estimate per gene pair | A likelihood-ratio test per orthogroup |
| A "positive" | dN/dS above 1, a threshold on an estimate | p_adj below your chosen level |
| Uncertainty | None attached | A p-value, corrected across orthogroups |

A significant BUSTED result and a pairwise dN/dS above 1 are not interchangeable, and one does not confirm the other.

## Workflow

1. **`phylo_prepare`** builds codon alignments and trees from genomes or CDS and an orthogroup table.
2. **`phylo_busted --prepared-dir`** runs BUSTED on them and corrects across orthogroups.

Codon alignments made elsewhere can go straight to `phylo_busted --alignment-dir`.

## Inputs

- **Genomes or CDS.** A genome table (`taxon`, `fasta`, `gff`) whose CDS dndsR extracts with `extract_cds` (longest valid isoform per gene, named by gene ID), or a CDS table (`taxon`, `cds`). Either may have an `orthofinder_name` column.
- **Orthogroups.** dndsR's long table (`orthogroup`, `taxon`, `gene_id`), OrthoFinder's `Orthogroups.tsv`, or OrthoFinder's `N0.tsv` (named by HOG). OrthoFinder's species columns must equal a taxon or its `orthofinder_name`, and its gene names must equal the CDS names, as when OrthoFinder is run on proteins from `extract_cds --export-proteins`.
- **A rooted species tree** in Newick. Taxon names use letters, digits, and underscores, starting with a letter. Branch lengths are ignored: HyPhy receives the topology only.
- **A foreground file**, tab-separated, decided before you see results:

  ```
  scope	taxa
  clade	spA,spB
  ```

  | Scope | Branches labelled |
  |---|---|
  | `terminal` | The terminal branch of each listed taxon |
  | `stem` | The branch leading to the clade formed by the listed taxa |
  | `clade` | That stem and every branch within the clade |

  Rows are combined into one hypothesis. Stem and clade taxa must form a clade in the species tree.
- **Optional gene trees**, `<orthogroup>.nwk` or OrthoFinder's `<orthogroup>_tree.txt`. Tips may be taxa, gene IDs, `taxon_geneID`, or OrthoFinder's `<proteome>_<geneID>`.

## Preparation

For each orthogroup, in order, stopping at the first exclusion:

1. **Single copy.** More than one gene from a taxon excludes the orthogroup (`multi_copy`). All genes stay in the member table, so duplication-aware analyses can be added later.
2. **CDS.** Every gene must have a CDS. A gene that `extract_cds` rejected excludes the orthogroup (`cds_invalid`, with its reason: length not a multiple of 3, or an in-frame stop); a gene with no CDS at all gives `gene_not_found`. One bad member excludes the whole orthogroup: nothing is dropped silently.
3. **Sequence checks.** Characters must be nucleotides or IUPAC codes, lengths multiples of 3, and stop codons only at the end. A final stop codon is removed and counted (`--terminal-stops exclude` excludes the orthogroup instead). Ambiguous codons and a missing start codon are reported, not changed.
4. **Minimum taxa** (`--min-taxa`, default 4).
5. **Trees and foreground**, below, before any alignment.
6. **Alignment.** Proteins are aligned with MAFFT (L-INS-i by default, one thread), and each residue is replaced by its own codon. Each aligned protein, without gaps, must equal the translation, and the codon alignment, without gaps, must equal the coding sequences exactly; otherwise the orthogroup fails (`backtranslation_mismatch`). Nothing is trimmed or masked.

## Outcomes

Every orthogroup ends with one status, and the reason and its basis:

| Status | Meaning |
|---|---|
| `tested` | HyPhy ran, and its output was read and checked |
| `excluded` | A rule kept it out |
| `failed` | A step errored: MAFFT, HyPhy, or a consistency check. Not in the correction family; investigate |

| `exclusion_basis` | Meaning | Reasons |
|---|---|---|
| `mathematical` | The requested branches cannot be separated from background history, or do not exist, in this orthogroup's topology | `foreground_unidentifiable`, `foreground_not_in_gene_tree` |
| `no_data` | No foreground taxon is present | `foreground_absent` |
| `policy` | Testable in principle, but excluded by a conservative rule | `multi_copy`, `too_few_taxa`, `foreground_partly_unidentifiable`, `foreground_incomplete` |
| `data_quality` | The sequences fail checks | `cds_invalid`, `internal_stop`, `length_not_multiple_of_3`, `invalid_character`, `split_codon_gap`, `empty_sequence`, `terminal_stop` |
| `input` | The inputs do not match | `gene_not_found`, `unknown_taxon`, `not_aligned`, `empty_alignment`, `gene_tree_missing`, `gene_tree_unreadable`, `gene_tree_taxa_mismatch` |
| `software` | (failed) | `alignment_failed`, `backtranslation_mismatch`, `hyphy_error`, `hyphy_output_invalid`, `label_mismatch`, `prepared_tree_mismatch` |

`stage` says whether it happened in `prepare` or `analysis`. `phylo_busted --prepared-dir` lists preparation's exclusions too, so one table accounts for every orthogroup.

## Missing taxa and identifiability

Branches are identified by the split of taxa they make, never by name. The species tree is pruned to the taxa present, and branches that then join become one branch. The hypothesis is never moved onto a different branch:

- **Mathematical.** If every requested branch present has joined a background branch, the requested history cannot be separated from background history: `foreground_unidentifiable`. HyPhy treats trees as unrooted, so the stem of a clade beside the root is always unidentifiable.
- **Policy.** If some requested branches remain separable and others have joined background branches, testing only the separable ones would test a narrower hypothesis: `foreground_partly_unidentifiable`. If a hypothesis row has no taxon present, or a terminal row lacks a listed taxon: `foreground_incomplete`.
- **Allowed and counted.** Missing members of a clade-scope row: the clade's remaining branches carry the hypothesis, and `fg_absent_branches` and the flag `foreground_branches_absent:<n>` count what is missing.

In practice: `stem` needs every member of the clade and at least one taxon of its sister lineage; `clade` tolerates missing members but needs the sister lineage; `terminal` needs the taxon's sister lineage. Taxa missing elsewhere do not matter.

## Gene trees and discordance

With gene trees, discordance with the pruned species tree is always measured: `gt_missing_species_splits` (species-tree branches absent from the gene tree), `gt_extra_splits`, `gt_missing_splits` (which), and `gt_discordant_near_foreground` (whether a missing branch is a foreground branch or touches one). Flags `gene_tree_discordant` and `gene_tree_discordant_near_foreground` are added. **Discordance never excludes an orthogroup.**

- `--gene-tree-mode topology` (default): HyPhy uses the gene tree. Requested terminal and stem branches, including each clade's stem, must exist in it, or the hypothesis cannot be placed (`foreground_not_in_gene_tree`, mathematical). Branches inside a clade-scope foreground are foreground.
- `--gene-tree-mode report`: HyPhy uses the species tree, and the gene tree is only compared.

An orthogroup without a gene tree, when gene trees are given, is excluded (`gene_tree_missing`) rather than mixed with species-tree orthogroups.

## Statistics

- `p_value` is HyPhy's: a 50:50 mixture of chi-square distributions with 0 and 2 degrees of freedom. When the likelihood ratio is 0, p is exactly 0.5, so many genes share that value.
- `p_adj` is Benjamini-Hochberg across every orthogroup tested in the run. Excluded and failed orthogroups are not in the family. Compare `family_size` with the number you started with.
- Each run tests one foreground hypothesis. Separate runs are separate families.

## Reproducibility and multiple starts

HyPhy optimizes from random starting values, and its result depends on them and on its thread count. dndsR always sets both (`--seed`, default 1; `--hyphy-cpu`, default 1), and records them; the same HyPhy build on the same kind of processor then gives identical results. HyPhy's arithmetic differs slightly between processors with and without AVX-512: on the examples, identical on Intel Skylake and AMD Genoa (AVX-512), and on AMD Rome and Milan (AVX2), with two of eight runs differing between the two groups, by up to 0.005 in the likelihood ratio ([evidence](validation/phylo-prepare.md#processors)). Record the processor with your results. MAFFT runs on one thread. A different seed can change a gene's likelihood ratio: on one simulated gene, from 7.0 to 8.8 over 20 seeds.

`--extra-seeds 2,3,4` runs independent further starts and reports, as diagnostics only: `n_starts`, `lrt_min_starts`, `lrt_max_starts`, the best alternative and null log-likelihoods found (`lnl_unconstrained_best`, `lnl_null_best`; where HyPhy skipped the constrained fit for lack of evidence, that start's unconstrained fit counts as a null point), the likelihood ratio between them (`lrt_best_fits`, `p_best_fits`), and whether the primary start found each best. **`p_value` and `p_adj` always come from `--seed`.** dndsR never chooses the largest likelihood ratio. Whether a multistart policy should change the reported test is for calibration, under rules fixed in advance.

## Alignment checks reported

These are reported in `qc_flags` and the sequences tables, but change nothing. Their thresholds are rules of thumb, not validated:

| Flag | Raised when |
|---|---|
| `ambiguous_codons:<taxon>` | Over 5% of the taxon's codons contain an ambiguous base |
| `gappy_sequence:<taxon>` | Over half the taxon's codons are gaps |
| `gappy_columns` | Over 30% of codon columns are gaps in more than half the sequences |
| `divergent_sequence:<taxon>` | The taxon's mean p-distance is over 3 times the median of all taxa's |
| `high_divergence` | A pairwise nucleotide p-distance exceeds 0.35 |
| `dissimilar_window:<taxon>` | In a window of 10 codons, at least 70% of the taxon's amino acids occur in no other sequence at that position: a sign of misalignment, frameshift, or a gene-model error |
| `no_start_codon:<taxon>` | The CDS does not start with ATG (preparation) |
| `terminal_stops_masked:<n>` | `phylo_busted` on supplied alignments masked n final stop codons |
| `foreground_branches_absent:<n>` | n requested branches have no taxa present (clade scope) |
| `gene_tree_discordant`, `gene_tree_discordant_near_foreground` | See above |
| `hyphy_refit_after_negative_lrt` | HyPhy reported a convergence problem and refitted |

Look at flagged alignments before interpreting their results. Alignment errors are a known cause of false positives in selection tests (for example, Markova-Raina and Petrov 2011, *Genome Res* 21:863; Jordan and Goldman 2012, *Mol Biol Evol* 29:1125).

## Limits and interpretation

- **Not calibrated.** Engineering correctness is established on synthetic data; false-positive rate and power of the whole pipeline are not.
- **Alignment errors, misassigned orthologs, gene-model errors, and recombination** can produce significant results. Recombination is not screened; alignments are not filtered.
- **Multinucleotide substitutions** can inflate positives; including them reduced BUSTED's positive calls by about a quarter across many empirical datasets (Lucaci et al. 2023, *Mol Biol Evol*). The default does not model them.
- **Gene tree discordance** (incomplete lineage sorting, introgression) can put substitutions on the wrong branches of a species tree and mimic branch-specific selection (Mendes and Hahn 2016, *Syst Biol* 65:711). It is reported, not corrected.
- **Ascertainment.** Only single-copy, identifiable, QC-passing orthogroups are tested. Duplicated genes, genes with problematic gene models, and genes missing from key taxa are under-represented, so conclusions apply to the tested set. The provenance files count orthogroups per taxon and exclusions by reason.
- A significant result is evidence for episodic positive selection under the model, not proof. Look at the alignment.

## Reporting

Report the dndsR version and commit, HyPhy and MAFFT versions, the settings, seed, and HyPhy threads (in the provenance files), the foreground file, the orthogroup source, the number of orthogroups tested, excluded by reason, and failed, and the correction family. Cite HyPhy (Kosakovsky Pond et al. 2020, *Mol Biol Evol* 37:295), BUSTED, BUSTED-S, MAFFT (Katoh and Standley 2013, *Mol Biol Evol* 30:772), and OrthoFinder if used.
