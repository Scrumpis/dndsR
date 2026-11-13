# Phylogenetic selection tests: methods

`phylo_busted` tests orthogroups of several species for episodic positive selection on chosen branches, using HyPhy's BUSTED. It is separate from dndsR's pairwise dN/dS: different inputs, a different question, and its own outputs. A tutorial is in [tutorials/phylogenetic-selection.md](tutorials/phylogenetic-selection.md).

**Status.** dndsR runs BUSTED and reproduces HyPhy run directly, exactly, on the bundled example ([evidence](validation/phylo-busted-wrapper.md)). The pipeline has **not** been calibrated as a selection test: false-positive rate and power under realistic conditions are not yet measured. Treat results as you would direct HyPhy results, with the caveats below.

## What BUSTED tests

BUSTED (Murrell et al. 2015, *Mol Biol Evol* 32:1365) fits a codon model in which each site has one of three dN/dS (ω) classes on the foreground branches, and asks whether the class allowed above 1 is needed. The null fixes it at 1. A significant result means: **at least one site, on at least one foreground branch, evolved with ω above 1.** It does not say which sites or branches, and it does not mean the gene as a whole evolved under positive selection.

By default, synonymous rates vary among sites (BUSTED-S, Wisotsky et al. 2020, *Mol Biol Evol* 37:2430), which reduces false positives from synonymous rate variation. This is HyPhy's default, and dndsR keeps every HyPhy default.

| | Pairwise dN/dS (`calculate_dnds`) | BUSTED (`phylo_busted`) |
|---|---|---|
| Input | Two genomes | One codon alignment per orthogroup, a species tree, a foreground |
| Output | A point estimate per gene pair | A likelihood-ratio test per orthogroup |
| A "positive" | dN/dS above 1, a threshold on an estimate | p_adj below your chosen level |
| Uncertainty | None attached | A p-value, corrected across orthogroups |

A significant BUSTED result and a pairwise dN/dS above 1 are not interchangeable, and one does not confirm the other.

## Inputs

- **Codon alignments**, one per orthogroup, named `<orthogroup>.fasta`. Sequences are named `taxon` or `taxon|gene_id`. Phase 1 needs alignments you have already made; building them from `extract_cds` output is planned.
- **A rooted species tree** in Newick. Taxon names use letters, digits, and underscores, starting with a letter. Branch lengths are ignored.
- **A foreground file**, tab-separated:

  ```
  scope	taxa
  clade	spA,spB
  ```

  | Scope | Branches labelled |
  |---|---|
  | `terminal` | The terminal branch of each listed taxon |
  | `stem` | The branch leading to the clade formed by the listed taxa |
  | `clade` | That stem and every branch within the clade |

  Rows are combined into one foreground. Stem and clade taxa must form a clade in the species tree.
- **Optional gene trees**, `<orthogroup>.nwk`, used instead of the pruned species tree.

## Orthogroups and missing taxa

Only single-copy orthogroups are tested: one sequence per taxon. Others are skipped as `multi_copy`.

The species tree is pruned to the taxa present. When taxa are missing, neighbouring branches join. dndsR never moves a foreground label onto a joined branch that would also contain background history. Such an orthogroup is skipped as `foreground_not_identifiable`, and `detail` names the branches involved.

- **`stem`** needs every member of the clade, and at least one taxon of the clade's sister lineage.
- **`clade`** tolerates missing members, and needs at least one taxon of the clade's sister lineage. Without it, the stem would include time before the clade split from its sister.
- **`terminal`** needs at least one taxon of that taxon's sister lineage.
- Missing taxa elsewhere in the tree do not matter.
- If no foreground taxon is present, the orthogroup is skipped as `foreground_absent`.

HyPhy treats trees as unrooted, so the stem of a clade beside the root cannot be told apart from the other root branch. That hypothesis is never identifiable; choose a hypothesis on a branch away from the root, or add an outgroup.

With gene trees, terminal and stem foreground branches, and each clade's stem, must exist in the gene tree, or the orthogroup is skipped as `foreground_discordant`. Inside a clade-scope foreground, any gene-tree branch is foreground.

## Alignment checks

dndsR does not repair or mask sequences, with one exception: a stop codon ending a sequence is replaced by gaps and counted, because HyPhy refuses stop codons (`--terminal-stops exclude` skips the orthogroup instead). These exclude an orthogroup:

`empty_alignment`, `unknown_taxon`, `multi_copy`, `not_aligned`, `length_not_multiple_of_3`, `invalid_character`, `split_codon_gap` (a codon with one or two gap characters, usually a frameshift or misalignment), `internal_stop`, `empty_sequence`, `too_few_taxa`.

These are reported in `qc_flags` and `<run>_sequences.tsv`, but change nothing. Their thresholds are rules of thumb, not validated:

| Flag | Raised when |
|---|---|
| `ambiguous_codons:<taxon>` | Over 5% of the taxon's codons contain an ambiguous base |
| `gappy_sequence:<taxon>` | Over half the taxon's codons are gaps |
| `gappy_columns` | Over 30% of codon columns are gaps in more than half the sequences |
| `divergent_sequence:<taxon>` | The taxon's mean p-distance is over 3 times the median of all taxa's |
| `high_divergence` | A pairwise nucleotide p-distance exceeds 0.35 |
| `dissimilar_window:<taxon>` | In a window of 10 codons, at least 70% of the taxon's amino acids occur in no other sequence at that position: a sign of misalignment, frameshift, or a gene-model error |
| `terminal_stops_masked:<n>` | n terminal stop codons were masked |
| `hyphy_refit_after_negative_lrt` | HyPhy reported a convergence problem and refitted |

Look at flagged alignments before interpreting their results. Alignment errors are a known cause of false positives in selection tests (for example, Markova-Raina and Petrov 2011, *Genome Res* 21:863; Jordan and Goldman 2012, *Mol Biol Evol* 29:1125).

## Statistics

- `p_value` is HyPhy's: a 50:50 mixture of chi-square distributions with 0 and 2 degrees of freedom. When the likelihood ratio is 0, p is exactly 0.5, so many genes share that value.
- `p_adj` is Benjamini-Hochberg across every orthogroup tested in the run. Skipped and failed orthogroups are not in this family. Compare `family_size` with the number of orthogroups you started with.
- Each run tests one foreground hypothesis. Several runs with different foregrounds are separate families; correct across them yourself if you compare them.

## Reproducibility

HyPhy optimizes from random starting values, and its result depends on them and on its thread count. dndsR always sets both (`--seed`, default 1; `--hyphy-cpu`, default 1), and records them. With the same HyPhy build, these give identical results. A different seed can change a gene's likelihood ratio: on one simulated gene, from 7.0 to 8.8 over 20 seeds. Genes near your significance threshold may change call. The evidence is in the [validation page](validation/phylo-busted-wrapper.md).

## Limits and interpretation

- **Not yet calibrated.** Wrapper correctness is established; false-positive rate and power of the whole pipeline are not.
- **Alignment errors, misassigned orthologs, gene-model errors, and recombination** can produce significant results. Recombination is not screened.
- **Multinucleotide substitutions** can inflate positives; including them reduced BUSTED's positive calls by about a quarter across many empirical datasets (Lucaci et al. 2023, *Mol Biol Evol*). HyPhy's default, kept here, does not model them; `--multiple-hits Double` does.
- **Gene tree discordance** (incomplete lineage sorting, introgression) can put substitutions on the wrong branches of a species tree and mimic branch-specific selection (Mendes and Hahn 2016, *Syst Biol* 65:711).
- **Ascertainment.** Only single-copy, identifiable, QC-passing orthogroups are tested. Duplicated genes, fast-evolving genes that fail alignment checks, and genes missing from key taxa are under-represented, so conclusions apply to the tested set.
- A significant result is evidence for episodic positive selection under the model, not proof. Look at the alignment.

## Reporting

Report the dndsR version and commit, HyPhy version, the settings, seed, and HyPhy threads (all in `<run>_provenance.json`), the foreground file, the number of orthogroups tested and skipped by reason, and the correction family. Cite HyPhy (Kosakovsky Pond et al. 2020, *Mol Biol Evol* 37:295), BUSTED, and BUSTED-S.
