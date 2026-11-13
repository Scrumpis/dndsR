# phylo_prepare: input preparation, end to end

## Question

From genomes and orthogroups, does `phylo_prepare` build exactly the intended codon alignments and labelled trees, exclude what it should for the stated reasons, and never change the hypothesis? Do the results then equal HyPhy run directly, including on gene trees?

This is **engineering correctness**. It is not evidence that BUSTED's selection calls on this pipeline have acceptable false-positive rates or power: that is phase 3 calibration, not done. Wrapper evidence from phase 1: [phylo-busted-wrapper.md](phylo-busted-wrapper.md).

## Software

dndsR 0.4.0.9000 on `dev` after `d293a8e`, in the commit adding this page. HyPhy 2.5.103 (Bioconda `hyphy-2.5.103-h526e2cb_0`); MAFFT 7.526 and R 4.6.2 with ape 5.8.1 and Biostrings 2.80.2, from the v0.4.0 image (`sha256:f1b42d9d…`). The same runs in a development image built locally from `tools/Dockerfile` at `d293a8e` with Apptainer (`sha256:5d6807590eedf6cd5a409c8b05d2c5fc0e3c8ded3b8e526460f303c70de311a5`, not published), which has HyPhy on its `PATH`. Linux amd64, Singularity CE 4.1.2.

## Inputs

`inst/extdata/phylo_genomes`, made by `tools/dev/validation/make-phylo-genomes.R` (base R, independent of dndsR, seed 20261114; rerunning it reproduced the files byte for byte). 26 files; sha256 of their sorted sha256 list (`find . -type f | LC_ALL=C sort | xargs sha256sum | sha256sum`): `a89de8659846f563c4b67f4cf59ada4f055c3d27c0d48448f91b5b37d3061e36`.

Six genomes, each one chromosome with genes on both strands, a third of them with two exons. Coding sequences were simulated on the species tree `((((spA,spB),spC),(spD,spE)),spF)` under a BUSTED-form model (ω 0.05, 0.5, 1 on 55, 25, 20% of sites; in three orthogroups the third class has ω 8 on the foreground: the clade of spA and spB with its stem). Gene IDs are unrelated to orthogroups. Orthogroups are given as dndsR's long table and as OrthoFinder's `Orthogroups.tsv`, with OrthoFinder-style species names (`spA_proteins`) and gene-tree tips (`spA_proteins_spA0014`).

| Orthogroup | Design | Expected |
|---|---|---|
| OG01-03 | Episodic selection on the foreground | prepared, tested |
| OG04-06 | No foreground difference; gene trees | prepared, tested |
| OG07 | Two spD genes, the second a diverged copy | excluded: `multi_copy` (policy) |
| OG08 | spC gene model with an in-frame TGA | excluded: `cds_invalid` (data quality) |
| OG09 | No spB gene | prepared, tested; spA's branch joined with the stem |
| OG10 | 2-codon insertion in spA; 3-codon deletion on the branch to spD and spE; true alignment stored | prepared, tested |
| OG11 | No spA or spB gene | excluded: `foreground_absent` (no data) |
| OG12 | No spC gene: the stem joins the background branch above | excluded: `foreground_partly_unidentifiable` (policy) |
| OG13 | 1-nt deletion in spE's gene model | excluded: `cds_invalid` (data quality) |
| OG14 | No spB or spC gene: nothing separable | excluded: `foreground_unidentifiable` (mathematical) |

| Gene tree | Discordance | Expected |
|---|---|---|
| OG04 | spD and spF swapped: species split {spD,spE} missing, not beside the foreground | tested, 1 missing split, not near |
| OG05 | {spA,spB,spC} and {spD,spE} missing, the first beside the stem | tested, 2 missing, near |
| OG06 | No spA+spB clade | excluded: `foreground_not_in_gene_tree` (mathematical) |

## Results

### Preparation

`dndsr phylo_prepare` in the v0.4.0 image and the development image: **every orthogroup's status, reason, and basis as designed.** Prepared 8, excluded 6, failed 0, in about 20 s with 4 threads.

| Check | Independent expectation | Result |
|---|---|---|
| Each prepared alignment, gaps removed, equals the designed CDS of its genes, final stop removed | `truth/designed_cds.fasta`, written by the generator | All 8 equal; gene sets equal the orthogroup table |
| OG10, with a codon insertion and deletion | `truth/OG10_true_codon_alignment.fasta` | Equal to the true alignment (after the start-codon column) |
| OG01 alignment, built through `extract_cds`, MAFFT, and back-translation | The same alignment built by `awk` from the designed CDS | Byte-identical |
| OG09 foreground | spA's branch joined with the stem of spA and spB | Only spA labelled; `species_tree_branches` = `mrca(spA,spB)+spA`; `foreground_branches_absent:1` |
| OrthoFinder `Orthogroups.tsv`, and the same as `N0.tsv` | The long table | Identical statuses and byte-identical alignments |
| CDS table instead of genomes | The genome route | Byte-identical alignments; OG08 and OG13 become `gene_not_found`, as no extract_cds QC record says why |
| Repeat run | First run | Byte-identical alignments and trees |
| Member table | Orthogroup table | Every gene listed, both OG07 paralogs with `copies` 2 and unused |

### Results against direct HyPhy

`tools/dev/validation/phylo-direct-hyphy-genomes.sh` runs HyPhy without dndsR: alignments from the designed CDS by `awk`, trees written by hand, JSON read by Python (stored as `truth/direct_hyphy.tsv`).

| Run | Topology | LRT | log L, unconstrained | dndsR against direct |
|---|---|---:|---:|---|
| OG01 | species tree | 6.439360017630861 | -1634.258983526616 | identical |
| OG05 | species tree (also gene-tree `report` mode) | 0 | -1365.805355941506 | identical |
| OG04 | gene tree | 0.5226564205786417 | -1390.842198892926 | identical |
| OG05 | gene tree | 0 | -1408.852804368081 | identical |

**Gene-tree mode, missing from phase 1, reproduces direct HyPhy to every digit.** Discordance counts and "near the foreground" equal the design for OG04 to OG06 in both modes; OG06 is excluded in `topology` mode and tested with flags in `report` mode.

### End to end

From genomes: `phylo_prepare`, then `phylo_busted --prepared-dir`, then `verify_phylo_example`: **11 of 11 checks pass** (the comparison with direct HyPhy, as revised under Processors), with the CLI in the v0.4.0 image plus Bioconda HyPhy, and in the development image alone. One results table lists all 14 orthogroups: 8 tested, 6 excluded at the `prepare` stage, 0 failed. The R workflow in the tutorial runs block by block.

Of the three orthogroups simulated with selection, OG03 had p_adj 0.0027, OG01 0.080, and OG02 0.28; the five without selection had p 0.31 or 0.5. One small example says nothing about power or false positives.

### Processors

The direct-HyPhy scripts (no dndsR code) were run on four node types, HyPhy 2.5.103 Bioconda build, seed 1, 1 thread:

| Processor | Instruction sets | Results |
|---|---|---|
| AMD EPYC 9654 (Genoa) | AVX2, AVX-512 | Reference |
| Intel Xeon Gold 6148 (Skylake) | AVX2, AVX-512 | Bit-identical to the reference, all 8 runs |
| AMD EPYC 7H12 (Rome) | AVX2 | Bit-identical to Milan; differs from the reference in 2 of 8 runs |
| AMD EPYC 7763 (Milan) | AVX2 | As Rome |

The two differences: alignment example OG0001, LRT 7.307789 against 7.307814 (p 0.012945 against 0.012945); genome example OG01, LRT 6.439360 against 6.444410 (p 0.01998 against 0.01993). Both classes are stored in each example's `truth/direct_hyphy.tsv`. Because of this, `verify_phylo_example` and the tests compare dndsR with HyPhy run directly **on the same machine**, on the stored independent inputs in `truth/direct/`, requiring agreement to 1e-9; the stored values are reported, not checked. Found when `check.sh` first ran on an AVX2 node and failed against values stored on an AVX-512 machine.

### The hypothesis is not changed between steps

| Guard | Test | Result |
|---|---|---|
| `phylo_busted --prepared-dir` refuses a species tree or foreground file that differs from what preparation used (MD5 against the preparation's provenance) | Foreground edited after preparation | Stops with an explanation |
| Each orthogroup's labelled tree must equal the one preparation wrote | `trees/OG01.nwk` edited | OG01 `failed: prepared_tree_mismatch`, before HyPhy runs; others unaffected |
| HyPhy's tested branches must equal dndsR's labels | Stand-in HyPhy reporting other branches (phase 1) | `failed: label_mismatch` |

### Failure conditions

All in `tests/testthat/test-phylo-prepare.R` and `test-phylo-tree.R`:

| Condition | Result |
|---|---|
| MAFFT exits nonzero | `failed: alignment_failed` (software), with a warning |
| MAFFT returns proteins altered by one residue | `failed: backtranslation_mismatch`, for every orthogroup it touched |
| A given `--mafft` or `--hyphy` path that does not exist | Stops; does not fall back to another copy |
| Neither or both of genome and CDS tables; a taxon not in the species tree; OrthoFinder columns matching no taxon | Stop with a message naming the problem |
| `--terminal-stops exclude` | Every orthogroup with genes ending in a stop is excluded as `terminal_stop` |
| Hypothesis rows: terminal row lacking a taxon, row with no taxa | `foreground_incomplete` (policy) |
| Stem with a member or sister lineage missing; terminal with sister missing; stem beside the root | `foreground_unidentifiable` (mathematical) |
| Clade with sister lineage missing | `foreground_partly_unidentifiable` (policy) |

### Multiple starts

With `extra_seeds = 2:3` on the phase 1 example's OG0001: three starts ran in separate folders; `lrt`, and so `p_value` and `p_adj`, equal the single-start result and direct HyPhy; `lrt_min_starts ≤ lrt ≤ lrt_max_starts`; the best unconstrained log-likelihood is at least the primary's. The diagnostic arithmetic is tested on constructed cases.

## Tests run

`tools/dev/run-tests.sh phylo` in the v0.4.0 image with Bioconda HyPhy: 5 files, 51 tests, 248 expectations, 0 failed, 0 skipped.

## Limits

- Synthetic data: equal codon frequencies, six taxa, 200 codons, clean gene models apart from designed errors. Real genomes have messier annotations, alignments, and orthology.
- MAFFT's accuracy on difficult alignments is not tested; one designed indel case is. No alignment filtering is applied, and none is validated.
- Agreement with direct HyPhy shows inputs reach HyPhy as intended. It does not show that BUSTED's p-values are calibrated here, or that QC thresholds are appropriate.
- Gene-tree discordance is reported, not modelled. Its effect on false positives is a phase 3 question.
- A development image built locally from `tools/Dockerfile` at `d24d183` (`sha256:7acf9ba81542bb15ae2311de6c8bd8ff20e611c6bc3f022ed8df97b18139f862`, not published) passed all 41 checks of `tools/dev/container-checks.sh`, including both phylogenetic workflows with HyPhy and MAFFT inside it, and the pairwise workflows. CI has not run.
