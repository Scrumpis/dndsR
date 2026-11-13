# Phylogenetic module: calibration design (P3)

Fixed before any simulation was run, and committed before the pilot. Changes after that point are listed at the end with their reason, and do not apply to results already produced. Approved scope: maintainer's P3 instructions, 2026-11-13. Context: [phylo-design.md](phylo-design.md).

## Questions

| | Question | Kind |
|---|---|---|
| Q1 | Under no foreground positive selection, and a correctly specified generating model, is BUSTED as run by dndsR at or below its nominal false-positive rate? | Pass/fail |
| Q2 | How much power is there, by proportion and strength of positively selected sites, gene length, and divergence? | Reported |
| Q3 | Under misspecification (multinucleotide substitutions, strong site variation, alignment error, indels and realignment, a discordant gene tree), how far does the false-positive rate move? | Reported, classified |
| Q4 | Across many genes with some under selection, does Benjamini-Hochberg keep the false discovery rate near its nominal level? | Pass/fail |
| Q5 | How do optimizer strategies affect null false-positive rates, call stability, failures, and runtime? | Reported; informs a recommendation |
| Q6 | How do `--multiple-hits Double` and BUSTED-E change false-positive rates and power? | Reported; informs a recommendation |
| Q7 | Do the QC flags identify the genes that cause false positives? Does alignment filtering, as a sensitivity analysis, help? | Reported |
| Q8 | How much do processor-dependent numerical differences change calls? | Reported |

No default changes on the strength of the pilot. Recommendations go to the maintainer.

## Simulator

`tools/dev/validation/phylo-sim.R`, base R, independent of HyPhy and of dndsR's package code. It shares only the 61-codon table and rate-matrix construction in `codon-model.R`, already used for the pairwise estimator benchmark.

- **Codon model.** A Goldman-Yang-style model: HKY nucleotide mutation (κ) with F3x4 codon frequencies; nonsynonymous rates multiplied by ω; synonymous rates multiplied by a site's synonymous rate class. Branch lengths are expected substitutions per codon at ω = 1 and synonymous rate 1. Transition probabilities by eigendecomposition, cached by (ω, synonymous rate, length).
- **Selection.** Each site draws one of three ω classes, kept across branches. On foreground branches, class 3 may have a different ω (the alternative); on background branches it keeps its null value. BUSTED models branch-site variation as a mixture at each site, so the generating and fitting models differ in detail even when "correctly specified". This affects power, not the null, because under the null no class exceeds 1 anywhere.
- **Synonymous rate variation.** Three site classes with multipliers 0.3, 1, 2.2 and weights 0.3, 0.5, 0.2, rescaled to mean 1. Under strong variation (misspecification M2): a gamma with shape 0.4 discretized into 8 classes, so BUSTED's 3 classes cannot fit it exactly.
- **Multinucleotide substitutions (M1).** Instantaneous changes at two positions of a codon, at rate δ times the product of the two single-change rates, then the same ω scaling.
- **Indels (M4).** Codon insertions and deletions on branches, at rate 0.02 events per unit of branch length each, lengths geometric with mean 2 codons. The true alignment is kept. The CDS then go through `phylo_prepare` (MAFFT and back-translation), so the realignment is the pipeline's own.
- **Alignment error (M3).** In one foreground taxon, a window of 12 codons replaced by random sense codons. This mimics a misannotated exon or a misaligned segment. It is applied to every gene in the cell.
- **Discordance (M5).** Data simulated on a gene tree in which Al and Cr are swapped, ((At,Cr),Al), as incomplete lineage sorting would make, and analysed on the species tree.

### Verification of the simulator, before the pilot

1. Each transition matrix is stochastic, reversible, and stationary at the codon frequencies (to 1e-10).
2. Two-taxon alignments of 20,000 codons, at ω 0.1, 0.5, 1, 2 and κ 2, give PAML codeml (M0, F3x4) estimates of ω and κ within 5% of the truth, and a branch length within 5%.
3. With δ > 0, the fraction of differences at two codon positions over a short branch (0.005 per codon, 100,000 codons) is within its 95% binomial interval of the fraction predicted by the rate matrix.
4. Indel bookkeeping: the true alignment, with gaps removed, equals each simulated CDS.

The pilot runs only if all four hold.

## Parameters from the real-genome benchmark (P3A)

These rules are fixed now, and applied once P3A finishes, before any simulation:

- **Topology:** the species tree `(((At,Al),Cr),(Es,Br))`, foreground the clade of At and Al, as in P3A.
- **Gene lengths:** L25, L50, and L75 are the quartiles of `n_codons` over P3A orthogroups tested with all five taxa, rounded to 10 codons.
- **Codon frequencies:** F3x4, from all prepared P3A alignments pooled.
- **Branch lengths:** the median, across P3A orthogroups tested with all five taxa, of each branch's length in HyPhy's unconstrained model. These are converted to the simulator's units with one factor c. c is chosen on 50 null genes (N0, L50) so that the median HyPhy tree length of the simulated data equals the empirical median within 5%. Only tree lengths are looked at in that step, never p-values. Divergence ×0.5 and ×2 scale all branches.
- **κ:** 2. P3A does not estimate it independently.

## Scenarios

ω classes are (ω1, ω2, ω3) with weights (p1, p2, p3).

| Cell | ω on background, and on foreground for classes 1-2 | Foreground class 3 | Notes |
|---|---|---|---|
| N0: null at the boundary | (0.05, 0.3, 1.0), (0.6, 0.3, 0.1) | 1.0 | Neutral sites on the foreground: the hardest null for BUSTED |
| N1: purifying only | (0.02, 0.15, 0.5), (0.6, 0.3, 0.1) | 0.5 | No neutral or positive sites |
| A(ω+, p+) | N0's classes 1-2, reweighted so class 3 has weight p+ | ω+ | ω+ in {3, 10}; p+ in {0.05, 0.2} (pilot) |

| Cell | Generating change | Null? |
|---|---|---|
| M1 | Multinucleotide substitutions, δ = 0.2 | Yes |
| M2 | Strong synonymous rate variation (gamma shape 0.4, 8 classes) | Yes |
| M3 | 12-codon misaligned window in At | Yes |
| M4 | Indels, realigned by `phylo_prepare` | Yes |
| M5 | Simulated on ((At,Cr),Al), analysed on the species tree | Yes |

Settings are dndsR's defaults unless stated: seed 1, 1 HyPhy thread, `--srv Yes`, `--multiple-hits None`, no error sink.

## Pilot (bounded)

All runs go through `phylo_busted` (and `phylo_prepare` for M4), as users run them, on node-local storage with compact results kept. One processor class (Intel Skylake, AVX-512) unless stated.

| Block | Cells | Genes per cell | Fits |
|---|---|---:|---:|
| Null, correctly specified | N0 at L25, L50, L75; N0 at L50 with divergence ×0.5 and ×2; N1 at L50 | 200 | 1,200 |
| Power | A(3, 0.05), A(3, 0.2), A(10, 0.05), A(10, 0.2), at L50 | 100 | 400 |
| Misspecification | M1, M2, M3, M4, M5 at L50 | 200 | 1,000 |
| Optimizer | On N0-L50 and A(10, 0.05): S1 and S2 (below) | 200, 100 | about 900 |
| Model settings | `--multiple-hits Double` on N0-L50 and M1; `--error-sink Yes` on N0-L50 and M3 | 200 | 800 |
| Processor | N0-L50 rerun on an AVX2 node (AMD Rome) | 200 | 200 |
| Genome-wide | One set of 500 genes: 450 null (N0, lengths L25-L75) and 50 alternative (A(10, 0.1) and A(3, 0.2)) | 500 | 500 |

About 5,000 fits. At 200 null genes, a 95% interval for a false-positive rate near 0.05 is about ±0.03, so the pilot can detect gross inflation, not small departures. The pilot reports and flags; pass/fail is for the full study.

### Optimizer strategies

- **S0, default:** seed 1, 1 HyPhy start.
- **S1:** HyPhy's own `--starting-points 5`.
- **S2:** dndsR's independent starts, seeds 1, 2, 3 (`--extra-seeds 2,3`). The likelihood-ratio statistic is recomputed from the best alternative and best null likelihoods found across starts (`lrt_best_fits`), never by taking the largest LRT.

For each: null false-positive rate; the share of genes whose call at 0.05 differs from S0; the share of starts where the best alternative or null was not found by the first start; negative-LRT refits; runtime.

## Analysis

- **Rates:** false-positive rate and power at α = 0.05 and 0.01, with Clopper-Pearson 95% intervals, per cell. Never averaged across cells.
- **Null distribution:** BUSTED's statistic is compared with HyPhy's reference distribution, a 50:50 mixture of a point mass at 0 and χ²₂. Reported: the share of LRT = 0 (0.5 expected at the boundary, more below it); the empirical distribution of positive LRTs against χ²₂; and p-value histograms with the mass at exactly 0.5 shown separately. Null p-values are not expected to be uniform.
- **Genome-wide:** BH at q = 0.05 and 0.1; observed false discovery proportion and power.
- **QC flags (Q7):** for each flag, in M3 and M4, its rate among false positives against among true negatives.
- **Filtering (Q7):** as a sensitivity analysis in validation scripts only, M3 and M4 rerun after dropping codon columns gapped in more than half the sequences, and windows flagged as dissimilar. Not added to the package.

## Criteria for the full study

Fixed now. Applied in the full study, not the pilot.

| Item | Pass | Concern | Fail |
|---|---|---|---|
| Correctly specified null cell, α = 0.05 (n ≥ 1,000) | 95% upper bound ≤ 0.075 | Point estimate above 0.05, upper bound ≤ 0.075 | Lower bound above 0.05 |
| The same at α = 0.01 | Upper bound ≤ 0.02 | | Lower bound above 0.01 |
| Genome-wide BH at q = 0.05 (pooled over 20 sets) | Pooled false discovery proportion's upper bound ≤ 0.075 | | Lower bound above 0.05 |
| Misspecification cells | Classified, not passed or failed: **robust** (upper bound ≤ 0.075), **inflated** (lower bound above 0.05), **severely inflated** (point estimate above 0.10) | | |
| Optimizer or model setting as a recommended alternative | No correctly specified null cell worse than under the default rule; fewer calls changing across seeds, or better power; runtime at most 3 times | | |

Power is reported with intervals, with no pass/fail. Low power is not an implementation failure. Inflation under misspecification is a property of the method in that condition, not a bug, and is reported per cell, never hidden in an average.

## Full study (proposed, not run)

Sized from the pilot's runtimes and variances, and presented to the maintainer before any launch: n = 1,000 per null cell; the full length and divergence grid; a 10-taxon tree; power over ω+ in {2, 3, 5, 10} and p+ in {0.02, 0.05, 0.1, 0.2}; δ in {0.1, 0.3}; 20 genome-wide sets; and the optimizer and model comparisons across all null cells.

## Changes after this design was fixed

(none yet)
