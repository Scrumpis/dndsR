# Phylogenetic module: calibration pilot (P3)

Simulated data only, plus the parameters taken from public RefSeq genomes (P3A). Design, fixed before any simulation: [phylo-calibration-design.md](../dev/phylo-calibration-design.md). The pilot reports and flags. Pass or fail belongs to the full study, which has not been run. **No default has changed.**

## Question

Does BUSTED, as `phylo_busted` runs it with its defaults, keep false positives near the nominal rate on data like the five-species Brassicaceae set? How much power does it have, and what breaks it?

## Conclusions

1. **Under correctly specified nulls, the false-positive rate is at or below nominal.** At α = 0.05, the rate was 0 to 0.5% in all 9 null cells and in the null genes of the genome-wide set, with 200 genes per cell. Every 95% upper bound is at most 2.8%. The test is conservative here: 68-88% of LRTs are exactly 0 (50% expected at the boundary), and positive LRTs are smaller than χ²₂.
2. **Two kinds of misspecification inflate false positives severely, and the defaults do not protect against them.** Multinucleotide substitutions (M1, δ = 0.2) gave 95%, and a 12-codon misaligned window in one foreground taxon (M3) gave 87.5%. Strong synonymous rate variation (M2) and a discordant gene tree (M5) did not inflate. Indels realigned by `phylo_prepare` (M4) were uncertain (6.5%, interval 3.5-10.9%).
3. **`--multiple-hits Double` removed the M1 inflation (95% to 1.5%), and BUSTED-E (`--error-sink Yes`) reduced M3's (87.5% to 13%, still inflated).** Neither changed the correctly specified null rate (0.5%). Their effect on power was not measured in the pilot.
4. **A single optimizer start often misses the best fit when selection is present.** With seeds 1, 2, and 3 (S2), the best alternative fit improved by more than 0.5 log-likelihood units in 25 of 100 genes, and the best null fit in 16. Taking the best alternative and best null fits (never the largest LRT) changed 9 of 100 calls at α = 0.05: 7 became significant, and 2 stopped being significant. One of those two moved from p 0.048 to 0.060. The other was significant only because its primary null fit had stalled (LRT 20.8, then 0.21). No call changed in the null cell. S2 took 2.8-3.2 times as long.
5. **Power is low for small fractions of selected sites.** ω+ = 10 on 20% of sites: 84%. ω+ = 10 on 5%: 34%. ω+ = 3: 4-12%, which cannot be distinguished from the null rate.
6. **Benjamini-Hochberg made no false discoveries on 500 genes (50 alternatives), and found 5 (q 0.05) or 8 (q 0.1).** One set cannot establish false discovery control.
7. **AVX2 against AVX-512 processors changed no call** in 200 null genes. LRTs were identical in 159 genes, with a largest difference of 0.47.
8. **The QC flag `dissimilar_window` marked 173 of 175 M3 false positives.** It marked 0-0.5% of clean genes, but also 8.5% of M4 genes, where indels, not errors, are the cause. It is a useful flag, not a filter.
9. **Under discordance (M5), 59 of 200 genes were excluded as `foreground_branch_zero_length`.** On these data, HyPhy shrinks the At+Al stem to 0, and dndsR excludes the gene rather than test a narrower foreground. That keeps the test honest, but excludes genes in a way that depends on discordance. The rate is reported per cell, and the policy is an open question.

## Method

- **Simulator:** `tools/dev/validation/phylo-sim.R`, verified before the pilot (design, "Verification").
- **Parameters:** `phylo-calibration-params.R` on P3A, fixed rules. Gene lengths are L25/L50/L75 = 280/420/600 codons. The tree is `(((At,Al),Cr),(Es,Br))`, with foreground the clade At+Al. Scale factor c = 7.360, chosen so that the median HyPhy tree length of simulated null genes (0.406) is within 5% of the real median (0.396), in 2 rounds that looked only at tree lengths. κ = 2, F3x4 from all prepared P3A alignments.
- **Runs:** `phylo-calibration-pilot.sb` runs every cell through `phylo_busted` (and `phylo_prepare` for M4), with dndsR's defaults unless the cell changes one: seed 1, 1 HyPhy thread, `--srv Yes`, `--multiple-hits None`. `phylo-calibration-summary.R` writes rates with Clopper-Pearson intervals, null distributions, paired comparisons, BH, and QC flags. Each run is identified by its directory name, as in the tables.
- **Software:** HyPhy 2.5.103 inside a development image built locally from `tools/Dockerfile` (`sha256:7acf9ba81542bb15ae2311de6c8bd8ff20e611c6bc3f022ed8df97b18139f862`, not published). Intel Xeon Gold 6148 (AVX-512); the processor cell ran on AMD EPYC 7H12 (AVX2).
- **Code provenance gap:** the R library the pilot loaded dndsR from was reinstalled during the run by development checks, so the cells ran on commits between `3a28b1f` and `8af8ca9`. Between those commits, the phylogenetic code changed only for gene trees (`f483869`), which no pilot cell uses. No result is affected, but the full study must use a library that nothing else installs into.
- **Not run:** the filtering sensitivity analysis for Q7 (dropping gappy columns and flagged windows in M3 and M4).

## Results

False-positive rate (null cells) or power, at α = 0.05, with 95% intervals. All genes were tested, with 0 failures, except where noted.

| Cell | Genes | Rate at 0.05 | 95% interval | Rate at 0.01 | Share of LRT = 0 | Classification |
|---|---:|---:|---|---:|---:|---|
| N0, L25 | 200 | 0.005 | 0.000-0.028 | 0 | 0.68 | |
| N0, L50 | 200 | 0.005 | 0.000-0.028 | 0.005 | 0.74 | |
| N0, L75 | 200 | 0 | 0-0.018 | 0 | 0.72 | |
| N0, L50, divergence ×0.5 | 200 | 0 | 0-0.018 | 0 | 0.79 | |
| N0, L50, divergence ×2 | 200 | 0 | 0-0.018 | 0 | 0.70 | |
| N1 (purifying only), L50 | 200 | 0.005 | 0.000-0.028 | 0 | 0.88 | |
| M1, multinucleotide (δ 0.2) | 200 | 0.950 | 0.910-0.976 | 0.875 | 0 | Severely inflated |
| M2, strong synonymous rate variation | 200 | 0.015 | 0.003-0.043 | 0 | 0.71 | Robust |
| M3, misaligned window | 200 | 0.875 | 0.821-0.917 | 0.730 | 0 | Severely inflated |
| M4, indels, realigned | 200 | 0.065 | 0.035-0.109 | 0.035 | 0.50 | Uncertain |
| M5, discordant gene tree | 141 of 200 | 0 | 0-0.026 | 0 | 0.67 | Robust, among tested genes |
| A(ω+ 3, p+ 0.05) | 100 | 0.04 | 0.011-0.099 | 0.02 | 0.50 | |
| A(3, 0.2) | 100 | 0.12 | 0.064-0.200 | 0.03 | 0.28 | |
| A(10, 0.05) | 100 | 0.34 | 0.248-0.442 | 0.18 | 0.21 | |
| A(10, 0.2) | 100 | 0.84 | 0.753-0.906 | 0.74 | 0.05 | |

Settings, on the same genes:

| Comparison | Rate, default | Rate, alternative | Calls changed | LRT identical | Runtime ratio |
|---|---:|---:|---:|---:|---:|
| `--multiple-hits Double`, N0-L50 | 0.005 | 0.005 | 0 of 200 | 140 | 1.5 |
| `--multiple-hits Double`, M1 | 0.950 | 0.015 | 187 of 200 | 0 | 1.1 |
| `--error-sink Yes`, N0-L50 | 0.005 | 0.005 | 0 of 200 | 141 | 1.3 |
| `--error-sink Yes`, M3 | 0.875 | 0.130 | 149 of 200 | 0 | 1.2 |
| S1 (`--starting-points 5`), N0-L50 | 0.005 | 0.005 | 0 of 200 | 139 | 1.1 |
| S1, A(10, 0.05) | 0.34 | 0.34 | 10 of 100 | 15 | 1.2 |
| S2 best fits (seeds 1-3), N0-L50 | 0.005 | 0.005 | 0 of 200 | | 2.8 |
| S2 best fits, A(10, 0.05) | 0.34 | 0.39 | 9 of 100 | | 3.2 |
| AVX2 processor, N0-L50 | 0.005 | 0.005 | 0 of 200 | 159 | |

S2's primary result is identical to S0's by construction (same seed). With `--starting-points 5`, the overall rate was unchanged but 10 individual calls moved in both directions. In S2, the first start was not the best alternative fit in 134 of 200 null genes. The gains were small there (median 0.006, largest 1.85 log-likelihood units).

Genome-wide set (450 null and 50 alternative genes, lengths L25-L75):

| q | Discoveries | False | Power |
|---:|---:|---:|---:|
| 0.05 | 5 | 0 | 0.10 |
| 0.10 | 8 | 0 | 0.16 |

## Real genomes (P3A)

Five public RefSeq genomes (*A. thaliana*, *A. lyrata*, *C. rubella*, *E. salsugineum*, *B. rapa*: GCF_000001735.4, GCF_000004255.2, GCF_000375325.1, GCF_000478725.1, GCF_000309985.2, listed in `tools/dev/validation/phylo-realgenome.sb`, which records their checksums), OrthoFinder orthogroups, then `phylo_prepare` and `phylo_busted` with defaults, foreground At+Al. Rerun after `f483869` with dndsR at `cb87647`, on an AMD node. These are checks of the pipeline, not findings about selection in these species.

| | Orthogroups |
|---|---:|
| OrthoFinder orthogroups | 23,092 |
| Tested (5 taxa / 4 taxa) | 9,864 (9,084 / 780) |
| Excluded: multi-copy / too few taxa (policy, prepare) | 11,490 / 1,563 |
| Excluded: foreground partly unidentifiable (policy, prepare) | 128 |
| Excluded: foreground branch zero length (policy, analysis) | 46 |
| Failed (HyPhy internal error, software) | 1 |

- 48% of LRTs are exactly 0 (p = 0.5). 11.5% have p < 0.05, and 322 (3.3%) have BH q < 0.05.
- **Alignment artifacts are concentrated among the significant genes.** Genes flagged `dissimilar_window` in At or Al are 4.6% of tested genes, but 148 of the 322 with q < 0.05. They are significant at 33% against 1.8% for the others. 162 of the 322 have a foreground ω above 100. In the two strongest flagged genes inspected, a 20-codon stretch of the At sequence shares 11-18% of amino acids with the other four species, while the rest of the gene shares 75-78%. That is the pattern of a misannotated exon, which the pilot's M3 cell showed gives false positives in 87.5% of genes. One unflagged gene with foreground ω of 2,066 has an apparent local frameshift in Al. A real-genome hit list from the defaults needs alignment review before any interpretation.
- **dndsR equals HyPhy run directly on real genes.** `tools/dev/validation/phylo-direct-hyphy-real.sh` writes each tree by hand and runs `hyphy busted` with the recorded settings, on the same image and processor class (AMD EPYC 9654, AVX-512). It was run on 41 orthogroups chosen at random by category with seed 20261114: 20 borderline (0.01 < p < 0.1), 5 significant with flags, 5 significant without, 5 with LRT 0, 5 with four taxa, and the one failure. For the 40 tested, the LRT agrees to 5e-14, the precision of the text table. The failure reproduces directly: HyPhy stops with an internal error in `ComputeBranchCache` on OG0020313, so it is HyPhy's, and dndsR correctly reports it as a software failure.
- In the first run (`3a28b1f`), 918 four-taxon orthogroups were excluded with an R error message as their reason. That was a dndsR bug, fixed in `f483869` and shown here after the fix.

## Limits

- 200 genes per null cell detect gross inflation only (±3% at 5%). One tree, five taxa, κ = 2, one foreground.
- The generating model is close to the model BUSTED fits, so the correctly specified null favors the method. Real genes violate it in ways not simulated here.
- The model settings were compared only on null cells, so their cost in power is unknown.
- The P3A parameters come from the first real-genome run, which wrongly excluded 918 four-taxon orthogroups (fixed in `f483869`). The parameters use five-taxon orthogroups only, so they are unaffected.

## Full study

Proposed to the maintainer with these results. It has not been run, and nothing is launched without approval. Sizes follow the design. The additions, marked +, come from the pilot. Cost uses the pilot's runtimes: about 65 core-seconds per null fit at L50, and 140 per fit under the alternative.

| Block | Cells | Genes per cell | Fits | Core-hours |
|---|---|---:|---:|---:|
| Null, correctly specified | N0 and N1 × L25, L50, L75 × divergence ×0.5, 1, 2 | 1,000 | 18,000 | 350 |
| Power | ω+ 2, 3, 5, 10 × p+ 0.02, 0.05, 0.1, 0.2 at L50; four corners also at L25 and L75 | 300 | 7,200 | 280 |
| Misspecification | M1 at δ 0.1, 0.2, 0.3; M2; M3; M4; M5 | 1,000 | 7,000 | 200 |
| + `--multiple-hits Double` and `--error-sink Yes`, each | 4 null cells, the 16 power cells at L50, M1, M3 | 1,000 / 300 | 21,600 | 850 |
| Optimizer: S2 best fits (seeds 1-3); S1 on a subset | 4 null cells, the 16 power cells at L50 | 1,000 / 300 | 26,400 starts | 900 |
| Genome-wide BH | 20 sets of 500 | 500 | 10,000 | 180 |
| 10-taxon tree | N0-L50; four power cells; M1; M3 | 1,000 / 300 | 4,200 | 230 |
| + Q7 filtering sensitivity, not added to the package | M3, M4, and the P3A genes flagged in At or Al, rerun after dropping flagged windows and gappy columns | | 3,000 | 80 |
| Processor | N0-L50 on AVX2 | 1,000 | 1,000 | 20 |
| **Total** | | | **about 98,000** | **about 3,100 (4,000 with reruns)** |

Conditions: dndsR installed once from one commit into a library nothing else writes to, one image by digest, and the same processor class within each comparison. Criteria as fixed in the design. The S2 runtime ratio in the pilot (3.2 under selection) already exceeds the design's limit of 3, so S2 can be recommended only if a cheaper variant (two seeds) meets it. Recommendations after the study go to the maintainer. They are not applied as defaults.
