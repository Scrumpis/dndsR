# Study L: lineage omega, its uncertainty, and its biases (design)

Status: **draft for the maintainer's review; not run.** Motivated by the audit's exploratory findings E1-E3 ([validation-audit.md](validation-audit.md)), which used the first lineage study's outputs. Those genes are the *development set*: anything chosen on them (an interval method, a flag threshold) is confirmed here on new genes with new seeds, and the first study's record is not changed.

## Questions

- **L1, coverage.** Does the 95% interval of `pooled_omega` cover its estimand at 95%, by gene-set size (10, 25, 50, 100, 200 genes), tree, and missingness?
- **L2, bias under site variation.** How large is the pooled omega's bias when omega varies among sites, as a function of branch length and of how much it varies, and does it differ between branches enough to distort comparisons between lineages?
- **L3, heterogeneous gene sets.** With gene-to-gene variation in omega (and in branch effects), does the interval cover the population value it targets?
- **L4, ascertainment.** When a lineage loses or duplicates genes selectively (so they are excluded as `multi_copy` or absent), how biased are unmatched lineage summaries, and do matched-orthogroup comparisons (`_lineage_pairs.tsv`) remove the bias of the *contrast*?
- **L5, consequence for a trait analysis.** If no lineage's omega is related to a trait, but the trait correlates with branch length, does a regression of log pooled omega on the trait across branches give false associations because of length-dependent bias (L2)? Answered from the bias curves measured in L2, propagated through a phylogenetic regression (no extra fits).

## Estimand

For branch *b*, the ratio of the expected number of nonsynonymous substitutions per nonsynonymous site to the expected number of synonymous substitutions per synonymous site, summed over the gene population: with the simulator's scaling (Q normalized at omega = 1, so synonymous rates do not depend on omega), this equals the opportunity-weighted mean of omega over genes and site classes, sum over genes of (p_k omega_k), computed exactly from the simulated parameters. This is what a ratio of pooled counts targets; it is not what a one-omega model fitted to each gene targets when omega varies among sites, and the gap between the two is the bias L2 measures.

## Cells

Trees: T8b, T8p, T16m (existing definitions, `phylo-lineage-sim.R`); 300 codons; one omega per gene unless stated.

| Cell | Site variation | Gene variation | Missingness | Tree | Genes |
|---|---|---|---|---|---|
| L-S0 | none | none | none | T8b | 1,000 |
| L-S1 | classes 0.05/0.3/1 at 0.6/0.3/0.1 (as A3) | none | none | T8b | 1,000 |
| L-S2 | gamma, shape 0.5, mean 0.2 (discretized, 8 classes) | none | none | T8b | 1,000 |
| L-G1 | gamma shape 0.5 | gene omega lognormal, sd 0.7 on the log scale; branch effect lognormal sd 0.5 per tree, fixed | none | T8b | 1,000 |
| L-P1 | gamma shape 0.5 | as L-G1 | none | T8p (short and long branches) | 500 |
| L-M1 | gamma shape 0.5 | as L-G1 | 4 of 16 taxa dropped per gene | T16m | 400 |
| L-A1 | gamma shape 0.5 | as L-G1 | on lineage t1, genes in the top third of omega are made `multi_copy` with probability 0.5 (removed); others unaffected | T8b | 1,000 |

Models fitted: M0 and free-ratio (`phylo_branch_omega`, defaults), then `phylo_summarize` on disjoint gene subsets of each size.

## Criteria (fixed before running)

- **Coverage (L1, L3)**: on L-S0 and L-G1, a 95% interval is *adequate* at a gene-set size if its coverage is at least 93% and the lower 95% bound of the coverage is at least 91%. The smallest adequate size becomes the documented minimum; smaller sets get a flag. Intervals are counted per branch and gene set; precision is reported with a cluster bootstrap over gene sets, because branches estimated from the same genes are dependent.
- **Bias (L2)**: median relative bias within plus or minus 5% on every branch class (short, middle, long), and differences in bias between branches within 5 points, for the pooled omega to be described as comparable across lineages without a branch-length covariate.
- **Ascertainment (L4)**: matched comparisons are *adequate* if the bias of the log ratio of two lineages on shared genes is within plus or minus 0.05.
- **Trait false positives (L5)**: at most 5% (with its interval) at alpha 0.05, for a trait correlated 0.5 with branch length.
- Candidate interval refinements (log-scale interval with t quantiles; a minimum-gene flag) are judged by the same coverage rule; a candidate is adopted only with the maintainer's approval, and the current interval stays reported beside it.

## Precision

L-S0 with 1,000 genes gives 40 disjoint sets of 25 genes times 12 branches (480 intervals; a coverage of 93% is estimated to about plus or minus 2.3 points before the cluster correction) and 10 sets of 100. The bias at a branch over 1,000 genes is estimated to about 1-2% (from E1's spread).

## Seeds and cost

New seeds 201-207, one per cell. Projected from measured times (M0 + free-ratio: about 134 core-seconds per gene on 8 taxa, 372 on T16m): about 6,500 genes on 8 taxa (240 core-hours) and 400 on T16m (41), plus summaries: **about 290 core-hours**. A pilot of 40 genes per cell (about 12 core-hours) checks the simulator's new site and gene distributions and the truth computation before the full runs.

## Not covered

Recombination, alignment error (study B), codon-usage variation between genes, real genomes (study E).
