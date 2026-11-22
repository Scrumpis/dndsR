# Study L: lineage omega, its uncertainty, and its biases (design)

Status: **approved by the maintainer 2026-11-21 (pilot first, then the full study if the pilot passes its gates, within about 290 core-hours); amended below before the pilot ran.** Motivated by the audit's exploratory findings E1-E3 ([validation-audit.md](validation-audit.md)), which used the first lineage study's outputs. Those genes are the *development set*: anything chosen on them (an interval method, a flag threshold) is confirmed here on new genes with new seeds, and the first study's record is not changed.

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

## Amendments before the pilot (2026-11-21, before any run of this study)

Made after the maintainer's approval, which asked that estimator validity come first: compare alternative, biologically interpretable estimators; establish whether any genome-wide lineage summary can be used in comparative trait models; and, if none can, document that rather than force a biased estimator into the stable API. Nothing below was chosen from results of this study.

**Implementation.** `phylo-lineage-sim2.R` (cells `LS0`, `LS1`, `LS2`, `LG1`, `LA1`; it reuses `phylo-lineage-sim.R` without changing it), `phylo-lineage-study2.sb` (one unit per cell), units `phylo-study-L-units.tsv` (seeds 201-207 as designed) and `phylo-study-L-pilot-units.tsv` (40 genes per cell, seeds 291-297: development data, never pooled with the main runs). L-P1 and L-M1 use cell `LG1` on T8p and T16m.

**Gene and branch effects.** omega of gene g, branch b, site class k is c_k x m_g x e_b: c_k the site classes (gamma: 8 equal classes, conditional means, mean 0.2), m_g lognormal with sd 0.7 on the log scale and mean 1, e_b lognormal with sd 0.5, one fixed draw per tree (seed 20261122). Population truth for a branch: 0.2 x e_b (m has mean 1). Each gene's own truth is written too (`gene_branch_truth.tsv`), so the estimand of a particular gene set (sample truth) and the population value are both available; L3 reports coverage of both.

**L-A1, the ascertainment mechanism.** A duplicated gene makes its orthogroup `multi_copy`, which removes it from every lineage: that changes the gene population of all lineages alike and so cannot bias a contrast between lineages. The mechanism that can is lineage-specific loss: genes whose m_g lies in the top third of its distribution lose t1 with probability 0.5 (an absent taxon, `n_absent`), so t1's unmatched summary is computed on fewer fast-evolving genes. L4 compares t1 against t3, unmatched and matched (not against t2: when t1 is absent, t2 merges with their stem, so t2 is estimated on t1's gene set too).

**Estimators compared** (all from the same M0 and free-ratio fits; criteria L1-L5 apply to each):

| Name | Definition | Interpretation |
|---|---|---|
| `pool` | current `pooled_omega`: (sum NdN / sum N) / (sum SdS / sum S), high_ds excluded | ratio of pooled substitution counts per site |
| `median` | median over genes of the branch's free-ratio omega | typical gene; its estimand is the median gene's omega (0.2 x e_b x exp(-0.245) in the cells with gene effects), not the mean |
| `rel` | `pool` of the branch divided by `pool` of M0 over the same genes | lineage omega relative to the tree-wide omega of its own gene set (cancels gene-set composition, which L4 tests) |
| `pool_adj` | `pool`, with log branch dS as a covariate in the trait model (L5 only) | remedy by adjustment rather than by a new estimator |

aBSREL's per-branch omega distribution was considered and not included: on the first study's B0 genes most branches fit a single class (for example `Rate Distributions [[0.21, 1]]`), so it has the same one-omega misspecification, and its baseline was already biased high (audit, E3).

**L5, made concrete.** On L-P1 (T8p, branch lengths 0.005-0.25) and L-G1, with disjoint gene sets of 100: 10,000 trait draws per gene set, each trait correlated 0.5 with the part of log branch length that is uncorrelated with the true log omega of the branches (so that the trait is exactly unrelated to the truth, even though one fixed draw of branch effects can correlate with length by chance), plus normal noise; ordinary least squares across branches of log estimate on the trait, with and without log branch dS as a covariate; false-positive rate at 0.05. The trait's noise varies independently among branches, so the test has its nominal size when the estimator has no length-dependent bias; the same regression on the true omegas is reported as the reference.

**Correction before the main runs (2026-11-21, from the pilot's development data).** The summary first compared `median` with the population mean, which builds in a bias wherever gene effects vary; `median` is now compared with its own estimand, as in the table. No threshold changed.

**Decision rule for the stable API (fixed now).** An estimator may be documented for comparative trait models only if it meets L2 (bias) and L5 (false positives at most 5%) on every cell with site or gene variation. If none does, but adjustment (`pool_adj`) meets L5, the documented use is the adjusted model, and the limitation is stated. If neither, the documentation says that no genome-wide single-number lineage omega from these models is safe for trait association, and the summary stays descriptive.

**Pilot gates (fixed now; the full study starts only if all hold).** (1) every unit completes, with at least 95% of genes tested; (2) the simulator does what is planned: mean of m_g within 3 standard errors of 1, the gamma classes average 0.2, and in L-A1 t1 is absent only in genes above the cutoff (a logical check; the share, planned 16.7%, is reported, since 40 genes cannot pin it down); (3) on LS0 the pooled omega lies within 10% of 0.2 on every branch class (the first study found 3-6%); (4) measured cost per gene projects the full study at no more than 290 core-hours; if it projects above, the cells are run in the order L-G1, L-P1, L-S1, L-S2, L-A1, L-S0, L-M1 and the remainder reported as not run, without reducing any cell's genes.

## Pilot outcome (2026-11-21, development data; gates as fixed above)

`phylo-study-pilot-gates.R L` on the 7 pilot units (40 genes each, seeds 291-297), 760ca03:

- **Passed:** every unit complete, 100% of genes tested; mean gene effect 1.045 (standard error 0.064, 160 genes); gamma classes average 0.2 exactly; in L-A1, t1 lost only above the cutoff (7 of 40 genes, 17.5%, planned 16.7%); projected cost 196 core-hours (budget 290; measured 101 core-seconds per gene on T8b, 142 on T8p, 323 on T16m).
- **Failed: gate 3**, pooled omega within 10% of 0.2 on every branch of LS0: 12 branches from 0.179 to 0.211, so t7 (0.179, -10.5%) is outside. Every branch's 95% bootstrap interval contains 0.2 (the intervals are 12-15% wide either side at 40 genes), and the mean over branches is 0.195 (-2.4%), in line with the first study's 3-6% at 300 genes. The gate, applied to point estimates of 12 branches at 40 genes, could fail by sampling alone. It failed as written, so **the full study is not started**; the maintainer decides whether to proceed. No threshold was changed.
