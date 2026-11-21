# Scientific validation audit

Started 2026-11-21, on `feat/lineage-selection`, under the maintainer's instruction of that date (`DECISIONS.md`): validate until the methods earn high confidence, and release v0.5.0 only when each method meets its acceptance criteria or its remaining limits are accepted by the maintainer. This page is the working record: the evidence matrix, the gaps that matter, and the studies run or proposed to close them. It summarizes existing evidence and links to it; it does not alter any earlier study's record.

McDonald-Kreitman is owned by the MK maintainer and stays experimental by intent; it is listed for completeness and not assessed here.

## How to read the matrix

Each method is judged on seven dimensions, kept separate:

- **Correctness**: dndsR computes what it says (exact or near-exact agreement with an independent calculation or the engine run directly).
- **Calibration**: a test's false-positive rate, or an interval's coverage, matches its nominal level under its own null, in simulation.
- **Power / accuracy**: detection rate, or estimation error, against known truth.
- **Robustness**: the same, when real-data problems are present (alignment and orthology errors, missing taxa, discordance, model misspecification).
- **Reproducibility**: same inputs and version give the same answer (threads, seeds, nodes).
- **Supported use** and **known limits**.

Status words: **Pass** (meets a stated criterion, with the evidence linked), **Fail** (misses it), **Insufficient** (not measured, or measured too narrowly to support the intended use), **n/a**. A pass is only as wide as the conditions tested.

## Evidence matrix

### Pairwise pipeline

| Method | Correctness | Calibration | Power / accuracy | Robustness | Reproducibility | Status for stable |
|---|---|---|---|---|---|---|
| CDS extraction (`extract_cds`) | **Pass**: gffread on a public genome, 0 differences; designed example exact ([cds-extraction](../validation/cds-extraction.md), [public-genome](../validation/public-genome.md)) | n/a | n/a | Partial: first exons with nonzero phase not trimmed (rare, [known issues](../validation/known-issues.md)); trans-spliced organelle genes stop with a message | **Pass**: byte-identical | Ready, with the known issue documented |
| Orthology (reciprocal best hits) | **Pass** as implemented | n/a | Against OMA: 3.6% disagreement (*A. lyrata*), 15% (*B. rapa*) ([orthology](../validation/orthology.md)) | **Insufficient**: disagreement rises with dS and polyploidy; OMA is not truth; 35-37% of pairs unscored | Pass | Ready as "one-to-one RBH", with the stated limits; not for polyploids without subgenome splitting |
| dN/dS estimation (Comeron default; NG, YN, others) | **Pass**: PAML, KaKs_Calculator, independent NG86 ([estimators](../validation/estimators.md)) | n/a (no test) | Comeron within about 15% of omega for omega 0.1-2 at dS up to about 0.3; compresses ratios at high divergence | Simulation model has one omega per gene, no recombination | **Pass** across threads | Ready within the stated dS range; no per-gene uncertainty (a bootstrap was tested and rejected: it misses bias) |
| Reliability flags (high dS, saturation) | Pass (unit tests) | n/a | dS above 3x median flags paralogs (85% disagreement with OMA) | One plant family | Pass | Ready as flags |

### Enrichment

| Method | Correctness | Calibration | Power | Robustness | Status for stable |
|---|---|---|---|---|---|
| `go_enrichment`, `term_enrichment`, `ipr_enrichment` (Fisher/hypergeometric) | **Pass**: exact hypergeometric, BH | **Fail**: false terms in up to 43-70% of runs under estimator-noise nulls linked to length and dS ([enrichment-confounding](../validation/enrichment-confounding.md)); warn on every run | Comparable to `covariate_enrichment` | Pair-dependent excess | **Fails calibration within its intended use.** Not to be described as validated; the maintainer's open question on defaults stands (DECISIONS, open questions) |
| `covariate_enrichment` (`pb_floor`) | **Pass**: independent recomputation | **Pass** on 4 plant pairs, 2 families, selection + confirmation + held-out pair (worst cell 5.0%, upper bound 11.2% at 100 runs) | 61-100% of spiked terms | **Insufficient** beyond: polyploids, other clades, quantile thresholds, annotation completeness | Supportable as validated for the tested scope; the 5.0% (2.2-11.2) worst cell is at the edge, see gaps |
| `dnds_summary`, contrasts, ideograms | Pass (base-R recomputation) | `dnds_state_contrast` zero cells and suppressed warnings open; gene level disabled | n/a | n/a | Descriptive tools ready; `dnds_state_contrast` exploratory |

### Phylogenetic module

| Method | Correctness | Calibration | Power / accuracy | Robustness | Reproducibility | Status for stable |
|---|---|---|---|---|---|---|
| `phylo_prepare` | **Pass** on designed genomes: alignments equal designed CDS, true indel alignment, exclusions as designed ([phylo-prepare](../validation/phylo-prepare.md)) | n/a | n/a | Realigned indels cause BUSTED's main residual false positives (8.3%); real OrthoFinder output **not tested** | Pass | Insufficient: real OrthoFinder input and alignment-error impact on lineage tests |
| BUSTED (`phylo_busted`, C1: double hits, 3 starts) | **Pass**: equals HyPhy directly ([wrapper](../validation/phylo-busted-wrapper.md)) | **Pass** under correct nulls, 10 cells (0-1.3%) ([2B](../validation/phylo-calibration-stage2b.md)); **Fail** under realigned indels (8.3%) and unmasked foreground corruption (FDP 77% after BH in a genome-like mixture) | Low for small selected fractions; little survives BH on 5 taxa | One tree (5 taxa), one foreground; masking is a sensitivity analysis only (decided) | AVX2/AVX-512 changed no call | **Insufficient** for general use: one tree; robustness failures documented but not resolved for other trees and real alignment error |
| PAML branch models (`phylo_branch_omega`) | **Pass**: equal to codeml directly within 1e-3 lnL ([phylo-lineage](../validation/phylo-lineage.md)) | `p_two_vs_m0`: consistent with nominal (3.3-7.0%), never shown calibrated. `p_fg_omega_gt1`: conservative below 1; **inconclusive at omega 1** (5.5%, 6.8%). `p_free_vs_m0`: inconclusive with missing taxa (8%, then 4%). BH-corrected: no false positive ([calibration](../validation/phylo-lineage-calibration.md)) | Two-ratio difference 98-100%; omega above 1 at 1.5: 78%; per-branch free-ratio omega 1.3-1.4-fold median error | **Insufficient**: divergence, length, long branches, starting values, alignment error, misassigned orthologs not tested | Pass (deterministic) | **Insufficient** |
| aBSREL (`phylo_absrel`) | **Pass**: identical to HyPhy | `p_raw`, `p_bh_lineage` **Pass** on 4 trees; gene `p_value` Pass on T5, consistent elsewhere | `p_bh_lineage` 93-96%, `p_holm` 72-79% (omega 10 on 10% of sites) | **Insufficient and high-risk**: alignment errors, which ruin BUSTED unmasked, are untested for aBSREL; FDR in a mixture untested | Pass | **Insufficient** |
| RELAX (`phylo_relax`) | **Pass** on the example; directly through HyPhy, 286/300 genes identical, 14 differ by optimizer non-determinism | **Fail**: 8.3-9.3% at 0.05 on 16 taxa, reproduced; 7.0% (n.s.) on 5 taxa; 5.0% on 8 | k = 0.5: about 60%; k = 2: about 37% | Not tested beyond the simulator | Optimum not reproducible to 1e-3 (5% of genes) | **Fails calibration** |
| Lineage summaries (`phylo_summarize`) | Pass (hand-built tables) | **Insufficient; exploratory evidence of a problem** (below): percentile interval covers 91% at 25 genes when the model is right; under site variation the pooled omega is 4-11% low and coverage falls to 70% at 100 genes | Pooled point estimate within 3-6% of truth with one omega per site | Gene heterogeneity, matched vs unmatched sets, single-copy ascertainment: **not tested** | Pass (seeded) | **Insufficient** |
| McDonald-Kreitman | Owned by the MK maintainer | | | | | Experimental by intent |

## Exploratory findings from existing outputs (2026-11-21)

Reanalyses of the first lineage study's stored outputs (`lineage-calib2`, at `dcf3231`) and the follow-up's (`lineage-follow1`). They were not prespecified, are labelled exploratory, and are not folded into those studies' results; each motivates a prespecified study below. Scripts: `tools/dev/validation/phylo-lineage-pooled-coverage.R` (E1), `phylo-lineage-relax-diag.R` (E2), `phylo-lineage-absrel-pooled.R` (E3).

**E1. Pooled omega intervals.** Disjoint gene sets of 25, 50, and 100 from the A0-A3 cells, interval as `phylo_summarize` computes it (percentile bootstrap of orthogroups, 1,000 resamples, high-dS genes excluded). Clopper-Pearson intervals below ignore the dependence among branches estimated from the same genes, so they are too narrow.

| Truth | Gene set | Intervals | Coverage | Median log2(estimate/truth) |
|---|---|---:|---|---:|
| One omega per site (A0, A1, A2; 5 trees) | 25 | 1,424 | 91.2% | about 0.01 |
| | 100 | 356 | 94.1% | about 0.01 |
| Sites vary, classes 0.05/0.3/1 (A3; 2 trees) | 25 | 216 | 84.7% (79.2-89.2) | -0.069 |
| | 100 | 54 | 70.4% (56.4-82.0) | -0.067 |

When sites vary, the truth lies above the interval in every miss (33 of 33 at 25 genes): the free-ratio model's single omega per branch underestimates the mean omega (4% low on branches of 0.1-0.4 substitutions per codon, 11% on the longest), and more genes only narrow the interval around a biased value. Real genes always vary among sites. With one omega per site, coverage is slightly under nominal at small gene sets (91%) and near it at 100. On T16m (missing taxa) 90%.

**E2. RELAX's statistic.** Over the 1,500 null genes of the three independent sets, the LRT's quantiles are 1.18-1.46 times chi-square(1)'s (median 0.66 against 0.45; 95th percentile 4.86 against 3.84), close to a scaled chi-square with factor 1.27; it is not a few outlying genes. Genes whose alternative fit has a small class (proportion below 5%) at omega above 2, absent from the simulation, have a larger excess (11.4% of 307) than the rest (7.8% of 1,193). Weak identification of the omega classes contributes, but does not explain the whole excess. Significant genes split between k below and above 1.

**E3. A candidate that failed: pooling aBSREL's branch dN and dS.** On the B0 cells (the same site classes as A3, true mean omega 0.22 on every branch), the ratio of summed `adaptive_dn` to summed `adaptive_ds` was 0.85-2.6 by tree, 3-12 times the truth: those columns are not on a per-site dN/dS scale, and HyPhy does not report the site counts needed to put them there. Not pursued. Incidentally, aBSREL's single-omega baseline per gene had a median of 0.27 against 0.22, biased high, while codeml's pooled free-ratio omega is biased low (E1): under site variation, what a "branch omega" estimates depends on the model, so study L defines its estimand explicitly.

## Gaps that matter, ranked

Ranked by how much each could change a scientific conclusion in the central use (branch-level estimates and tests related to LTR activity and genome size), then by how cheaply it can be resolved.

1. **Lineage omega uncertainty and bias (E1).** The pooled omega is the quantity a comparative model of genome size would use. Its interval undercovers, and it is biased under site variation, which is universal in real genes. Study L below.
2. **Robustness of the lineage tests to alignment and orthology error.** BUSTED's dominant real-data failure (corrupted foreground segments: 77% false discoveries after BH) was never measured for aBSREL or the branch models, and the first lineage studies simulated perfect alignments. Study B below.
3. **RELAX calibration (Fail).** Study R below.
4. **Branch models across divergence, length, and tree shape**, including the omega = 1 boundary at more conditions and the optimizer's starting values and bounds. Study P below.
5. **Single-copy ascertainment.** Lineages with more duplication lose more orthogroups; pooled omega is then computed on a different gene set per lineage. Matched-orthogroup comparisons (`_lineage_pairs.tsv`) address it in part; never measured. Part of study L.
6. **A realistic multispecies public dataset end to end**, including real OrthoFinder output. Study E below. Agreement with biological expectation is not evidence of correctness; this checks that the workflow runs and that diagnostics behave on real data.
7. **Covariate enrichment at its edge** (worst cell 5.0%, 2.2-11.2): more runs in that cell would show whether it holds. Lower priority for the lineage use case.
8. **BUSTED beyond five taxa**: the lineage use case runs on larger trees. Folded into study B.

Fails that are documented and not proposed for more work now: the Fisher-type enrichment commands (calibration fails by design; the remedy exists as `covariate_enrichment`; whether they stay is the maintainer's open question).

## Studies

Each study is designed in full (question, hypotheses, nulls and alternatives, criteria with precision, seeds, cost, pilot) in its own design file before it runs, and reported in its own results file. Budgets below are projections from measured run times (`lineage-calib2/summary/timing.tsv`: per gene and model set, codeml M0 + free-ratio on 16 taxa about 990 core-seconds, on 8 taxa about 134; two-ratio + null on 8 taxa about 36; RELAX on 16 taxa about 44; aBSREL on 16 taxa about 64).

| Study | Question | Status | Projected cost |
|---|---|---|---|
| R | Does a parametric bootstrap of RELAX's LRT, simulated from HyPhy's own fitted null, control size, and at what power and cost? | Pilot done (bootstrap reproduces the inflation); main study proposed: [relax-bootstrap-study.md](relax-bootstrap-study.md) | Pilot 4 core-hours; main about 930 |
| L | Coverage and bias of pooled omega under site variation, gene heterogeneity, missingness, tree shape; interval methods; matched vs unmatched; ascertainment; consequence for a trait regression | Proposed: [lineage-summary-study.md](lineage-summary-study.md) | About 290 |
| B | aBSREL and branch-model tests under alignment corruption, realigned indels, a misassigned paralog, missing taxa, and discordance, including false discoveries after correction in a genome-like mixture, with and without masking; BUSTED on the same larger trees | Proposed: [lineage-robustness-study.md](lineage-robustness-study.md) | About 250 |
| P | Branch models across divergence (x0.25, 1, 4), length (100, 300, 1,000 codons), foreground omega (0.5, 1, 2), trees; starting values and bounds | Design after B | About 265 |
| E | Public multispecies dataset end to end (genomes varying in size and TE content), with real OrthoFinder | Design last | About 600 |

## Checkpoints

### 2026-11-21, audit start

- **Tested:** nothing new was simulated; existing outputs were reanalysed (E1, E2), and the RELAX bootstrap's mechanics were checked on one development gene: RELAX's `--save-fit` leaves its result unchanged (LRT 0.2843813938652602, equal to dndsR's), the null can be rebuilt from the saved fit by constraining `relax.K = 1` (log-likelihood within 0.0007 of RELAX's own null fit), and HyPhy simulates from it.
- **Learned:** E1 and E2 above.
- **Unresolved:** everything marked Insufficient or Fail.
- **Stable release:** not supported by the evidence for any phylogenetic method; the pairwise pipeline and `covariate_enrichment` are supportable within their tested scope.
- **Next:** the RELAX pilot's result, then study L's design.

### 2026-11-21, after the RELAX pilot and interval checks

- **Tested:** the RELAX bootstrap pilot (16 development genes, 20 replicates each, 4 core-hours); interval methods for pooled omega (percentile, BCa, log-normal) on the development outputs; aBSREL-based pooled omega (E3).
- **Learned:** HyPhy's own simulation from the fitted RELAX null reproduces the inflated LRT (mean 1.25, SE 0.10; 7.2% above 3.84), so the excess is RELAX's, and a bootstrap can in principle correct it; the refitted null matches RELAX's to 0.026 log-likelihood units; about 60 core-seconds per replicate. The choice of bootstrap interval does not fix pooled omega's undercoverage at small gene sets (89% / 92% / 94% at 10 / 25 / 50 genes for all three methods). aBSREL's branch dN and dS cannot be pooled as omega (E3).
- **Changed:** the methods page and skill now state the pooled-omega interval's undercoverage and the bias under site variation, labelled exploratory.
- **Unresolved:** all Fail and Insufficient rows; RELAX's correction is not validated (a pilot on development genes is not evidence of calibration).
- **Stable release:** not supported for any phylogenetic method.
- **Next, for approval:** study R (about 930 core-hours), study L (about 290), study B (about 250); P and E to be designed after them.

