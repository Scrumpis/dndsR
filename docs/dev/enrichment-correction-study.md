# Enrichment correction: development and independent validation (design)

Status: **approved by the maintainer 2026-11-23 (`DECISIONS.md`; budget 300-600 core-hours); this design is fixed before any run.** It implements [enrichment-correction-plan.md](enrichment-correction-plan.md). Development (phase D) is exploratory and uses only data already used in development. Confirmation (phase C) uses new public genome pairs and new seeds, never examined before its criteria are committed. No default changes from this study; any proposed default or interface change goes to the maintainer for review before it is merged into stable functionality.

## Questions

- **Q1. Correction.** Which specification of a length- and dS-adjusted enrichment test controls the complete-null probability of any false term and the mixed-null false discovery rate, for GO (hierarchical, overlapping terms), InterPro and custom gene-set families, at `dNdS > 1` and at quantile thresholds, with adequate power?
- **Q2. The adjustment itself.** Does adjusting for dS, which is part of the ratio, bias the test or remove genuine signal? How do model specifications (additive splines, interaction, annotation density, floor) behave under correlated predictors, missing or noisy covariates, extreme dS, and divergence?
- **Q3. `q0.95`.** Why did `covariate_enrichment` fail at `q0.95` (11.0%), and does any specification fix it? If not, the limitation stays documented.
- **Q4. Comparators.** How do GOseq (Wallenius approximation with a probability weighting function) and GOglm (per-term logistic regression with covariates), adapted to a dN/dS-based selection, compare?
- **Q5. Domain diagnostic.** Can a diagnostic computed from the input identify runs where the unadjusted tests are miscalibrated, with measured sensitivity and specificity? An unflagged run is never described as free of confounding.
- **Q6. Default.** Does the correction show calibration superior to the unadjusted tests with adequate power on independent data, so that making it recommended or default is justified?

## Candidates (phase D; all computed on the same labels in each run)

| Name | Model for P(positive) | Term test |
|---|---|---|
| `fisher` | none | hypergeometric (the current `term_enrichment`, `ipr_enrichment`, classic GO) |
| `go_default` | none | the real `go_enrichment()` default (`weight01`), every fifth run |
| `pb_floor` | logistic, `ns(log length, 3) + ns(log dS, 3)`, floor 0.1 × positive rate (current `covariate_enrichment`) | Poisson-binomial upper tail |
| `pb_floor05`, `pb_floor20` | the same with floors 0.05 and 0.2 | Poisson-binomial |
| `pb_gam` | `mgcv::gam` with `te(log length, log dS)` (nonlinear interaction, penalized), same floor | Poisson-binomial |
| `pb_rich` | `pb_floor` plus `ns(log number of terms, 2)` (annotation density) | Poisson-binomial |
| `pb_len` | `ns(log length, 3)` only (dS not adjusted; addresses Q2) | Poisson-binomial |
| `goseq_len`, `goseq_both` | GOseq 1.58.0 probability weighting function on length, and on the fitted `pb_floor` probability as its bias covariate | Wallenius (GOseq's default) |
| `goglm` | per term, logistic regression of the call on the term indicator plus `ns(log length, 3) + ns(log dS, 3)` | one-sided likelihood-ratio test of the term coefficient |
| `cmh_10x10` | none (strata of length and dS deciles) | exact stratified hypergeometric |
| `oracle` | the true P(positive) where the generator knows it | Poisson-binomial (benchmark for power and for Q2; not a candidate) |

GOseq and GOglm are comparators, adapted here: their RNA-seq validation (length bias of differential expression) is not assumed to transfer. GOseq runs from a separate conda environment (`envs/comparators`: bioconductor-goseq 1.58.0, r-biasedurn 2.0.12, R 4.4.3), not a dndsR dependency. GOglm is implemented here from its published description (Mi et al. 2012), as no maintained package exists.

## Generators with known truth

On each pair's real gene table (genes, lengths, observed dS, GO with ancestors, InterPro, and a custom family below). For every scenario the term-level truth is known: in nulls every term is null; in spiked scenarios a term is true if spiked, otherwise null (terms overlapping spiked genes beyond chance are reported separately, as before).

- **Exact labels** (`enrichment-direct-labels.R`, extended): N0 (constant probability); N1 at three strengths (the fitted length-dS dependence on the logit scale × 0.5, 1, 2); N1x (with a length × dS interaction); N1r (dependence on annotation density too); N1c (length and dS made more correlated by a rank-preserving transform, Q2).
- **Term-linked covariates (null for the scientific question):** chosen terms whose genes get a lower true dS (×0.5) and no omega change, so their call rate rises only through the estimator: a fair test must not call them (Q2).
- **Estimator labels** (`sim-genome.R`): `long`, `bylen`, and `bylen` at divergence ×0.5 and ×2 (dS scaled), with dndsR's real Comeron estimator, so labels carry the estimator's own dependence on length and dS, extreme dS included.
- **Planted enrichment:** spikes of 2%, 5% and 20% of genes in 5 terms (20-300 genes) to omega 3 (estimator labels) or to positive (exact labels); and **signal with dS**: spiked terms whose genes also have lower true dS (Q2: does adjustment remove genuine signal?). Power is compared with the oracle's.
- **Missing or noisy covariates (Q2):** the tests given length with 20% multiplicative noise, or dS from a different estimator (NG86 from the same alignments), against labels generated from the true ones.
- **Term families:** GO (all three ontologies, ancestors propagated: hierarchical and overlapping), InterPro, and a **custom family**: 500 random gene sets (sizes 5-500, so sparse terms are included) plus 100 overlapping sets, each the union of 2 or 3 InterPro entries' gene sets (amended before any run: InterPro clan data are not part of the inputs, and unions give the intended overlap and nesting with the InterPro family).
- **Thresholds:** `dNdS > 1`, `q0.9`, `q0.95`.

## Phase D: development (exploratory)

**Data:** the development pairs already used: *A. thaliana* against *A. lyrata*, *Capsella rubella*, *Brassica rapa* (Brassicaceae), and the Solanaceae pair (`SlSt`). The AtAl table is the recorded one (`56321306d1b52a23`). The AtCr, AtBr and SlSt tables are no longer on scratch; they are rebuilt with the current dndsR pipeline from the public inputs whose hashes `enrichment-confounding.md` records (inputs checked against those hashes), and the rebuilt tables are checked against the recorded summaries (number of pairs, median dS, genes above 1); any difference is reported. They remain development data. Seeds 980,000 and up (new, but on development pairs). Runs: 200 per complete-null cell and threshold, 100 per mixed cell, at AtAl and AtBr for every scenario; AtCr and SlSt for the main nulls.

**Selection rule (fixed now, before any phase D run).** Eligible: a candidate whose complete-null share of runs with any false term is at most 5% (point estimate) in every phase D null cell, at `dNdS > 1`; among those, the one with the highest median power relative to the oracle in the planted cells; ties within 2 points go to the simpler model (`pb_floor` simplest, then `pb_len`, `pb_rich`, the floors, `pb_gam`, `goglm`, `goseq_*`). Separately, the same rule applied at `q0.95` tells whether any candidate is eligible there (Q3). The chosen specification is then frozen for phase C.

**Q2 answered in phase D, reported as exploratory:** bias of the adjusted test (complete-null and term-linked-covariate nulls), loss of genuine signal (power against the oracle in the signal-with-dS cells, and against `pb_len`), stability (correlated predictors, noisy covariates, divergence). If adjusting for dS removes genuine signal beyond what the oracle loses, that is reported as a limitation of the method, with the cells where it occurs.

## Implementation (between phases; reviewed before merge into stable)

One shared internal implementation (the frozen specification) used by `covariate_enrichment`, and by `go_enrichment`, `term_enrichment` and `ipr_enrichment` through a new argument `adjust = c("none", "length_ds")`, default `"none"`. With `"length_ds"`, each command's own term family and universe are kept (for `go_enrichment`, adjustment applies to the classic all-ancestor test; topGO's elimination algorithms cannot use a covariate model, and the output says so), and the result gains `expected`, `p_value_adjusted` and `p_adj_adjusted` beside the unchanged unadjusted columns. The unadjusted output stays byte-identical. Tests: an independent reimplementation (Python, `covariate-refit.py` extended to the frozen specification: refits the splines and recomputes the tails) on a public table; hand-computable cases for the Poisson-binomial tail; the unadjusted columns' byte identity.

## Phase C: confirmation (criteria fixed now)

**Data (new, never examined):** three public plant pairs outside Brassicaceae and Solanaceae, chosen before the annotation-coverage gate below:

| Pair | Lineage | Notes | Fallback if the gate fails |
|---|---|---|---|
| *Glycine max* (GCF_000004515.6) - *Phaseolus vulgaris* (GCF_051622355.1) | Fabaceae, **polyploid** (soybean's recent whole-genome duplication) | orthology control below | *Vigna radiata* (GCF_000741045.1) for *Phaseolus* |
| *Sorghum bicolor* (GCF_000003195.3) - *Setaria italica* (GCF_000263155.2) | Poaceae (monocots) | | *Brachypodium distachyon* (GCF_000005505.3) for *Setaria* |
| *Prunus persica* (GCF_000346465.2) - *Malus domestica* (GCF_042453785.1) | Rosaceae, higher divergence | extreme dS | *Medicago truncatula*-*Glycine max* |

**Annotation, uniform across pairs:** GO and InterPro for the query species from its UniProt reference proteome (GO cross-references, including electronic annotations, and InterPro), mapped to RefSeq gene IDs through UniProt's RefSeq cross-references; GO ancestors from GO.db. **Gate (fixed):** at least 60% of the query's genes in the dN/dS table carry a GO term and 60% an InterPro entry after mapping, else the fallback. Mapping rates and annotation density are reported.

**Orthology control for the polyploid:** dndsR's reciprocal best hits as users run them, plus OrthoFinder on both proteomes; the primary analysis uses all RBH pairs, and a sensitivity analysis restricts to pairs whose orthogroup has the expected 1:2 (*Phaseolus*:*Glycine*) structure with the RBH partner in it. Results for both are reported.

**Scenarios:** all phase D generators with new seeds (990,000 and up): 300 runs per complete-null cell (Wilson interval about ±2.5 points at 5%), 150 per mixed and planted cell.

**Criteria for the frozen adjusted test, at `dNdS > 1`, per pair:**

- **C1 complete nulls:** share of runs with any false term at most 5% with the upper Wilson bound at most 7.5%, in every complete-null cell (N0, N1 at three strengths, N1x, N1r, N1c, term-linked covariates, `long`, `bylen` at three divergences), for GO, InterPro and the custom family separately.
- **C2 mixed nulls:** mean false-discovery proportion at most 5% (upper bootstrap bound reported) in every planted cell.
- **C3 power:** reported against the oracle and against Fisher; in N0 (where Fisher is valid), power non-inferior to Fisher's (loss under 5 points, interval lower bound above −10).
- **At `q0.95` and `q0.9`:** the same C1-C2, reported; the adjusted test is described as validated at a quantile threshold only if it meets them on every pair.

**The adjusted test is confirmed** if C1-C3 hold on all three pairs (the polyploid on its primary analysis, the sensitivity analysis reported). **Default rule (Q6):** if confirmed, and the unadjusted test fails C1 in at least one cell of at least one pair, the report recommends making `adjust = "length_ds"` the default of the three commands, for the maintainer's review; it is not made the default in this study.

**Domain diagnostic (Q5), validated on phase C runs:** the likelihood-ratio test of whether the call depends on `ns(log length, 3) + ns(log dS, 3)`, at α = 0.05, computed on every run. *Positive condition:* the unadjusted test's family-wise rate in that cell exceeds 5% (lower Wilson bound above 5%). Reported: sensitivity (share of runs in positive cells that are flagged) and specificity (share of runs in negative cells not flagged), and the unadjusted test's false-term rate among unflagged runs. **Criterion:** among unflagged runs in every cell, the unadjusted test's rate at most 5% (upper bound at most 7.5%). If it fails, the diagnostic is reported as unreliable, and the unadjusted tests carry a warning, not a pass/fail flag. Whatever the result, the output never says an unflagged run is free of confounding.

**Independent implementation checks:** the frozen specification against the Python reimplementation on one full phase C table (all terms, both thresholds; agreement within 1e-8 in fitted probabilities and 1e-10 relative in tails), and GOseq's own implementation as an external comparator.

## Empirical examples

On the phase C pairs' real labels, the adjusted and unadjusted term lists are reported side by side as examples. They are not calibration evidence.

## Cost

Exact-label runs: about 80 core-seconds each for the current candidates, plus GOglm (one fit per term) and GOseq; estimator-label runs add the codon simulation and estimator. Measured on a pilot of 20 runs per generator on AtAl before the main runs, and recorded here. Projection before measuring: phase D about 150-250 core-hours, phase C about 250-350, within the approved 300-600. Intel nodes only for consistency.

## Records

Every run records the seed, generator, pair, commit, image digest, comparator environment, and the full term table; summaries per cell with intervals go to `docs/validation/data/`. Development and confirmatory results are reported separately: `docs/validation/enrichment-correction-development.md` and `enrichment-correction-confirmation.md`.

## Setup record (before any phase D or phase C result)

- Comparators: `envs/comparators` (micromamba; bioconductor-goseq 1.58.0, r-biasedurn 2.0.12, R 4.4.3) plus geneLenDataBase 1.42.0 installed from Bioconductor 3.20's tarball (SHA-256 `f731cfa1ff83a38d`), because the conda data package left it uninstalled. GOseq runs with dndsR's library paths cleared.
- GOglm is computed as the logistic score test of the term indicator given the covariates (null model fitted once per universe; asymptotically equivalent to the per-term likelihood-ratio test, and tractable for thousands of terms).
- Inputs downloaded by `enrich-corr-download.sb` (`enrich-corr-genomes.tsv`), hashes in the run folder's `sha256.tsv`. Development genomes match the hashes recorded in `enrichment-confounding.md` (*C. rubella*, *B. rapa*, tomato, potato; the potato GFF after removing its plastid sequence, NC_008096.2, as the recorded "nuclear GFF").
- `sim-genome.R` gained `tlds`, `spikeds_F_W` and a divergence factor; existing scenarios regenerate byte for byte (`bylen` on 300 genes and `spike_0.2_3` on 3,000 genes, seed 4242, old against new).
- Exact-label quantile counterparts: with binary labels a quantile threshold cannot be applied, so the exact-label generator sets the positive share directly (observed, 5%, 10%) as the counterpart of `dNdS > 1`, `q0.95` and `q0.9`; estimator-label runs apply the real thresholds.
- Cost pilot (development seeds 995,000+, AtAl, 17 blocks): about 120 core-seconds per exact-label run (all candidates including GOseq and GOglm) and 490 per estimator-label run (simulation plus three thresholds). Phase D on AtAl as launched: 4,250 runs, about 208 core-hours; with the three other development pairs, phase D projects to about 290 core-hours, above the design's 150-250 projection; with phase C, the study projects to the top of the approved 300-600. Recorded before any result.

