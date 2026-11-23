# What dndsR results establish: a method-by-method limitations and diagnostics matrix (preliminary)

Status: **preliminary assessment for the maintainer's review (2026-11-23); no code, default or output changes.** Requested as deliverable 2 of the maintainer's message of 2026-11-23 (interpretation). It consolidates dndsR's own validation and the primary literature. Every entry carries its basis:

- **[V]** demonstrated in dndsR's own validation (report linked);
- **[L]** established in the published literature (reference given; *verified* means the source or its abstract was read in this session, *secondary* means seen only through a citing source);
- **[U]** a plausible concern not yet tested by us, never to be presented as a validated QC criterion;
- **[D]** a diagnostic dndsR already computes reproducibly from its outputs (name of the column or file).

It is not a QC rule set. A diagnostic listed here does not mean a result is wrong when it fires, or right when it does not.

## Release-ready or candidate methods

### Pairwise dN/dS (`calculate_dnds`; Comeron default; NG, YN, LWL, LPB through KaKs_Calculator)

| | |
|---|---|
| Question answered | The ratio of nonsynonymous to synonymous divergence between two sequences, averaged over the whole gene and the whole time since divergence. |
| Assumptions | Correct one-to-one orthology; correct CDS and alignment; dS not saturated; one ω summarising all sites and the whole lineage pair. |
| Known failure modes | Underestimates ω at high divergence and compresses ratios [V: [estimators](../validation/estimators.md)]; orthology errors rise with dS and polyploidy (3.6-15% disagreement with OMA, 35-37% of pairs unscored) [V: [orthology](../validation/orthology.md)]; first exons with nonzero phase not trimmed [V: known issues]. A ratio above 1 for a whole gene is rare even under selection, because purifying selection on most sites dilutes it [L, textbook; not re-checked here]. Between closely related lineages or within species, ω also reflects segregating polymorphism, not only fixed differences (Kryazhimskiy and Plotkin 2008, *PLoS Genet*) [L, not re-checked in this session]. |
| Diagnostics [D] | `dS`, high-dS and saturation flags; reliability flags; RBH status. |
| Warn or withhold | dS outside the validated range (about 0.3 for Comeron's accuracy) [V]; polyploid comparisons without subgenome-aware orthology [V]. |
| Not justified | "Gene X is under positive selection" from dN/dS > 1 alone; "no selection" from dN/dS < 1; comparisons of ω between pairs at very different divergence. |

### Enrichment (`go_enrichment`, `term_enrichment`, `ipr_enrichment`, `covariate_enrichment`)

| | |
|---|---|
| Question answered | Fisher-type: whether the share of genes above the dN/dS threshold differs between a term and its universe. `covariate_enrichment`: whether it differs beyond what gene length and dS predict. |
| Assumptions | Genes independent; annotation unbiased with respect to the call; for selection claims, that length and dS are the main confounders. |
| Known failure modes | Confounding by length and dS gives 33-70% complete-null false-term rates for the Fisher-type tests at one pair [V: [study N](../validation/enrichment-error-calibration.md), [confounding](../validation/enrichment-confounding.md)]; `covariate_enrichment` validated at `dNdS > 1` on 4 plant pairs (worst cell 1.2%, 0.7-2.1) but not at `q0.95` (11%) [V]. Length bias in categorical enrichment is established for RNA-seq (Young et al. 2010, *Genome Biol*; Mi et al. 2012, GOglm) [L, verified]. Annotation completeness, gene dependence (tandem families) [U]. |
| Diagnostics [D] | `pos_total`, `n_tested`, `expected`, `obs_over_expected` (`covariate_enrichment`); the calibration notice file. |
| Warn or withhold | Any Fisher-type term offered as evidence of selection [V]; few positive genes (low power) [V]; quantile thresholds with `covariate_enrichment` [V]. |
| Not justified | "Function F evolved under positive selection" from an enriched term; "no functional enrichment" from no significant term. |

### BUSTED (`phylo_busted`, recommended C1: double hits, three starts)

| | |
|---|---|
| Question answered | Whether at least one site on at least one foreground branch evolved with ω > 1 (episodic diversifying selection somewhere in the foreground). Not which sites or branches; not that the gene as a whole is selected (Murrell et al. 2015, *MBE*) [L, cited in dndsR docs]. |
| Assumptions | Correct alignment and homology, correct tree and foreground, codon model adequacy (with SRV and multiple hits as configured). |
| Known failure modes | One corrupted foreground segment: 89-96% false positives on 5, 8 and 16 taxa [V: stage 2B, study B]; realigned indels 6-8% [V]; multinucleotide changes without the double-hit option 97% [V; Lucaci et al. 2023, *MBE*, L]; selection and alignment error both look like a small class of high ω, which is why errors are hard to separate [L: Minus the Error, eLife reviewed preprint, verified]. |
| Diagnostics [D] | `dissimilar_windows`, `qc_flags`, the `Alignment check` warning and its provenance record; `n_starts`, `lrt_best_fits` and the agreement of starts; fitted rate classes (`fg_omega_*`, `fg_prop_*`). |
| Warn or withhold | A significant gene with a flagged foreground window [V]; a gene significant only when unmasked [V]. |
| Not justified | Adaptation to a particular environment or trait; that the significant gene is the cause of a phenotype; absence of selection from a nonsignificant result (power is low for small selected fractions [V]). |

### BUSTED-E (`phylo_busted --error-sink Yes`; experimental)

| | |
|---|---|
| Question answered | BUSTED's question after an error-sink class (ω ≥ 100, at most 1% of sites by default) is allowed to absorb residual alignment error [L: Minus the Error, verified; V: HyPhy 2.5.103 `BUSTED.bf` lines 207-223]. |
| Assumptions | Residual errors are rare and small; the heuristic sink can separate them from selection [L: the authors call the error component heuristic, missing many real errors and flagging some real variation, verified]. |
| Known failure modes | A whole non-homologous foreground segment: 24-28% (5 taxa) and 32% (8 taxa) false positives, rising with segment length [V: stage 2B, busted-e pilot, [forensics](../validation/busted-e-forensics.md)]; power loss of 14-25 points on clean data [V: stage 2B]. |
| Diagnostics [D] | Error-sink class weight and ω in the results (`fg_prop_1`, `fg_omega_1`). |
| Not justified | That BUSTED-E protects against misassembled, misannotated or non-homologous segments; any claim from BUSTED-E alone without the masking comparison. |

### aBSREL (`phylo_absrel`)

| | |
|---|---|
| Question answered | Per branch, whether some sites on that branch evolved with ω > 1; gene `p_value`: whether any tested branch did (Smith et al. 2015, *MBE*) [L, cited]. |
| Assumptions | As BUSTED; for `p_bh_lineage`, a lineage chosen in advance. |
| Known failure modes | Corrupted segments (99-100% per gene) [V]; realigned indels: gene `p_value` 21-26%, false discoveries in 11 of 13 and 24 of 29 lineages, every correction failing [V: study B, [multiple testing](../validation/lineage-multiple-testing.md)]; the window profile misses realignment error (0.3-1.7% flagged) [V]. HyPhy's `p_holm` lacks Holm's running maximum (Hochberg-like) [V: code audit]. Sitewise and branch-site tests are sensitive to the aligner (Markova-Raina and Petrov 2011, *Genome Res*: 50-55% false positives even with PRANK in 12 *Drosophila* genomes) [L, verified]. |
| Diagnostics [D] | as BUSTED; `lineage_family_size`, `n_branches_tested`. |
| Not justified | Counts of selected branches summed across lineages as calibrated; per-branch claims on realigned gapped genes without a sensitivity analysis. |

### PAML branch models (`phylo_branch_omega`)

| | |
|---|---|
| Question answered | Whether the foreground's ω differs from the background's (`p_two_vs_m0`); whether the foreground ω exceeds 1 (`p_fg_omega_gt1`), both averaged over all sites. |
| Known failure modes | `p_fg_omega_gt1` inconclusive at ω = 1 (5.5%, 6.8%) [V]; per-branch free-ratio ω biased low under site variation (4-11%) [V: audit E1; study L running]; corrupted segments (45-51%) [V]. Robust to realigned indels (2.7-5.3%) [V]. |
| Not justified | "Relaxed constraint" or "adaptation" from a higher foreground ω alone; ω > 1 being undetectable is not absence of positively selected sites. |

### RELAX (`phylo_relax`)

| | |
|---|---|
| Question answered | Whether selection intensity (k) differs between test and reference branches (Wertheim et al. 2015, *MBE*) [L, cited]. |
| Known failure modes | HyPhy's `p_value` anticonservative on 16 taxa (8.2-9.3%), a finite-sample property of the asymptotic test, reproduced in HyPhy directly [V: [stage 1](../validation/relax-calibration-stage1.md)]; a calibrated scaled p-value passes stage 1, confirmation running. |
| Not justified | Relaxation as evidence of trait loss; k on its own without its interval and the calibration caveat. |

### Lineage summaries (`phylo_summarize`)

| | |
|---|---|
| Question answered | Descriptive pooled ω per lineage and its interval. |
| Known failure modes | Interval undercovers and pooled ω biased low under site variation, more on long branches [V: audit E1]; full study L running (bias, coverage, false trait associations). |
| Not justified | Use in comparative trait regressions until study L reports [V]. |

## Experimental or future methods (not release-ready)

McDonald-Kreitman (owned by the MK maintainer; experimental by intent); a PRANK alignment mode (confirmation running); the adjusted enrichment mode (development running); bootstrap-calibrated RELAX p-values (confirmation running). None of these is described here as validated.

## Cross-cutting

- **Significant result ≠ adaptation.** Every test above detects a pattern of substitution; linking it to a trait, environment or function needs independent evidence [L, general; V for the alignment-error cases].
- **Nonsignificant ≠ absence.** Power is limited and measured only in simulation [V: power tables].
- **Implementation fidelity ≠ statistical validity.** dndsR reproduces HyPhy and codeml to 1e-12 or better wherever checked [V], which establishes the wrapper, not the inference.
- **Masking**: in study B's mixtures masking removed almost no raw detection of selected genes (BUSTED 36 of 72 clean selected genes both ways; aBSREL 110 against 107 of 116); the drop in BH-adjusted power came from BH's adaptivity once the invalid discoveries were removed [V, exploratory reanalysis, [forensics](../validation/busted-e-forensics.md)]. Masking can remove genuine signal when selected sites form one block (20 points in stage 2B) [V].
