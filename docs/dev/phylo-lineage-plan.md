# Lineage-level selection: architecture review and plan

Review only, 2026-11-20. Nothing here is implemented, and none of it is in v0.5.0. The maintainer set the order on 2026-11-20: PAML branch models, then lineage summaries, aBSREL, and RELAX (`DECISIONS.md`). The motivating question: compare lineage-level selection with LTR activity and genome expansion or contraction across species, accounting for branch length, sampling, and phylogenetic nonindependence. Read with [phylo-design.md](phylo-design.md).

## What already exists and can be reused

Checked by reading `R/phylo_*.R` at `3943ffc`.

| Need | Existing piece | Reusable as is? |
|---|---|---|
| Orthogroups to codon alignments, QC, provenance | `phylo_prepare`, `R/phylo_alignment.R` | Yes |
| Species tree, pruning per orthogroup, branch identity by bipartition, identifiability | `R/phylo_tree.R` (`.phylo_map_foreground`, branch map with `taxa_side`, `species_tree_branches`) | Yes; branch identity by bipartition is what lets a branch be followed across orthogroups |
| Isolated HyPhy runs with fixed seed and threads, JSON reading, extra starts | `R/phylo_hyphy.R`, `.phylo_start*` in `R/phylo_run.R` | Yes, for HyPhy methods |
| Method registry (settings, arguments, parser, unit of test) | `.hyphy_methods` | Yes; the registry already has a `unit` field (`gene` or `branch`) |
| Per-orthogroup runner, reuse, statuses, correction, provenance | `.phylo_run` | Partly: it writes one result row per orthogroup; a branch-unit method needs a branch-level table and its own correction |
| Simulation and calibration machinery | `tools/dev/validation/phylo-sim.R`, stage 2 scripts | Partly: built for one five-taxon tree and gene-level truth |

## What is actually missing

- **`phylo_absrel`.** A registry entry (aBSREL is in HyPhy 2.5.103) and an exported function are small. The real work is a branch-level results table keyed by orthogroup and species-tree branch, a correction policy (aBSREL corrects within a gene by Holm-Bonferroni; across genes is a separate family question), and mapping HyPhy's branch names back to species-tree bipartitions, including branches merged by pruning.
- **`phylo_relax`.** RELAX needs two branch sets, test and reference. The foreground file and `.phylo_write_newick` label one set (`FG`), so a second label, identifiability of both sets after pruning, and the parse of `K` are new. RELAX's K is sensitive to model misspecification, including multinucleotide changes, which needs calibration like BUSTED's.
- **`phylo_branch_omega`.** PAML `codeml` branch models (free-ratio, or two-ratio with marked lineages) are a second engine: a control-file writer, PAML's `#1` tree labels, an output parser, its own optimizer variability (multiple starting ω), and PAML in the image and `environment.yml`, which is a new dependency needing approval. HyPhy 2.5.103 has no built-in per-branch ω command (`FitModel.bf` fits one global ω). aBSREL's output may carry per-branch baseline ω, which would avoid PAML, but this has not been checked.
- **`phylo_summarize`.** Aggregating per-branch results across orthogroups needs only the branch-level table. Missing: a stable key per species-tree branch, rules for orthogroups where that branch was merged or absent, and diagnostics per branch: number of orthogroups, total dS or synonymous changes, saturation, and how many estimates hit bounds.

## Scientific constraints to design for, not after

- **Branch length.** Per-branch ω is noisy on short branches (few synonymous changes) and biased on long ones (saturation, the same compression as pairwise Comeron at high dS). Summaries must carry dS and substitution counts, and flag or weight, never average raw ratios.
- **Sampling.** Each branch is estimated from a different set of orthogroups, because of missing taxa and pruning. Branch comparisons must report and, where possible, condition on that set.
- **Ascertainment.** Only single-copy orthogroups are tested. Genome expansion and contraction are about gene families that are rarely single-copy, so lineages with more duplication are under-represented. This bears directly on the LTR and genome-size question and must be stated or addressed (duplication-aware orthology is an open v1.0 question).
- **Nonindependence.** Branch-level traits share ancestry. dndsR should export a tidy branch table (species-tree branch, ω and its uncertainty, dS, counts) that fits phylogenetic comparative tools such as phylogenetic generalised least squares in established R packages (for example `phylolm`, `caper`), rather than reimplement comparative regression. LTR activity and genome size per branch (for example from ancestral reconstruction) come from outside dndsR.

## HyPhy alternative to PAML branch models (checked 2026-11-20)

Read in HyPhy 2.5.103's own analysis files (`share/hyphy/TemplateBatchFiles`), not run:

- **aBSREL's baseline stage** fits "MG94xREV with separate omega for each branch": one dN/dS per branch, no site-to-site variation, a local model (`terms.local`), and writes each branch's "Baseline MG94xREV omega ratio" to its JSON. That is the closest HyPhy analogue of PAML's free-ratio model (`model = 1`).
- It is **not the same model.** MG94 and GY94 (codeml) define ω on different scales and differ in codon frequencies (MG94's nucleotide-position frequencies against codeml's `CodonFreq` options), so the same data can give different per-branch ω. Equivalence would have to be shown, not assumed, on simulated alignments with known per-branch ω and on real orthogroups, with agreement criteria fixed in advance.
- PAML's **two-ratio and multi-ratio branch models** (`model = 2`, with marked lineages) have no standard HyPhy counterpart: `FitModel.bf` and `SingleOmega.bf` fit one ω for the whole tree, and `TestBranchDNDS.bf` is a legacy interactive script.

So PAML stays first, as decided. The aBSREL baseline is a candidate cross-check, used as a substitute only after an equivalence study.

## Recommended order (a recommendation, not a decision)

1. **`phylo_branch_omega`**, after deciding the engine: PAML as specified (independent, standard for branch models; new dependency) or HyPhy if aBSREL's baseline ω suffices (no new engine). Either way, validate against the other implementation on simulated data with known per-branch ω, the same way pairwise estimators were checked against PAML and KaKs_Calculator. This is the quantity the comparative question needs most.
2. **`phylo_summarize`**, on top of 1: aggregation, diagnostics, and the exported branch table. Validated by independent recomputation of the summaries and by simulations where branch ω differs by design.
3. **`phylo_absrel`**: the branch-level table and runner change also serve 1 and 2. Calibrate in simulation with branch-level truth, including multinucleotide changes and alignment errors, as BUSTED was.
4. **`phylo_relax`**: two-set hypotheses, then calibration of K.

Before 1: decide the engine, the phylogenies to calibrate on (more taxa than the five-taxon tree; ideally resembling the study's), and whether duplicated genes must be included. Each new method stays experimental until calibrated.
