# Study E: a public multispecies dataset, end to end (design draft)

Status: **draft for the maintainer's approval; not run; species and accessions to be verified before approval.** Gap 6 of the [audit](validation-audit.md), and the final integrated test the maintainer asked for before release: the whole workflow, on real public genomes that vary in genome size and transposable-element content (the lineage use case), with real OrthoFinder output, including the pairwise estimators, orthology, and enrichment.

## What this study can and cannot show

It can show that every documented command runs on real data at realistic scale; that its diagnostics behave (flag rates, exclusions, failures, warnings); that results are reproducible (reruns, thread counts, image); that independent computations agree where they must (consistency checks below); and what it costs. **It cannot validate biological conclusions, and agreement with biological expectation (for example, more relaxed constraint in larger genomes) is not counted as evidence of correctness.** Calibration and power come from the simulation studies.

## Dataset (candidates; to be confirmed against NCBI before approval)

Brassicaceae, which the existing public fixtures use (`docs/validation/public-genome.md`), so that results connect with the pairwise and BUSTED evidence: *Arabidopsis thaliana*, *A. lyrata*, *Capsella rubella*, *Eutrema salsugineum*, *Brassica rapa* (already used), plus candidates with small and large, TE-rich genomes, for example *Schrenkiella parvula*, *Arabis alpina*, and *Thlaspi arvense*, whichever have public chromosome-level assemblies with gene annotations under terms that allow use. Genome sizes and TE content are taken from the assemblies' own reports, not asserted here. Polyploid histories (*Brassica*'s triplication) are kept deliberately: they exercise `multi_copy` exclusion and single-copy ascertainment on real data.

## Workflow (every step as documented for users; commands recorded)

1. CDS extraction for every genome (`extract_cds`), with organelle removal where needed.
2. Pairwise: `calculate_dnds` of *A. thaliana* against each species, default estimator and the alternatives in [estimators.md](../validation/estimators.md); `append_annotations`; `covariate_enrichment`, `go_enrichment`, `term_enrichment`, `ipr_enrichment` on the *A. thaliana*-*A. lyrata* pair.
3. Orthology: OrthoFinder (version pinned in the image) on all proteomes; `phylo_prepare` from its output with `orthogroup_format = "orthofinder"`.
4. Lineage: `phylo_branch_omega` (M0, free-ratio, two-ratio with a foreground chosen before the run, for example the large-genome lineage), `phylo_absrel`, `phylo_relax` (flagged), `phylo_busted`, `phylo_summarize` with the species tree.
5. The skill's own agent test on the same commands.

## Checks (fixed before running)

- **E1, completion**: every step finishes; every exclusion has a documented reason; failures are listed with causes.
- **E2, consistency between independent computations** (these must agree whatever the biology): pairwise dS between *A. thaliana* and each species against the sum of free-ratio branch dS along the path between them in the same orthogroups (median ratio within 0.8-1.25; differences explained by the models' codon frequencies are reported); RBH orthologs against OrthoFinder single-copy orthogroups (share of pairs shared, against the 3.6-15% OMA disagreement already measured); `pooled_omega` recomputed outside dndsR from the branch tables (exact).
- **E3, engine equivalence on real data**: 50 orthogroups chosen by seed before the run, codeml and HyPhy run directly on dndsR's inputs; statistics within 1e-3 (as in `phylo-lineage.md`).
- **E4, reproducibility**: a second run of every command, and at 1 and 8 threads, with byte-identical tables (or differences listed and explained, as for RELAX's optimizer).
- **E5, diagnostics on real data**: rates of `high_ds`, `few_synonymous_changes`, `dissimilar_window`, `multi_copy`, and missing taxa per species, against the simulation studies' rates; any flag above 20% of genes is examined.
- **E6, cost**: core-hours per step, for the documentation.

## Cost and seeds

Projected from measured times (codeml M0 + free-ratio about 134 core-seconds per gene on 8 taxa; aBSREL about 20; RELAX about 22; OrthoFinder and alignment about 50 core-hours): about 10,000 single-copy orthogroups, **about 600 core-hours**, most of it codeml. Seeds: dndsR's defaults (the workflow as users run it); E3's sample seed 20261201.

## Data handling

Public genomes only; results that are aggregate summaries of public data may be committed as validation evidence, as the existing public-genome reports are. Nothing from unpublished data enters this study.
