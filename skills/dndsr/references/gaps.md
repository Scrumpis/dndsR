# Known gaps

Current open issues are in `docs/validation/known-issues.md`. For an agent, the most important are:

- Real-genome validation covers plant genomes. Other genomes may behave differently.
- `covariate_enrichment` has been calibrated on 4 plant pairs from 2 families. Polyploids and other clades are untested, and `q0.95` exceeded 5% in one setting, so calibration is not established for it. A user can run the simulation on their own table: `docs/validation/enrichment-confounding.md`.
- Docker and macOS are untested.
- A release's image is built after its tag, so a tag's documents name the previous digest. The v0.4.0 tag's documents name v0.3.0's image, which has the same dependencies (digest in `docs/validation/stack-comparison.md`), and dndsr-install puts the clone's code ahead of the copy in that image. `dev` names v0.4.0's image.
- `dnds_state_contrast --level gene` is disabled.
- The McDonald-Kreitman test is development code, not in v0.4.0. It has run on designed examples and the bundled example only, not on real population data. It uses bi-allelic SNPs only, ignores allele frequencies, counts codons only, and drops non-1:1 pairs, so it is not the full MK framework (no DoS weighting, no multiallelic handling). It does no SNP filtering itself, and it excludes 2- and 3-fold sites from both classes, so `Ps`/`Ds` are smaller and `alpha` is undefined for more genes than in tools that count those sites.
- `phylo_prepare` and `phylo_busted` are experimental development code, not in v0.4.0, and no published image includes HyPhy. They reproduce HyPhy run directly on synthetic examples; their false-positive rate and power have not been measured. Only single-copy orthogroups are tested, and alignments are not filtered. HyPhy's result depends on its random seed: on one simulated gene, the likelihood ratio ranged 7.0 to 8.8 over 20 seeds. `docs/validation/phylo-busted-wrapper.md`.
