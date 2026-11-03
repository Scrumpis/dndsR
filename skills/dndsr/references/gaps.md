# Known gaps

Current open issues are in `docs/validation/known-issues.md`. For an agent, the most important are:

- Real-genome validation covers plant genomes. Other genomes may behave differently.
- `covariate_enrichment` has been calibrated on 4 plant pairs from 2 families. Polyploids and other clades are untested, and `q0.95` exceeded 5% in one setting, so calibration is not established for it. A user can run the simulation on their own table: `docs/validation/enrichment-confounding.md`.
- Docker and macOS are untested.
- A release's image is built after its tag, so a tag's documents name the previous digest. v0.4.0's documents name v0.3.0's image, `sha256:feb5e4b1cebe574bb0a1962d0589a3325ddf5d12fa599def72d2b2ddd24d7e11`, tested with this code; dndsr-install puts this code ahead of the copy in that image.
- `dnds_state_contrast --level gene` is disabled.
