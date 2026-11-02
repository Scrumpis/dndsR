# Known gaps

Current open issues are in `docs/validation/known-issues.md`. For an agent, the most important are:

- Real-genome validation covers plant genomes. Other genomes may behave differently.
- `covariate_enrichment` has been calibrated on 4 plant pairs from 2 families. Polyploids and other clades are untested. A user can run the simulation on their own table: `docs/validation/enrichment-confounding.md`.
- Docker and macOS are untested.
- The published image, `v0.2.0`, predates this code. It has been tested with it. No version-matched image exists yet.
- `dnds_state_contrast --level gene` is disabled.
- No DOI exists yet.
