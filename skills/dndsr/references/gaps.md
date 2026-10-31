# Known gaps

Current open issues are in `docs/validation/known-issues.md`. For an agent, the most important are:

- dndsR has been tested on a designed example, not validated on real genomes.
- Docker, Apptainer, macOS, and native installation are untested.
- The published image, `v0.2.0`, predates this code. It has been tested with it. No version-matched image exists yet.
- `dnds_state_contrast --level gene` is experimental.
- `go_enrichment` and `dnds_ideogram` have been confirmed to run, not independently validated.
- No DOI exists yet.
