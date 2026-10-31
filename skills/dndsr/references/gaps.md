# Known gaps

Later development phases should resolve these. None is an audit result.

## Needed for the skill

- A small example dataset with an independently established expected result.
- Validated runs of every analysis command, with input and output descriptions.
- Output column definitions and interpretation guidance grounded in the audit.
- Tested Docker, Apptainer, and native paths.
- A released image matching the released package, and a compatibility table.
- A citation and DOI.

## Noticed while bootstrapping

- The README uses `dnds_contrasts` and `dnds_state_contrasts`; the CLI has `dnds_contrast` and `dnds_state_contrast`.
- The README references `dndsR_usage_vignette.Rmd`, which does not exist. `inst/extdata/dndsR-vignette.Rmd` names example files that are not in the repository.
- `dndsr --help` shows `--comparison-file <int>`, although it takes a path, and omits `--output-dir`.
- The `tools/dndsr-install` help uses the image tag `0.2.0`; the published tag is `v0.2.0`.
- The README setup step `bash source ~/.bashrc` is not a valid command.
- The image's `lme4` (a suggested dependency) and `roxygen2` are absent.
- The Dockerfile's base image, `ORTHOLOGR_REF`, and `MAMBA_VER` float by default.
