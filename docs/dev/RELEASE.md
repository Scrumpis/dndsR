# Release process

Each release links one git tag, one `software-version` in the skill, one image digest, and one Zenodo DOI. Every publication step below needs the maintainer's explicit authorization.

## Before release

1. `git fetch origin`, with `dev` clean and pushed, and CI green on `dev`, including the container workflow.
2. Update `Version` in `DESCRIPTION`, `metadata.software-version` in the skill, `version` in `CITATION.cff`, and `NEWS.md`.
3. Regenerate docs with roxygen2. Run `tools/dev/test-safeguards.sh`, `tools/dev/run-tests.sh`, `tools/dev/test-tutorial.sh`, and `tools/dev/test-agent.sh claude`.
4. Build the package and run `R CMD check` on the tarball.

## Publish, once authorized

1. Open a pull request from `dev` to `main`, and merge once its checks and review pass.
2. Build the image at the merge commit, and verify it by running the example in it, as in `.github/workflows/container.yaml`:

   ```
   docker build -f tools/Dockerfile -t scrumpis/dndsr:vX.Y.Z .
   ```

3. Push the image, and record its digest.
4. Replace the digest in `skills/dndsr/`, `docs/tutorials/`, `README.md`, and `INSTALL.md`. Replace `--branch dev` in clone commands with `--branch vX.Y.Z`. Update `references/verification.md`, and rerun the tutorial test.
5. Tag `vX.Y.Z`, publish the GitHub release, with `NEWS.md` as notes, and have Zenodo archive it. Add the DOI to `CITATION.cff`, `README.md`, and the skill.

The tag must contain step 4, so the image digest is recorded in a commit made after the image is built. The image's own copy of dndsR differs only by that documentation commit, unless the release process changes.

## CRAN (later)

Known obstacles: `orthologr`, `metablastr`, and `rdiamond` are GitHub-only and `Remotes` is not accepted on CRAN, so they need an `Additional_repositories` route or a CRAN release upstream. `inst/fonts/ArialBold.ttf` has unconfirmed rights. `inst/extdata/interpro_entry.list.tsv` is 2 MB, against a package size limit of about 5 MB. Examples must not need external tools.
