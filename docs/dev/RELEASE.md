# Release process

Each release links one git tag, one `software-version` in the skill, one image digest, and one Zenodo DOI. Every publication step below needs the maintainer's explicit authorization.

## Before release

1. `git fetch origin`, with `dev` clean and pushed, and CI green on `dev`, including the container workflow.
2. Update `Version` in `DESCRIPTION`, `metadata.software-version` in the skill, `version` in `CITATION.cff`, and `NEWS.md`.
3. Regenerate docs with roxygen2. Run `tools/dev/test-safeguards.sh`, `tools/dev/run-tests.sh`, `tools/dev/test-tutorial.sh`, and `tools/dev/test-agent.sh claude`. Repeat the old-against-current comparison in `tools/dev/validation/` on real genomes if CDS, pairing, annotation, or enrichment code changed. Keep unpublished results outside Git.
4. Build the package and run `R CMD check` on the tarball.
5. On `dev`, point clone commands at `--branch vX.Y.Z`, and confirm versions match, so that the tag's own documents are right. Tag promptly after merging, because the clone command does not work until the tag exists.

## Publish, once authorized

1. Complete `docs/dev/GITHUB-SETUP.md`. Open a pull request from `dev` to `main`, and merge once its checks and review pass.
2. Tag `vX.Y.Z` on `main`. `publish.yaml` checks the tag, waits for approval, builds, tests, and pushes that one version tag, recording the digest. You never build on a laptop. It never touches `latest`.
3. In a follow-up commit through a pull request, replace the image digest in `skills/dndsr/`, `docs/tutorials/`, `README.md`, and `INSTALL.md` with the published one. Update `references/verification.md`, and rerun the tutorial test.
4. Publish the GitHub release, using `docs/dev/release-notes-<version>.md`, with the manifest and skill tarball that the workflow uploaded, and have Zenodo archive it. Add the DOI to `CITATION.cff`, `README.md`, and the skill.

The published image's digest cannot be inside the commit it was built from, so tutorials in the tag still show the previous digest until step 3. The release notes should give the new one.

## After v0.3.0

The maintainer re-added the correctness notice and a DOI badge to `README.md` after v0.3.0 (PR #3), and both stay in v0.4.0. Remove the notice only when the maintainer decides. Keep `docs/migrating-to-0.3.md`, `docs/validation/`, `NEWS.md`, and the v0.3.0 release notes permanently.

## CRAN (later)

Known obstacles: `orthologr`, `metablastr`, and `rdiamond` are GitHub-only and `Remotes` is not accepted on CRAN, so they need an `Additional_repositories` route or a CRAN release upstream. `inst/extdata/interpro_entry.list.tsv` is 2 MB, against a package size limit of about 5 MB. Examples must not need external tools.
