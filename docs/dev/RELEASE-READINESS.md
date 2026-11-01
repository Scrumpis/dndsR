# Release readiness: v0.3.0 candidate

**Recommendation: NO-GO until CI passes on the commits after `b59cf23`.** See the evidence below and `GITHUB-SETUP.md` for steps only the maintainer can do. Unpublished real-genome results are kept outside Git; this report covers them only in general terms.

## CI at `b59cf23`, the first push

| Workflow | Result | Cause or note |
|---|---|---|
| `skill` | Pass | |
| `install` (clean conda, example) | Pass | |
| `container`: real `docker build` | **Pass** | First real Docker build of this Dockerfile |
| `container`: example in the image | Fail | Very likely a workflow bug: Docker created the mounted folder as root, and the container ran as another user. Logs need a login, so this was inferred. Fixed |
| `R-CMD-check`, all platforms | Fail before any package check | `upgrade: true` is no longer accepted by the setup action. Fixed |
| `R-CMD-check` rollup job | Reports under the name the ruleset requires | Good sign for the required check |

Later results replace this table.

## Blocking the release

1. **Green CI on the fixed commits,** including macOS and Windows `R CMD check`, which have never run on this code. Expect some findings.
2. **Maintainer decisions** below.

## Decisions for the maintainer

See the session report and `GITHUB-SETUP.md`: ruleset changes, Docker Hub setup, and final release approval. Release notes: `release-notes-0.3.0.md`. Published `dev` history stays as is, by the maintainer's decision.

## Not done or not tested

- macOS and Windows runtime use.
- R from CRAN and Bioconductor failed on the tested HPC module, because of `rtracklayer` and OpenSSL. Conda is the recommended route.
- **Codex:** in its sandbox, Singularity could not start (no `/dev/fuse`), so the agent unpacked the image and ran it under Bubblewrap. It reported passing all 14 checks, which our disk checks could not confirm, so this is not a clean pass.
- Docker Hub push: cannot be tested without publishing.
- Isoform effects on dN/dS were not compared end to end on real data. CDS were.
- Open items: `docs/validation/known-issues.md`.
