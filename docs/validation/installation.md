# Installation and container validation

Only what actually ran. Scripts are in `tools/dev/validation/` and `tools/dev/`.

| Route | Environment | How | Result |
|---|---|---|---|
| Published image with `tools/dndsr-install` | Singularity CE 4.1.2, Linux amd64 HPC | `tools/dev/test-tutorial.sh`, every tutorial command, clean `HOME` | Pass, 14 of 14 |
| Same | Apptainer 1.5.4 (rootless, from conda-forge), Linux amd64 HPC | Same, with only Apptainer on `PATH` | Pass, 14 of 14 |
| Candidate image with dndsR inside | Built with Apptainer from the Dockerfile's `RUN` steps, base image pinned; R 4.6.1, orthologr 0.4.3 | `dndsr example`, `extract_cds`, `calculate_dnds`, `verify_example` in the image | Pass, 14 of 14. Example dN and dS match the reference exactly, even with a newer R and Bioconductor |
| Conda (`tools/environment.yml`) | R 4.5.3 from conda-forge and bioconda, micromamba 2.5.0, clean `HOME` | `tools/dev/validation/test-conda.sb` | Pass, 14 of 14, CLI and R |
| R from CRAN and Bioconductor on an HPC R 4.5.1 module | Clean library, following `INSTALL.md` | `tools/dev/validation/test-native.sb` | **Fail:** `rtracklayer` did not compile against the module's OpenSSL, and `rsvg` needs librsvg. Conda is recommended instead |
| Fresh Claude agent | Clean clone, plain request | `tools/dev/test-agent.sh` | See `docs/dev/RELEASE-READINESS.md` |
| Fresh Codex agent, sandboxed | `workspace-write` sandbox, no network, no `/dev/fuse` | `tools/dev/test-agent.sh` | Partial: Singularity could not start in its sandbox, so it unpacked the image and ran it under Bubblewrap. It reported all 14 checks passed and gave correct caveats, and our on-disk checks could not confirm this because of the different install path. It found two defects, both fixed |
| `R CMD check` | Published image | `R CMD build`, `R CMD check --no-manual` | See `docs/dev/RELEASE-READINESS.md` |

## Not tested

- **Docker itself.** No Docker or Podman is available here. The Dockerfile's steps were run through Apptainer, but `docker build` has not parsed the file, `.dockerignore` was not applied, and a real build needs CI.
- macOS and Windows.
- GitHub Actions: this machine has no GitHub credentials, so nothing has been pushed. `.github/workflows/` has never run on this code.
- Offline use in full. `ipr_enrichment` now falls back to its bundled entry list, tested by blocking the network; parent-child mode still needs a tree file offline. The installer needs the image and clone already present.

## Building the candidate image without Docker

`tools/dev/validation/image.def` builds the same steps. Requires a user namespace, and Apptainer for `--fakeroot` without a subuid mapping.

```
python3 tools/dev/validation/dockerfile-steps.py tools/Dockerfile > steps.sh
sed -e "s#{{ REPO }}#$PWD#" -e "s#{{ STEPS }}#$PWD/steps.sh#" tools/dev/validation/image.def > image.def
apptainer build --fakeroot dndsr-candidate.sif image.def
```

Build from a clean checkout such as a `git archive` folder, so local files like `.audit` are not copied in.
