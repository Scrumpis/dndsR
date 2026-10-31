# Release readiness: v0.3.0 candidate

**Assessment: a good candidate, not yet ready to release.** The code on `dev` passes every local check, but CI has never run on it, its container has never been built, and it has not been run on real genomes. Publishing also needs the maintainer's decisions below.

## Evidence

| Check | Result |
|---|---|
| `R CMD check --no-manual`, in the image, `lme4` absent | Status OK |
| Package tests | 44 tests, 204 expectations, 0 failed, 0 skipped |
| `dndsr verify_example` | 14 of 14 |
| Tutorial commands in a clean `HOME` | Pass |
| Fresh Claude agent, clean clone, plain request | 7 of 7 script-checked criteria |
| Git safeguards | 88 of 88 |
| Estimators against PAML, KaKs_Calculator, and simulated truth | No error found; bias documented |

## Blocking the release

| Item | Why | Owner |
|---|---|---|
| Push `dev`, then get CI green | CI on `dev` and the new workflows have never run. macOS, Windows, and other R versions are untested | Maintainer pushes, then fix findings |
| Build and test the new image | The Dockerfile has not been built, since no Docker is available here. `container.yaml` does it on push | CI |
| Run on real genomes | Everything so far is on a designed example, so problems of scale and messy annotations are untested. Comparing a pair of real genomes between `c8b534f` and this version would also show how much the CDS bug changed real results | Maintainer's data |
| Decide the items below | They need the maintainer | Maintainer |

## Decisions for the maintainer

1. **Default estimator.** Comeron stays the default, and it reads low in simulation (true ω of 2 gives a median of about 1.45 to 1.9). YN-type methods were closer, but need KaKs_Calculator in the image and would change every result. Recommendation: keep Comeron for v0.3.0, add KaKs_Calculator to the image, and decide before v1.0.0.
2. **Announcing the CDS bug.** Results from v0.1.0 and v0.2.0 can be wrong. Recommendation: add a notice to those release pages and contact known users. Check any manuscripts or talks that used earlier results. The earlier version also gave some genes wrong-genome annotations.
3. **The bundled Arial font.** `inst/fonts/ArialBold.ttf` is proprietary and redistribution rights are unconfirmed, which matters for public archiving and CRAN. Recommendation: replace it with an openly licensed font such as Liberation Sans before Zenodo archiving.
4. **Gene-level `dnds_state_contrast`.** A code audit found design problems. It is marked experimental. Recommendation: exclude it from claims in v1.0.0 unless redesigned, because the right design is a scientific choice.
5. **`main` ruleset.** The workflow now reports a check named `R-CMD-check`, which should satisfy the requirement, but this is unconfirmed until a pull request runs. The required approving review cannot come from the pull request's author, so a sole maintainer needs a bypass or second reviewer.
6. **Publishing.** Merging to `main`, tagging, pushing an image, the GitHub release, and Zenodo all await authorization. See `RELEASE.md`.
7. **`.audit` location.** It is on scratch with the clone, so it depends on regular touching. Copying it to durable storage would protect the record: `tools/dev/setup.sh --audit-dir <path>`.

## Not done

- A Codex end-to-end run: an unsandboxed Codex run was refused by this session's safety controls, and the sandboxed mode now in `test-agent.sh` has not been run. Codex does discover the skill.
- Docker, Apptainer, macOS, Windows, and from-scratch native installation are untested.
- `go_enrichment` p-values and `dnds_ideogram` content are not independently validated, only confirmed to run.
- A plugin or marketplace package for the skill.
- Items in `docs/validation/known-issues.md`.
