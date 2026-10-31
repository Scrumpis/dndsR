# Progress log

Durable checkpoints for resuming autonomous work. Newest first. Evidence is in `docs/validation/`. The plan is at the end.

## Environment

- Runtime files live outside the clone, in a `dndsr-runtime` folder beside it: image, `SHA256SUMS`, R library, reference tools (PAML 4.10.10 and KaKs_Calculator 2.0.1 in a micromamba env). Recreate them as described in `RUNBOOK.md` if purged.
- Scratch is purged after 45 days, and the maintainer may recursively touch files. **File dates mean nothing here: identify everything by content hash.**
- The maintainer pushes. Agents report the push command.

## Plan

Done: audit, core fixes, statistical fixes, CLI and failure handling, example and self-check, documentation, skill, container and installer, CI, release materials. Remaining work needs the maintainer or CI: see `RELEASE-READINESS.md`.

## Checkpoints

Commit SHAs are in `git log`. Each subject below is a commit subject.

### Add end-to-end agent test and release materials
A fresh Claude agent in a clean clone passed 7 of 7 script-checked criteria. `R CMD check`: Status OK. Tests: 44, 204 expectations, none failed or skipped. Safeguards: 88 of 88. An unsandboxed Codex run was refused by session safety controls, so `test-agent.sh` now uses Codex's sandbox and that run remains to be done.

### Make dndsR easy to install, verify, and understand
`dndsr example` and `verify_example`; docs; skill; pinned Dockerfile (not built here); installer fixes; CI. Pinned orthologr gave byte-identical results to the image's on the example. Tutorial passes in a clean `HOME`.

### Fix annotation leakage and statistical errors in downstream steps
From the scientific audit; each fix has a test. Estimator comparison: `docs/validation/estimators.md`.

### Make failures visible and stop the CLI misreading its input
From the software audit.

### Fix CDS assembly for minus-strand and multi-isoform genes
Fixed the CDS assembly bug (`docs/validation/cds-extraction.md`), added the designed example, independent NG86 reference, first tests, input hashes, and hash-based audit archiving.
