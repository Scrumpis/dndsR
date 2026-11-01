# Decisions and open questions

## Decisions

| Date | Decision | Evidence or reason |
|---|---|---|
| 2026-10-31 | `dev` is the working line. Planned path: audit, merge to `main`, release v0.3.0, then v1.0.0 for publication. | Maintainer instruction. `DESCRIPTION` on `dev` is already 0.3.0.000. |
| 2026-10-31 | Agents commit to `dev`; the maintainer pushes. | Maintainer instruction. The clone machine also has no GitHub credentials. |
| 2026-10-31 | `skills/dndsr/` is canonical, with symlinks from `.claude/skills` and `.agents/skills`. | Claude Code and Codex docs both support symlinked skill directories. Avoids duplicate content. |
| 2026-10-31 | No plugin manifest yet. | Nothing released to package. |
| 2026-10-31 | Git hook sources are tracked in `tools/dev/githooks`; `tools/dev/setup.sh` copies them to `.git/hooks`. | `tools/dev/test-safeguards.sh` showed that a tracked `core.hooksPath` disappears when `main` is checked out, leaving `main` unprotected. No previous hooks existed to preserve. |
| 2026-10-31 | `extract_cds` uses the longest valid isoform per gene instead of `cdsBy(by = "gene")`, and excludes CDS with in-frame stops or lengths not divisible by 3. | The old call produced mis-ordered and chimeric sequences, and dN/dS 1.615 instead of about 0.11 for a designed gene: `docs/validation/cds-extraction.md`. Longest isoform is a common convention, and `group_by = "tx"` remains. |
| 2026-10-31 | dndsR stays a real R package fit for CRAN one day, as well as a containerized CLI, for people and agents equally. | Maintainer instruction. CRAN needs GitHub-only `orthologr` to be optional, so it moved to `Suggests` with a clear runtime error. |
| 2026-10-31 | Provenance uses content hashes, never file dates. | Scratch is purged after 45 days and the maintainer may recursively touch files. |
| 2026-10-31 | The skill documents only the Singularity install and CLI discovery as verified. | See `skills/dndsr/references/verification.md`. |

| 2026-10-31 | Comeron stays the default estimator. | Changing it alters every result. Evidence of its bias is documented: `docs/validation/estimators.md`. Final decision is the maintainer's. |
| 2026-10-31 | Pin orthologr to `06c3098` and drop the patch. | Upstream absorbed the patch, which broke the Dockerfile. Example results are byte-identical. |
| 2026-10-31 | The image now includes dndsR. | Simpler for users. Not built here; CI builds it. |
| 2026-10-31 | Agents may delete files they generated in the same session without asking. | Maintainer instruction. The guard still protects `.git` and `.audit`. |
| 2026-10-31 | Deleted unused helpers, the stale vignette, and the orthologr patch. | No references in code; the vignette named files that did not exist; the patch can no longer apply. Git history keeps them. |

| 2026-11-01 | Correctness notice goes in the v0.3.0 release notes, `docs/validation/`, the migration guide, and a README notice removed in the next release. No emails, pinned issue, or Docker Hub notice. | Maintainer decision: pre-releases, with no evidence of external adoption. |
| 2026-11-01 | Published `dev` history is not rewritten. Unpublished results stay out of current files and future commits. | Maintainer decision: exposure was limited to aggregate figures. |
| 2026-11-01 | Gene-level `dnds_state_contrast` stops with an error. | Confirmed by running that it gives no meaningful result. |

## Open questions

Decisions for release are in `RELEASE-READINESS.md`.

| Question | Why it matters | Owner |
|---|---|---|
| Can the `main` ruleset be satisfied? It requires 1 approving review and a check named `R-CMD-check`, but the workflow's jobs are named like `ubuntu-latest (release)`. | May block the v0.3.0 merge. Not yet tested against a pull request. A sole maintainer also cannot approve their own pull request without a bypass. | Maintainer |
| Should v0.3.0 get a new image? The only image predates this code and was built for v0.2.0. | The skill must link a release to one digest. | Maintainer |
| Where should `.audit` live for the long term? | Scratch is purged after 45 days unless touched. Hashes in `.audit/manifest.sha256` detect loss of content but cannot restore it. | Maintainer |
