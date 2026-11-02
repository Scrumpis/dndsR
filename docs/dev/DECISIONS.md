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
| 2026-11-01 | BH covers terms meeting `min_total` and `max_prop`; `min_pos` is report-only. GFF comments only at line start. | Maintainer approved. False-positive rates and before/after values are in `docs/validation/enrichment.md`. |
| 2026-11-01 | Superseded below. `ipr_enrichment` keeps `adjust_scope = "per_type"` as its default, documented with its rate across tables. | Not a defect: each table is near nominal, and changing a default was not part of the approval. Maintainer may change it. |
| 2026-11-02 | `ipr_enrichment` defaults to `adjust_scope = "global"`, for Fisher and parent-child; `per_type` stays selectable. | Maintainer decision. A list across tables should be corrected across them. Parent-child's elimination decisions stay within each type, and the reported `p_adj` is global. |
| 2026-11-02 | `go_enrichment` defaults to `p_adjust = "BH"` and `adjust_scope = "global"`; `per_ontology` and `none` are selectable. | Maintainer asked for correction by default unless a compelling reason opposed it. topGO's authors do not oppose it. The family is one across ontologies, as for InterPro: lower error (9.5% against 20.5% under `strata`), and tables can be pooled. |
| 2026-11-02 | No term-specific or length/dS adjustment for GO's excess under `strata`. Reported instead. | Maintainer: no inadequately validated correction; confounding analysis stays on the v1.0.0 roadmap. |
| 2026-11-02 | GO at 9.5% under `strata` accepted as a documented limitation for v0.3.0. No terms removed, no thresholds changed, defaults kept. | Maintainer decision: quantified, and a correction needs confounder modeling and biological validation for v1.0.0. After v0.3.0, priorities are length and dS confounding, independent estimator benchmarking, and empirical calibration of enrichment. |
| 2026-11-01 | No correction for length or dS in v0.3.0. | Maintainer decision: stays on the v1.0.0 roadmap until validated. |
| 2026-11-02 | `covariate_enrichment` (Poisson-binomial tail, logistic on `ns(log length, 3) + ns(log dS, 3)`, floor a tenth of the positive rate) added as an opt-in command. **Awaits maintainer approval as a new public interface.** No default changes. | Selected by a rule fixed before results, confirmed on new seeds and a held-out pair: `docs/validation/enrichment-confounding.md`. |
| 2026-11-02 | No degrees-of-freedom or floor argument on `covariate_enrichment`. | Other settings were the ones that failed calibration; an argument would invite tuning after results. |
| 2026-11-02 | No bootstrap interval feature. | Coverage measured and not established: `docs/validation/estimators.md`. Roadmap condition not met. |
| 2026-11-02 | `go_enrichment`, `ipr_enrichment`, and `term_enrichment` keep their defaults and are not removed. | Changing them is the maintainer's decision. The excess of false terms is documented. |
| 2026-11-02 | Candidate selection and enrichment thresholds stay as they are. No terms removed. | Standing maintainer instruction. |
| 2026-11-02 | `covariate_enrichment` accepted as a documented public command, with limits stated. `go_enrichment` and `term_enrichment` warn on every run and record it, with defaults unchanged. Next version is v0.4.0. | Maintainer decision. |
| 2026-11-02 | The notice goes in a `_notice.json` beside the tables, not in a new column. | Tables stay byte-identical, so existing scripts and comparisons still work. |
| 2026-11-02 | `pos_threshold` text that cannot be used now stops `term_enrichment`, `ipr_enrichment`, and `covariate_enrichment`, which compared it as text and found no positive genes. `covariate_enrichment` honors quantiles, as `go_enrichment` does. | A silent wrong result. No valid input changes. Regression tests added. |
| 2026-11-02 | `q0.95` is not recommended for a primary list. | 11.0% at AtAl InterPro `strata5`, above the rule: `enrichment-confounding.md`. |
| 2026-11-02 | Correctness notice and DOI badge stay in the README for v0.4.0. | Maintainer added them after v0.3.0, and asked that neither be lost. |
| 2026-11-02 | The audit archive moved to `~/.local/share/dndsr-audit`, mode 700. The original stays as `.audit.moved.*`. | Maintainer instruction. 27 files hash-identical in both copies. |
| 2026-11-02 | v0.4.0's documents point at the v0.3.0 image's digest. | A tag's documents cannot name an image built from it. dndsr-install puts this code ahead of the copy in the image. |

## Open questions

Decisions for release are in `RELEASE-READINESS.md`.

| Question | Why it matters | Owner |
|---|---|---|
| Should `go_enrichment` and `term_enrichment` change defaults, now that they warn? | Up to 65-70% false-term rates at one pair. | Maintainer |
| v0.3.0's release page shows `<from the publish workflow>` for its digest. | It is `sha256:feb5e4b1cebe574bb0a1962d0589a3325ddf5d12fa599def72d2b2ddd24d7e11`. | Maintainer |
| The Docker Hub `latest` tag still points at the v0.2.0 image, and `publish.yaml` never moves it. | Someone who pulls `latest` gets older dependencies. | Maintainer |
