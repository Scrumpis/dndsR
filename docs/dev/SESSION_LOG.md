# Session log

Summaries only. Full transcripts are local: see `RUNBOOK.md`. Newest first.

## 2026-11-02: v1.0.0 validation, overnight

- **Session:** `61d9f4da-7419-4c58-bb02-9ebdad0c13e6`
- **Done:** estimator benchmark, RBH against OMA, enrichment methods compared by simulated evolution under a prespecified rule, `covariate_enrichment`, confirmation on new seeds and a held-out tomato and potato pair, `calculate_dnds` notices, faster `append_annotations` and `ipr_enrichment` with identical output, public-genome tutorial run from a commit, publication checklist.
- **Found:** `go_enrichment` defaults report false terms in 43% of simulated runs at AtAl, against the 9.5% documented; `extract_cds` stopped on trans-spliced plastid *rps12*, now with a clear error; short lists move with annotation release.
- **Disclosed:** OMA's API returned 403 to Python's default user agent, so requests used an identifying user agent with pacing. PMC's reCAPTCHA was not circumvented; Europe PMC served that text. See `docs/validation/orthology.md`.
- **Failed or not done:** the first filter calibration array failed: `singularity --env` splits on the comma in the expression, and every task exited 1. It was fixed and rerun; 1,400 tasks completed. Nothing pushed, merged, tagged, or published. CI has not run on these commits. Chenopodium was not run.

## 2026-11-02: final enrichment defaults

- **Session:** `61d9f4da-7419-4c58-bb02-9ebdad0c13e6`
- **Start:** `40ff597`. **Software:** `eed1096`.
- **Done:** `ipr_enrichment` global and `go_enrichment` BH by default, with selectable alternatives. Null models rerun. Independent recalculation, including GO. Unrelated outputs byte-identical. Earlier evidence preserved.
- **Found by testing:** GO's default is 9.5% under `strata`, mostly GO:0005576, which length and dS largely reproduce; reported, not corrected. An ontology with no positive genes crashed topGO; fixed. `man/` usage lines were generated from an older install and showed the old defaults, while the code was right; regenerated.
- **Errors of the validation itself:** the end-to-end job and parent-child null shared an output folder, which stopped 10 of 100 `strata` runs; left out and counted.
- **Then:** maintainer accepted the GO limitation and pushed `f2da1f1`. All four workflows passed there.
- **Failed or not done:** nothing merged, tagged, or published. No biological comparison of terms. Not measured: false discovery with real associations.

## 2026-11-01: enrichment correction

- **Session:** `61d9f4da-7419-4c58-bb02-9ebdad0c13e6`
- **Start:** `5d53ced`. Maintainer approved the three fixes.
- **Done:** BH before reporting filters; GFF comment parsing; GO crash. Null models rerun with real functions, before and after. Independent recalculation, GFF recovery check, and evidence in `docs/validation/enrichment.md`.
- **Found by testing:** 73 pairs lost terms, not 71, which had counted one column at a time. Default `ipr_enrichment` is still 9–17% across tables, because it corrects per type; documented. `entries_source = "local"` still needs network unless release detection is off; recorded as open. An early parent-child run would have counted that network error as a clean run: stopped, fixed, and rerun.
- **Failed or not done:** nothing pushed, merged, tagged, or published. CI has not run on this code. No biological comparison of terms.

## 2026-11-01: enrichment validation and faster tests

- **Session:** `61d9f4da-7419-4c58-bb02-9ebdad0c13e6`
- **Done:** enrichment calibration on the public pair, with real functions and three nulls. Tests run in parallel, with `tools/dev/check.sh`.
- **Found:** `min_pos` before BH is badly anti-conservative; GFF `#` truncation; `include_definition = FALSE` crash. None is fixed: `dev` is the release source, and these need a decision.
- **Failed or not done:** nothing pushed, merged, tagged, or published. The literature comparison of the top terms has not been done.

## 2026-11-01: v0.3.0 release readiness

- **Session:** `61d9f4da-7419-4c58-bb02-9ebdad0c13e6`
- **Done:** public-genome validation, estimator characterization, Docker build and Apptainer fixes for CI, and release records. All workflows pass at `f84521f`.
- **Failed or not done:** the Codex agent test stays partial and was not rerun. Nothing was merged, tagged, or published by agents.

## 2026-10-31: autonomous development toward v0.3.0

- **Session:** `61d9f4da-7419-4c58-bb02-9ebdad0c13e6` (same as the bootstrap)
- **Done:** see `PROGRESS.md`. The major results are the CDS assembly bug in every release and the annotation leak, both fixed with tests.
- **Delegated:** read-only scientific and software audits. Their reports are in the transcripts. Findings were checked before fixing; some remain code-review findings only (`docs/validation/known-issues.md`).
- **Failed or not done:** see `RELEASE-READINESS.md`. Nothing was pushed, merged, tagged, or published.

## 2026-10-31: Phase 0 bootstrap

- **Session:** `61d9f4da-7419-4c58-bb02-9ebdad0c13e6`
- **Start:** fresh clone, clean `dev` at `c8b534f`
- **Done:** git safeguards, Claude Code hooks and permissions, audit trail, four development agents, end-user skill and tutorial, validation tooling, and project docs. No scientific code changed.
- **Verified:** see `skills/dndsr/references/verification.md` and `tools/dev/test-safeguards.sh`.
- **Failed or not done:** this machine has no GitHub credentials, so a push dry-run failed; the maintainer now handles pushes. No analysis was run.
- **Found by testing:** the first hook design left `main` unprotected once it was checked out, missed merges from `main`, and blocked `git pack-refs`. All were redesigned and retested. A real subagent run also showed that plan mode and a frontmatter hook left the auditors able to write through Bash, so a project-level hook now enforces read-only Bash; the probe left one empty file, which was removed. The guard also falsely denied a command that wrote documentation mentioning a blocked flag; those files were written with the file tools and the behavior is documented.
