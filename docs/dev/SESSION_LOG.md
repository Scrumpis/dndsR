# Session log

Summaries only. Full transcripts are local: see `RUNBOOK.md`. Newest first.

## 2026-11-13: phylogenetic module, P2

- **Session:** `83fbd26d-6d94-4764-aee9-e34ca0f4eb72`
- **Decided by the maintainer after P1:** fixed seed and threads stay default, with extra starts as diagnostics only; HyPhy defaults kept, BUSTED-E experimental; strict identifiability, with mathematical and policy exclusions distinguished; gene-tree discordance reported, never excluding; version 0.4.0.9000; development image built and tested locally, nothing published.
- **Done:** `phylo_prepare` (genomes or CDS and long or OrthoFinder orthogroups to codon alignments, labelled trees, QC, member and exclusion records, provenance, manifest); `phylo_busted --prepared-dir`, `--extra-seeds`, `--gene-tree-mode`; outcome bases and stages; a six-genome synthetic example with designed failures; direct-HyPhy reference for species and gene trees; tests; docs, tutorial, and skill.
- **Verified:** all outcomes as designed; alignments equal the designed CDS and a true indel alignment; gene-tree mode equals direct HyPhy exactly; the tutorial ran in a clean `HOME`; the workflow passed in a development image built locally with HyPhy.
- **Found:** HyPhy's results differ slightly between processors with and without AVX-512 (2 of 8 example runs; LRT by up to 0.005), and are bit-identical within each class. Found when `check.sh` ran on an AVX2 node against references stored on an AVX-512 one. Example checks and tests now run HyPhy directly on the same machine. Also: an explicit `--hyphy` or `--mafft` path that did not exist fell back to `PATH`; fixed. An image build on scratch exceeded the file quota and failed; builds now unpack on node-local storage. One synthetic case was first classed wrongly in its design (OG12: only the stem is lost, so it is partly unidentifiable), caught before running.
- **Not run here:** `check.sh`'s safeguards and tutorial on AMD Rome nodes, where the host's git crashes with an illegal instruction (a node issue, recorded before); they pass on Intel nodes.
- **Failed or not done:** P3 not started. CI has not run. Nothing pushed, merged, tagged, or published.

## 2026-11-13: phylogenetic module, P0 and P1

- **Session:** `83fbd26d-6d94-4764-aee9-e34ca0f4eb72`
- **Done:** proposal approved by the maintainer (HyPhy engine, BUSTED first). P0 decisions in `DECISIONS.md` and `phylo-design.md`. P1: `phylo_busted`, `phylo_example`, `verify_phylo_example`, method-agnostic tree, alignment, engine, and runner layers; synthetic example; tests; container and conda definitions with HyPhy 2.5.103 (not built or published); user methods, tutorial, outputs, skill.
- **Verified:** `phylo_busted` equals HyPhy run directly on independently prepared inputs, to every digit, for all 4 testable example orthogroups. Tutorial commands ran in a clean `HOME` on the v0.4.0 image with micromamba HyPhy. `tools/dev/check.sh`: package tests 136, 0 failed, 0 skipped; tutorial and skill pass.
- **Found by running HyPhy:** unseeded runs differ; thread count changes results; given branch lengths are starting values; rooted trees are unrooted; stop codons fail a run. Over 20 seeds one gene's LRT ranged 7.0 to 8.8. HyPhy 2.5 is in maintenance mode, and HyPhy 3 is too new to use.
- **Errors along the way:** the first gene-tree rule rejected discordance inside a clade-scope foreground, which the design allows; caught by a test and fixed. A Slurm job pointed at node-local `/tmp`; resubmitted with the example in the repository. The safeguard test failed before commit because its fixture combines committed code with the uncommitted skill.
- **Failed or not done:** P2 onward not started. No image built; nothing pushed, merged, tagged, or published.

## 2026-11-06: after the v0.4.0 release

- **Session:** `a897d1b6-0934-4e10-b07a-c791f9d940eb`
- **Read:** the maintainer merged PR #4 (`957cfe1`), tagged and released v0.4.0, and `publish` pushed `sha256:f1b42d9d`, whose revision label is `957cfe1`. Zenodo version DOI 10.5281/zenodo.23903253.
- **Done:** RELEASE.md step 3: the published digest replaces v0.3.0's in user documents and the skill. The tutorial passed with Singularity 4.1.2 and Apptainer 1.5.4, and the package tests passed in the image (102 tests, 0 failed, 0 skipped). The image's R packages match the v0.3.0 image's except dndsR. The Chenopodium equivalence check found its software identical to `v0.4.0`.
- **Failed or not done:** the release has no attached assets. Nothing was pushed, merged, tagged, or published by agents.

## 2026-11-02: v0.4.0 candidate

- **Session:** `765c203d-5ba8-48a9-ad70-6681331d05db`
- **Done:** refs and rulesets read; audit archive moved and hash-checked; v0.4.0 candidate with `covariate_enrichment` accepted, warnings on `go_enrichment` and `term_enrichment`; independent reproduction with a Python refit; image stack comparison.
- **Found:** a quantile `pos_threshold` was compared as text in `covariate_enrichment`, `term_enrichment`, and `ipr_enrichment`, silently giving no positive genes; fixed. The two published images choose slightly different pairs, because tiny E-values are read differently: `docs/validation/stack-comparison.md`. `q0.95` exceeds the calibration rule in one setting.
- **Errors of the checks:** the first refit diverged and reported up to 100% false terms; R was right, and it was fixed and rerun. A git crash from a module built for another CPU left a core dump that the public-safety check caught. Three safeguard tests failed because `dev` now contains `main`, and a version mismatch in the fixture; the test fixture was corrected, and no hook changed.
- **Then:** the maintainer pushed `5983405`, and all four workflows passed. The public-genome tutorial passed at that commit on its third try, after a UniProt outage stopped the first two.
- **Failed or not done:** nothing was merged, tagged, or published by agents. Chenopodium work is outside this repository.

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
