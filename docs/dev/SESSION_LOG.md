# Session log

Summaries only. Full transcripts are local: see `RUNBOOK.md`. Newest first.

## 2026-11-23: study B complete

- **Session:** `36176ebf-ed04-471c-8141-7042748b4b80` (continued)
- **Done:** remaining arms summarised and added to the study B report: 16-taxon E0, E1 and BE3h for aBSREL and branch models, every BUSTED arm, the BUSTED mixture arms; direct reruns of the five cells outside tolerance (all agree); BE3h-T16 true-alignment diagnostic (calibrated). Docs, phylo-methods and skill updated with the 16-taxon and BUSTED findings.
- **Verified:** 12 regenerated units byte-identical to their earlier units; exclusion script with default arguments reproduces the core mixture table byte for byte; `tools/dev/check.sh` on Intel: 251 tests, 0 failed, 7 skipped (engines not in the image).
- **Not done:** maintainer decisions on the pilot and multiple-testing proposals; full L and R confirmation held; push.

## 2026-11-22 (evening): alignment pilot and multiple-testing study reported

- **Session:** `36176ebf-ed04-471c-8141-7042748b4b80` (continued)
- **Done:** alignment pilot ([report](../validation/alignment-pilot.md)): PRANK codon mode cut aBSREL's `p_raw` false positives on BE3h genes from 7.3% to 3.6% and gene `p_value` from 18.7% to 2.7% with no power loss, MACSE to 5.4% and 9.3%; neither nominated, both fail the preregistered cost criterion (76-83 and 25-29 times MAFFT); aligner disagreement is not a usable per-gene flag. Multiple-testing study ([report](../validation/lineage-multiple-testing.md)): on error-free mixtures `p_bh_lineage` and gene `p_adj` control false discoveries; on realigned indels every aBSREL output fails (FDP 47-60%), while `p_adj_two_vs_m0` controls (1.4%). Documentation of which output answers which question; `p_holm`/`p_value` wording corrected in roxygen, `man/phylo_absrel.Rd` (edited by hand to match, roxygen2 not in the image), a code comment, and the provenance `testing.within_gene` string.
- **Verified:** pilot null genes byte-identical to study B's BE3h (150/150); `tools/dev/check.sh` on an Intel node 251 tests, 0 failed, 7 skipped (an AMD node gave `Illegal instruction` from git in the safeguard test, and 16 GB ran out of memory).
- **Corrected:** I had told the maintainer PRANK was within the cost rule, taking MAFFT to cost about 8 core-seconds per gene; it measured 0.64-0.75.
- **Not done:** study B remaining arms (array 21578592, 8 of 15 tasks still running); proposals in both reports await the maintainer; full L and the R confirmation held; push.

## 2026-11-22 (later still): maintainer's rulings on study B; three runs launched; corrections

- **Session:** `36176ebf-ed04-471c-8141-7042748b4b80` (continued)
- **Done:** rulings recorded (`6d108a9`); remaining study B arms amended (BE3h for unrun indel arms) and launched; alignment pilot designed, smoke-tested and launched; p-value audit and multiple-testing mixture design, launched; `Alignment check` warning with tests and docs (`820abf3`).
- **Verified:** `tools/dev/check.sh` 251 tests, 0 failed, 7 skipped (engines not in the image); `test-phylo-lineage-engines.R` with HyPhy and codeml, 0 failed.
- **Found, and corrected:** (1) `p_bh_lineage` was judged by its per-test rate; in BE3h it made false discoveries in 11 of 13 lineage families and the gene `p_adj` 30 in 300 null genes. (2) Study B's mixture used aBSREL's `foreground` column (TRUE for every branch with all branches tested) as truth: per-branch FDP 66.3% was 74.4%, power 16% was 53%; classifications unchanged. (3) HyPhy's Holm omits the running maximum, so the gene `p_value` is Hochberg's global test; documented.
- **Not done:** results of the three runs; full L and the R confirmation (held); push.

## 2026-11-22 (later): study B core report, BE3h

- **Session:** `36176ebf-ed04-471c-8141-7042748b4b80` (resumed from `49fd39f`)
- **Done:** BE3h (E3 with homologous insertions, seed 304) summarized; direct HyPhy reruns of its 20 most significant genes; aBSREL on its true alignments as a diagnostic; core report completed (`docs/validation/lineage-robustness-core.md`) with tables in `docs/validation/data/`.
- **Found:** BE3h `p_raw` 7.1% and gene `p_value` 21% after realignment, 2.3% and 3.0% on the true alignment: the excess is alignment error, not gaps. Direct reruns 260 of 260 agree (5e-14). `p_bh_lineage` 0.8%, `p_two_vs_m0` 5.3%.
- **Not done:** the four recommendations (maintainer); 16-taxon units and BUSTED arms (held); which aligner step causes the error, not tested; push; nothing merged, tagged, or published. Documentation only, so no package tests were rerun.

## 2026-11-22: validation studies N, L, B, R stage 1

- **Session:** `d6fcae14-1433-4af7-9831-a319a98cd839` (continued)
- **Done:** designs for N, L, B, R stage 1 fixed before running, and drafts for P and E (`760ca03`); study N run and reported; L and B pilots with gates; B paralog generator fixed and re-piloted; R stage 1 run and reported with its go/no-go.
- **Verified:** fast gate on intel18 at each commit, 243 tests, 0 failed (6 skips: slow or fast-mode).
- **Found:** the enrichment excess is confounding, not FDR failure; RELAX's excess is HyPhy's, present in 2.5.62/103/104, finite-sample, shared within a dataset; a shared scaled chi-square or pooled bootstrap p restores size with size-adjusted power kept; L gate 3 failed narrowly; B costs about 4 times its budget, mostly BUSTED.
- **Not done:** full L (held), B's 16-taxon and BUSTED arms (held), R confirmation (not approved), P and E (drafts); push; nothing merged, tagged, or published.

## 2026-11-21 (later): lineage follow-up study, calibration flags

- **Session:** `d6fcae14-1433-4af7-9831-a319a98cd839` (continued)
- **Done:** audit of the first report (aBSREL gene `p_value` was claimed calibrated on every tree; true on T5 only, corrected); follow-up designed and fixed before running (`phylo-lineage-followup.md`), smoke-tested, run (110.6 core-hours of 120), and reported; `calibration` column, warning, and provenance entry for flagged tests; `n_merged`/`n_absent` in `_lineages.tsv`; three-kinds-of-output section and single-copy limitation in the docs.
- **Verified:** `test-phylo-lineage.R` and `test-phylo-lineage-engines.R` in the candidate image with this code, 0 failed, 0 skipped; fast gate on intel18, 243 tests, 1,030 expectations, 0 failed (6 skips: slow or fast-mode).
- **Found:** RELAX's excess reproduced and not explained by dndsR, the simulator, or starts; `p_fg_omega_gt1` 5.5-6.8% at omega 1; free-ratio excess not reproduced. HyPhy's `--starting-points 10` errors in 2.5.103; RELAX optima differ by up to 0.6 in LRT on 14 of 300 genes between identical runs (R1 failed as specified).
- **Not done:** push, pull request, CI, candidate image rebuilt from the final commit; nothing merged, tagged, or published.

## 2026-11-21: lineage-level analyses implemented and calibrated

- **Session:** `d6fcae14-1433-4af7-9831-a319a98cd839` (continued)
- **Done:** `phylo_branch_omega`, `phylo_absrel`, `phylo_relax`, `phylo_summarize` with CLI, tests, tutorial, methods, skill, outputs; PAML in the image; calibration designed in advance and run in full; report written.
- **Verified:** engine comparisons against codeml and HyPhy run directly; candidate image (`0c7c813`) container checks and 248 package tests, 0 failed, 0 skipped; tutorial commands; 300 public Brassicaceae orthogroups; calibration 27 of 27 units.
- **Found:** the engine test's `scale = 1` was ignored by testthat 3e (relative tolerance); fixed at `4ba2f1d`. The first calibration launch filled the scratch file quota (about 140 files per orthogroup); the package now keeps about 25, or none with `keep_runs = FALSE`. RELAX and the free-ratio test with missing taxa are anticonservative in simulation.
- **Not done:** push, pull request, CI on the branch, version choice; nothing merged, tagged, or published.

## 2026-11-20 (final): v0.5.0 release preparation

- **Session:** `d6fcae14-1433-4af7-9831-a319a98cd839` (continued)
- **Done:** maintainer decisions recorded; fixture notice; clone references to `v0.5.0`; CI now fails on any skipped test in the release image and runs KaKs_Calculator and the recommended BUSTED settings in it; interim-image and g03 documentation from the agent test; channel message 11; HyPhy checked for a branch-model equivalent (aBSREL baseline MG94 per-branch ω is the closest, not equivalent without a study).
- **Verified:** agent test 7 of 7 at `010ab1e` (6 of 7 at `1ea7faf`, a wording miss); `check.sh` with the tutorial; the new container checks against this tree on Slurm.
- **Not done:** push (the guard refuses agent pushes, even when asked); CI on the candidate; the MK `gaps.md` fix.

## 2026-11-20 (later): gene IDs, name checks, CI locale fix, v0.5.0 preparation

- **Session:** `d6fcae14-1433-4af7-9831-a319a98cd839` (continued)
- **Done:** genes named by GFF3 ID (`f529235`); sequence-name checks (`c2e7233`); regression locale fix after the first `regression.yaml` run failed (`92471e1`); `p_best_fits_adj` and calibration docs; v0.5.0 version and labels; release notes draft; lineage-level plan (review only); channel message 10.
- **Verified:** gate with `--rcmd-check --regression` (236 tests, 0 failed, 4 skipped; regression byte-identical); BUSTED tests with HyPhy (3 of 3); the recommended settings on the example with an independent BH check; whole genomes for the gene-ID fix (CDS, dN/dS, annotated tables byte-identical; one QC column for 84 chloroplast transcripts); CI green at `3943ffc`.
- **Found:** safeguards block allowed pushes while a version bump is uncommitted (the test clones the last commit and copies the working tree's skill), so rerun them after committing a bump. A test run on the login node was stopped at the maintainer's request.
- **Not done:** no merge, tag, or release; `test-agent.sh` not rerun.

## 2026-11-20: PR #9-#10 integration, --unique-pairs across commands, tiered testing

- **Session:** `d6fcae14-1433-4af7-9831-a319a98cd839`
- **Synced:** fast-forward to `9a5c758` (PRs #9, #10), hooks reinstalled.
- **Done:** PR #9 closes message 7, item 2. PR #10's `--unique-pairs` broke `append_annotations` on the skipped row and dropped reversed MK rows; fixed per the maintainer's choices (`d310f77`). Tiered testing with a public *A. thaliana*-*A. lyrata* fixture, regression runner, skip reasons, and CI (`6fe718d`). Channel message 9.
- **Verified (v0.4.0 image):** `check.sh --rcmd-check --regression`: tests 222, 0 failed, 0 errors, 4 skipped (HyPhy, KaKs_Calculator absent); R CMD check, tutorial, safeguards, skill ok; fixture byte-identical. Whole genomes, v0.4.0 against `9a5c758`: CDS, QC, dN/dS, annotated tables byte-identical on all three pairs. Fixture: 503 of 503 shared pairs identical to the whole-genome run; v0.4.0 reproduces it.
- **Found:** gene IDs switch between locus tag and symbol with the GFF's name uniqueness, and annotations then go missing (open question). 24 GB is too little for `check.sh` on 8 workers (OOM in `test-threads.R`); runbook says 48 GB. Host `python3` dies with "Illegal instruction" on some compute nodes; the fixture comparison runs in R in the image.
- **Not done:** `regression.yaml` not yet run on GitHub. No end-to-end run with two SNP files for one pair. No push, merge, tag, or release.
- **Cleanup (approved):** scratch file quota; see STATUS.

## 2026-11-15: R CMD check in the local gate; feat/* safeguard tests

- **Session:** `5e80dd8b-3ab2-4679-a25e-00f0d7dc3366`
- **Done:** `tools/dev/check.sh --rcmd-check`, with `tools/dev/rcmd-check.sh`, closes channel message 3's second point. Five `feat/*` pre-push cases added to `test-safeguards.sh`. Channel message 8; messages 3 and 5 marked answered.
- **Verified, after `fda5dbe`, v0.4.0 image:** `check.sh --rcmd-check` on Slurm: tests 196, 0 failed; imports, skill ok; R CMD check `Status: OK`. Safeguards and tutorial crashed there, because `git` itself dies with "Illegal instruction" on those compute nodes, so they were rerun on the login node. Tutorial OK. Safeguards: 96 passed, 1 failed, the known stale installed `pre-push` copy.
- **Then:** the maintainer pushed `9ed47af`. At their request, agents now rerun `tools/dev/setup.sh` after every sync (git-policy, `CLAUDE.md`, `AGENTS.md`, RUNBOOK), only when the hook sources match `origin/dev`. Ran it here: safeguards 97 passed, 0 failed. The release plan for the next session is at the top of `STATUS.md`.
- **Not done:** message 7, items 1 and 2, still wait on the MK side. No merge, tag, or release.

## 2026-11-14 (late): stage 2 decisions, report audit, real-genome rerun, readiness

- **Session:** `b6da01bb-9705-4b1a-b515-5e35c72e11b0` (continued)
- **Done:** the maintainer approved the stage 2 guidance; recorded in `DECISIONS.md` and carried into `docs/phylo-methods.md`, the tutorial, and the skill. The stage 2B report was audited against its tables and D1-D9: one wrong claim corrected (C3 + F1 *was* tested on realigned indels, in 2A only), plus smaller fixes. Channel message 7 reviews PR #8 for integration. v0.4.0 and `dev` rerun on real genomes. Next-release readiness written.
- **Verified:** CI green at `b65bf8b` on every workflow. Real genomes (public pair and two unpublished pairs, same image, data read-only): every CDS, QC, dN/dS, and annotated table byte-identical to v0.4.0.
- **Not done:** no release, merge, tag, or image. MK stays experimental. The maintainer reruns `tools/dev/setup.sh` for the stale `pre-push` hook copy.

## 2026-11-14: stage 2B main blocks, B2 launch, PR #7 integration, comparison-file headers

- **Session:** `b6da01bb-9705-4b1a-b515-5e35c72e11b0` (continued)
- **Done:**
  - Stage 2B: the main blocks finished (101 of 101 units, 1,194.5 core-hours). D4 failed in one of 12 audited cells, so B2 was launched with 3 starts, as prespecified and approved.
  - CI fixes and MK documentation corrections; a review of the MK code, in the channel `docs/dev/channels/mcdonald-kreitman.md`.
  - PR #7 reviewed after merge.
  - Fixed `calculate_dnds` ignoring a header row in comparison files (v0.2.0 to v0.4.0), with regression tests.
  - Ownership and the git sync rule reconciled, at the maintainer's direction. Merging, not rebasing, now applies when a commit's identity is cited.
  - `DECISIONS.md` records corrected: MK choices are the owner's, not maintainer rulings.
  - Miles Roberts added to `CITATION.cff`, `.zenodo.json`, and `Authors@R`.
- **Verified:**
  - `check.sh`: 193 tests, 0 failed, 0 skipped.
  - The new header tests fail on the old code and pass on the new.
  - CI passed at `6b3245f`.
- **Then:** B2 finished (14 of 14, 281.1 core-hours; stage 2B 1,475.5 of 1,806), and the stage 2A/2B report was written (`docs/validation/phylo-calibration-stage2b.md`). PR #8 (MK allele-aware, the MK owner's) was merged on GitHub; its records were tidied. A Windows-only test failure in the new header tests (path separators) was fixed in the tests.
- **Not done:** PR #8's MK science was not reviewed (the MK owner's area); its fail-fast test, requested in channel message 4, is in. The safeguard check reports this clone's pre-push hook copy as stale after PR #8: the maintainer reruns `tools/dev/setup.sh`. Nothing was pushed, merged to `main`, tagged, or published.

## 2026-11-14: P3 stage 2B smoke test and budget guard

- **Session:** `b6da01bb-9705-4b1a-b515-5e35c72e11b0`, resuming `83fbd26d-6d94-4764-aee9-e34ca0f4eb72`
- **Done:** a budget guard for stage 2B: each unit checks the core-hours used, the estimates of units running, and its own against a cap before it starts; the submit script checks projected cost against the 1,806 ceiling. Per-unit estimates are in `phylo-stage2b-estimates.tsv`.
- **Verified:** the smoke test at `0aa137d` finished 101 of 101 units with 0 failed, and its summary and D4 audit ran. The guard was tested on a fake root (it passes under the cap, stops and writes `BUDGET_STOP` over it, and does nothing without a cap). The tables script reproduces every committed table byte for byte. Estimates: 1,458 core-hours for the main blocks, and 120 or 361 for B2.
- **Found:** the smoke test used 386 allocated core-hours, because 48 cores were held for 4 genes per cell.
- **Failed or not done:** the main blocks are not launched, pending the maintainer's decision on how the smoke test counts against the budget. Nothing pushed, merged, tagged, or published.

## 2026-11-14: P3 calibration stage 2A and pairwise reliability flags

- **Session:** `83fbd26d-6d94-4764-aee9-e34ca0f4eb72`
- **Approved by the maintainer:** a staged calibration study with conditions (configurations, realistic corruption, filtering, starts, reproducibility, a bounded batch first), and reliability warnings for pairwise dN/dS. No default, release, tag, or merge changes without a separate approval.
- **Done:** stage 2 plan and tooling committed before any job (`b208db5`, `1b3aaf3`); stage 2A, 108 units and 21,270 fits from a frozen commit, library, and image on one processor class, reported in `docs/validation/phylo-calibration-stage2a.md`; per-pair reliability flags in a separate `_dnds_reliability.tsv` (`3da3075`), with the dN/dS table byte-identical.
- **Verified:** 0 failed fits; every configuration calibrated under correct nulls (0-1.3%); the empirical corruption profile agrees with dndsR's window counts for all 48,765 real sequences; fast check 165 tests, 0 failed.
- **Found:** BUSTED-E costs much power (A(10, 0.05): 39% to 8-11%); BUSTED-E leaves 24-28% false positives from one realistic corrupted segment, and masking flagged windows leaves 0-0.7% with no measurable power loss; double hits fix multinucleotide changes; one start understates BUSTED-E's power (39% to 53% with 3); HyPhy can report LRT 0 after a refit while its likelihoods give up to 2.1.
- **Failed or not done:** stage 2B awaits a decision, since the rules as written exceed the budget; 10-taxon parameters not defined; nothing pushed, merged, tagged, or published.

## 2026-11-13 to 14: P3 pilot, real-genome benchmark, and pairwise dN/dS methods

- **Session:** `83fbd26d-6d94-4764-aee9-e34ca0f4eb72`
- **Done:** P3 pilot (25 cells, about 5,000 BUSTED fits) and the real-genome benchmark on five public Brassicaceae genomes, reported in `docs/validation/phylo-calibration-pilot.md` with a full-study proposal. At the maintainer's request (a collaborator asked), `calculate_dnds --dnds-method` NG, LWL, LPB, MLWL, YN, MYN through KaKs_Calculator 2.0.1, validated in `estimators.md`; `verify_example` refuses other methods' results. McDonald-Kreitman tests recorded as a later to-do.
- **Verified:** pilot rates and paired comparisons; dndsR equals HyPhy run directly on 40 real orthogroups (LRT to 5e-14), and the one HyPhy failure reproduces directly; methods identical to KaKs_Calculator run directly and with no zero conflicts against PAML yn00 on 22,149 real pairs; a development image with KaKs_Calculator passed all 41 container checks; `check.sh` 162 tests, 0 failed, 0 skipped, tutorial passed.
- **Found and fixed:** a gene-tree bug wrongly excluded 918 four-taxon orthogroups (`f483869`); the pilot's AVX2 job died on `grep` under `set -e` (`cb87647`); the pilot summary never read the multi-start columns (`ed1bf93`); my first zero rule for KaKs_Calculator gave dN = 0 for identical proteins, wrong for 27 real YN pairs (`ff2f758`, before any push); KaKs_Calculator's cap of 99 is now NA; the validation script misread yn00's `nan`.
- **Found, reported:** with defaults, multinucleotide changes (95%) and misaligned windows (87.5%) give severe false positives; a single optimizer start misses optima under selection (9 of 100 calls change with best fits over 3 seeds); in real genomes, 148 of 322 BH-significant genes carry a foreground dissimilar-window flag.
- **Provenance gap:** `check.sh` reinstalled the library the pilot loaded dndsR from during the run; only unused gene-tree code changed. The full study needs a dedicated library.
- **Failed or not done:** the full calibration study awaits approval; Q7 filtering sensitivity not run; no default changed; nothing published, tagged, or merged.

## 2026-11-13: phylogenetic module, P2

- **Session:** `83fbd26d-6d94-4764-aee9-e34ca0f4eb72`
- **Decided by the maintainer after P1:** fixed seed and threads stay default, with extra starts as diagnostics only; HyPhy defaults kept, BUSTED-E experimental; strict identifiability, with mathematical and policy exclusions distinguished; gene-tree discordance reported, never excluding; version 0.4.0.9000; development image built and tested locally, nothing published.
- **Done:** `phylo_prepare` (genomes or CDS and long or OrthoFinder orthogroups to codon alignments, labelled trees, QC, member and exclusion records, provenance, manifest); `phylo_busted --prepared-dir`, `--extra-seeds`, `--gene-tree-mode`; outcome bases and stages; a six-genome synthetic example with designed failures; direct-HyPhy reference for species and gene trees; tests; docs, tutorial, and skill.
- **Verified:** a development image built from `d24d183` passed all 41 `container-checks.sh` checks; package tests 153, 0 failed, 0 skipped (Intel node), phylogenetic tests on AVX2 and AVX-512 nodes; `R CMD check` OK; all outcomes as designed; alignments equal the designed CDS and a true indel alignment; gene-tree mode equals direct HyPhy exactly; the tutorial ran in a clean `HOME`; the workflow passed in a development image built locally with HyPhy.
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
