# Testing strategy

Three tiers. Each change runs the cheapest tier that can detect what it could break, and nothing heavier. Decided by the project maintainer, 2026-11-20 (`DECISIONS.md`).

| Tier | What | Time | Runs |
|---|---|---|---|
| 1. Unit and integration tests | `tests/testthat/` (222 tests): functions, the CLI, the bundled synthetic example end to end, McDonald-Kreitman, phylogenetic commands, estimators | 3 min on 8 CPUs (`--fast`: under 1 min) | Every commit (`tools/dev/check.sh`), every push and PR in CI |
| 2. Real-genome fixture regression | `tools/dev/regression.sh`: the full pipeline on 542 + 549 public Arabidopsis genes, compared byte for byte with a committed baseline, in the pinned image | about 2 min on 8 CPUs | Locally with `check.sh --regression` for changes to `R/`; in CI (`regression.yaml`) on PRs and `dev` pushes that touch `R/`, `DESCRIPTION`, `NAMESPACE`, or `tools/Dockerfile`; before every release |
| 3. Full-genome validation | `tools/dev/validation/` on whole genomes (the public pair and private pairs), results outside Git | about 40 min on 32 CPUs per set, through Slurm | Only for a release, a change to a core algorithm (CDS extraction, pairing, estimation, annotation), a substantial dependency update, or a scaling concern. Never automatically, and never for a change to unrelated functionality |

## Choosing the tier

- Docs, skill, or tests only: tier 1 is enough. CI skips R CMD check for docs-only pushes.
- R code that can reach pipeline output, such as a new option with its default unchanged: tiers 1 and 2. If tier 2 passes byte for byte, a whole-genome rerun adds nothing for that change.
- A change meant to alter results (an approved fix, a new default): tiers 1 and 2, with the tier 2 baseline updated and the before and after values recorded (`.claude/rules/scientific-integrity.md`), and tier 3 before the release that carries it.
- A new image or dependency versions: tier 2 in that image, which CI's `container.yaml` runs as an informative step on every image build, and tier 3 before release.

## Reading tier 2

`regression: PASS` means every output except logs and provenance is byte-identical to the baseline. On `FAIL`, the script prints the files that differ and, for the dN/dS table, pairs lost or new and the largest change per column; the work folder is kept for inspection.

The baseline is pinned to the image in `.github/workflows/regression.yaml` (the v0.4.0 image). Compare two versions in the same image; a different image is a dependency comparison, not a code one.

**Changing the baseline** (`tools/dev/regression.sh --update`) is a change in recorded output. Do it only for an approved change in results, in the commit that makes the change, with the reason and the before and after summary in the commit message and `NEWS.md`.

## Skipped tests

`tools/dev/run-package-tests.R` lists every skipped test with its reason (`SKIP <file>: <test>: <reason>`), and `check.sh` prints the list. In the v0.4.0 image, 4 are skipped: 3 BUSTED tests (HyPhy is not in that image) and the PAML comparison of other estimators (KaKs_Calculator is not in it). CI's `container.yaml` runs the same tests in a freshly built image, which has both, so they run there. Report skips with their reasons; never call a run with skips "none skipped".

## Coverage, and what each tier reuses

| Area | Evidence |
|---|---|
| CDS extraction, isoform choice, exclusions | `test-example-pipeline.R`, `test-cds-quality.R`, `test-select-isoform.R`, `test-gff.R` against designed truth; tier 2 on real gene structures; tier 3 against gffread |
| Comeron dN/dS | Example against independent NG86 (`ng86.py`); tier 3 against PAML and KaKs_Calculator (`docs/validation/`) |
| Other estimators (NG, LWL, LPB, MLWL, YN, MYN) | `test-dnds-methods.R`: method checks, zero cases, and parsing everywhere; agreement with PAML yn00 only where KaKs_Calculator is installed (CI image). Not in tier 2 |
| Reliability flags | `test-dnds-reliability.R` (hand counts); tier 2 |
| Comparison files, headers, MK columns | `test-comparison-file.R` |
| Reverse comparisons and `unique_pairs` | `test-unique-pairs.R` (the rule, including McDonald-Kreitman inputs that differ between rows); `test-unique-pairs-pipeline.R` (extract, calculate, annotate, and summarise end to end, the kept pair equal to a run of that row alone, and a missing output reported without the flag) |
| Annotation | `test-example-pipeline.R`; tier 2 on real RefSeq attributes |
| Enrichment, contrasts, summaries | `test-enrichment-*.R`, `test-covariate-enrichment.R`, `test-statistics.R`, `test-example-contrast.R` |
| CLI | `test-cli.R`; `tools/dev/test-tutorial.sh`; tier 2 runs through the CLI |
| McDonald-Kreitman (owned by the MK maintainer) | `test-mk.R` on designed sequences, with hand-derived counts and statistics |
| Phylogenetic module (BUSTED) | `test-phylo-*.R`; direct HyPhy runs where HyPhy is installed |

## Known gaps

- No public population SNP data is in any tier, so McDonald-Kreitman is tested on designed examples only. Adding real data is the MK maintainer's call.
- A comparison file with two SNP files for one genome pair is tested through the skip rule, not by producing both MK tables end to end.
- Tier 2 runs Comeron only; other estimators need KaKs_Calculator, which the pinned image lacks.
- Tier 2 cannot catch a problem that appears only at genome scale (memory, time, or pairing among many paralogs). Those are tier 3 questions.
