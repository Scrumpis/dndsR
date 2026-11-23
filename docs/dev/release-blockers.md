# Release blockers, documentation items, and deferrable work (interim, 2026-11-23)

Status: **interim classification for the maintainer; the consolidated readiness matrix follows when the running studies report.** It assumes the release ships the pairwise pipeline and enrichment as stable, and BUSTED and the four lineage commands as **experimental** (the scope decided on 2026-11-21). Shipping any phylogenetic or lineage method as stable would turn several documentation items below into blockers; that is marked. No preregistered result or criterion is changed here; exploratory findings are labelled.

## Genuine blockers (must be resolved before release)

| Item | Evidence | Smallest defensible remedy |
|---|---|---|
| Fisher-type enrichment presented as evidence of selection | Complete-null false-term rates 33-70% at one pair; FDP 20-34% ([study N](../validation/enrichment-error-calibration.md)) | Either the validated correction (study running, [design](enrichment-correction-study.md)), or, if the release cannot wait, an explicit restricted domain now: outputs and docs state the tests measure call-term association, not selection, with `covariate_enrichment` (validated at `dNdS > 1`, worst cell 1.2%) as the route for selection questions. The notice and docs already say most of this; the domain statement must be made explicit in the command help and skill |
| Final release testing | 0-skip real-engine suite passed on this branch on 2026-11-23 (251 tests, 1,210 expectations; commit `7dd5a0d`); CI on the pushed branch not yet run; candidate image predates the current code | Push, CI with `DNDSR_FAIL_ON_SKIP=1`, candidate image rebuilt from the final commit with HyPhy, PAML and KaKs_Calculator, the 0-skip suite in that image |
| MK `gaps.md` correction | Channel message 11 | The MK maintainer's edit (not ours) |
| Approvals | `RELEASE.md` | Merge, tag, image and release need the maintainer's explicit go |

## Documentation only (for an experimental phylogenetic and lineage release)

| Item | Status |
|---|---|
| Masking versus BH effect | **Done** (`b09af32`): preregistered results kept; BH adaptivity separated from genuine masking loss (about 20 points, block-selected sites) |
| `p_holm`, gene `p_value`, `p_bh_lineage` families | Done (`c40487a`) |
| aBSREL under realigned indels (every output fails FDR control) | Documented; the alignment confirmation (running) decides whether to propose an optional PRANK mode. *Blocker if aBSREL were to be stable* |
| RELAX calibration | Confirmation failed (`fa90a79`): RELAX stays flagged, k descriptive. Docs must state that neither correction was confirmed. *Blocker if RELAX were to be stable* |
| BUSTED-E | Not a remedy for non-homologous segments; forensics documented; mention in BUSTED docs that the error sink is capped at 1% of sites and targets residual errors |
| Lineage summaries (`pooled_omega`) | Study L running; until it reports, docs keep the summaries descriptive and not for trait regressions. *Blocker if offered for trait models* |
| F-all (all gene-by-branch) multiple testing | Stage 2 running; until then docs say no output controls it |
| Supplied alignments: in frame is not homologous | Done (`e9537ad`) |
| Limitations matrix's "not justified" statements per method | To be carried into the docs and skill (no code) |
| PRANK seed | Diagnostic done; only relevant if a PRANK mode is proposed |

## Deferrable (does not affect a release decision)

- The evidence record (layer B) beyond the existing flags; the agentic layer and the methods paper ([proposals](interpretation-architecture.md)).
- BUSTED-E's 41 nested likelihood shortfalls (queued reruns may close them); any upstream report.
- A misspecification-aware RELAX correction (a new method).
- `covariate_enrichment` at quantile thresholds (documented limitation unless the correction study resolves it).
- Studies P (PAML grid) and E (public end-to-end): needed before calling the branch models or the lineage workflow stable; not for an experimental release. Cheaper designs to be reviewed with the matrix.
