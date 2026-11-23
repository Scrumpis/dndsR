# BUSTED-E under a corrupted foreground segment: forensic assessment (exploratory)

Status: **exploratory analysis of existing results plus small direct HyPhy reruns, requested by the maintainer on 2026-11-23; not confirmatory.** It explains the BUSTED-E pilot's failed gate 1 ([design](../dev/busted-e-pilot.md)); the pilot's preregistered outcome stands unchanged. Findings and open questions are separated. Part of the checking was done independently by a second agent (OpenAI Codex CLI 0.171.0, read-only sandbox), which recomputed the key table from the raw files; its transcript is kept outside the repository with the run folder.

## The result being explained

300 null genes, 8 taxa, foreground t1, t2 and their stem; one foreground sequence (t1) per gene carries a contiguous segment from a different gene (non-homologous), with lengths drawn from stage 2B's empirical flagged windows. At 0.05: BUSTED 95.0%, BUSTED-E 32.3% (27.1-37.9), BUSTED on masked alignments 0.7%; on error-free genes BUSTED-E 0% (0-1.2).

## Findings

**1. Not a dndsR integration error in the settings passed (verified).** The HyPhy JSON of each gene records `error-sink: 1`, `multiple-hit: Double`, `srv: Yes` and HyPhy 2.5.103; the foreground tree is `((t1{FG},t2{FG})Node2{FG},…)`; dndsR's provenance records `error_sink = "Yes"`, three starting points, seed 1. The fitted test-branch distribution has four classes, the first being the error sink (ω 100). Codex checked these from the archived JSON and provenance independently. (Whether HyPhy run directly gives the same statistics on these inputs is being checked: *open, below*.)

**2. The error sink's design limits (verified).** In HyPhy 2.5.103's `BUSTED.bf` (lines 207-223) the error class has ω ≥ 100 (`--error-sink-bound`) and a maximum weight of 0.01 (`--error-sink-weight`). The method's authors describe it as a heuristic for *residual* alignment errors, expected to miss many real errors and to flag some real variation, with the error weight well below 1% on average in their data ([Minus the Error, eLife reviewed preprint](https://elifesciences.org/reviewed-preprints/106921); [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.11.13.620707v2.full.pdf)).

**3. False positives rise steeply with segment length (verified independently).** Segment lengths are 10-121 codons (median 14) in 300-codon genes:

| Segment length (codons) | Genes | BUSTED-E positive | BUSTED positive |
|---|---:|---:|---:|
| ≤ 10 | 37 | 5% | 89% |
| 11-20 | 204 | 21% | 95% |
| 21-40 | 52 | 88% | 98% |
| > 40 | 7 | 100% | 100% |

Codex reproduced every cell exactly. In BUSTED-E's positives the error class mostly gets no weight (median 0; 17 of 97 at the 1% cap), and the ordinary ω > 1 class takes a share close to the segment's share of the gene (median 6.6% against 7.3%).

**4. BUSTED-E's alternative fit is often not at its maximum (new, found by the independent review).** BUSTED-E contains BUSTED as a limit (error weight towards 0), so its unconstrained log-likelihood should not be lower than BUSTED's on the same data. It is lower by more than 0.01 in **117 of 300 genes** (by more than 1 in 22), against 2 of 247 for the constrained fit. Among the 188 genes that BUSTED calls and BUSTED-E does not, 89 have this deficit. So part of BUSTED-E's lower false-positive rate here may be **optimisation that stops short of the alternative's maximum**, not error absorption; a fully optimised BUSTED-E could give more false positives than 32%, not fewer. (Example: `BE0_T8b_0027`, BUSTED −2766.29 with a 2.6% class at ω 76; BUSTED-E −2767.95 with no class above 1 and LRT 0.)

**5. Masking did not suppress true signal in study B's mixtures; the reported "power cost" is BH's adaptivity (reanalysis of existing results).** In the BUSTED mixture (116 selected genes): clean selected genes have identical LRTs masked and unmasked, and the same 36 of 72 have raw p ≤ 0.05 in both; selected genes with background-only corruption gain under masking (24 to 27). Yet BH-adjusted detection falls from 30 to 16, because the unmasked run's 32 invalid discoveries loosen BH's threshold for everyone. aBSREL's gene `p_value` behaves the same way (raw: 110 against 107 of 116; BH-adjusted 107 against 82, with 494 unmasked discoveries, 78% false). Study B's non-inferiority comparison of adjusted power therefore credited the unmasked arm with power that exists only because of invalid discoveries. The preregistered classifications stand as recorded; the documentation's wording ("masking costs 12 and 22 points of power") describes BH family composition, not lost signal, and should be corrected (proposal to the maintainer). Masking can still remove genuine signal when the selected sites form one block (20 points in stage 2B).

## Which explanation fits (maintainer's categories)

| Explanation | Assessment |
|---|---|
| dndsR implementation or integration error | Settings, foreground and parsing verified; direct-HyPhy equivalence pending (below). No evidence for it so far. |
| Wrong HyPhy configuration, model, foreground, interpretation | Configuration matches the documented BUSTED-E; foreground is the intended clade; p-values read from `test results`. No evidence. |
| Optimisation or convergence | **Evidence for it** (finding 4): the alternative fit falls short of BUSTED's in 39% of genes. Its direction would raise, not lower, the false-positive rate once fixed. |
| Failure under conditions outside the method's assumptions | **Largely this.** The corruption is a whole non-homologous segment of 10-121 codons in one sequence, which the authors' design (residual error, a class capped at 1%) does not claim to absorb; false positives are low for the shortest segments and near total for the longest. |
| Limitation within the intended scope | Possible for the shortest segments (5% at ≤ 10 codons, near the nominal rate); not established. |
| Unrealistically severe test | No: segment lengths come from windows flagged in real public alignments (stage 2B), so such errors occur in real data; but they are a different error class from the residual misalignment BUSTED-E targets. |

## Open questions (being checked, or proposed)

1. **Direct HyPhy equivalence and optimisation** (running, exploratory, about 8 core-hours): the 20 most significant BUSTED-E false positives and 10 nonsignificant controls rerun in HyPhy 2.5.103 directly on dndsR's files with (a) dndsR's exact settings and seed, (b) seed 2, (c) error-sink weight cap 5%, (d) 20%. (a) tests the wrapper; (b) the optimiser; (c-d) whether the cap is what forces the segment into the ω > 1 class.
2. **Optimisation of BUSTED-E in general** (proposal): how often its alternative fit falls short of BUSTED's on clean and selected data, and whether more starts or starting from BUSTED's fit removes it; affects BUSTED-E's power and size everywhere, not only here.
3. **Masking and power on clean, block-selected genes** (existing stage 2B evidence; a confirmation would need new runs).

## What this means for dndsR now

- BUSTED-E is not a substitute for masking or for review of flagged alignments when the error is a non-homologous segment; dndsR's documentation should say so with these numbers.
- BUSTED-E's optimisation behaviour (finding 4) is worth reporting upstream once confirmed by the direct reruns, with a minimal reproducible example.
- The documented "power cost of masking" should be reworded (finding 5).
