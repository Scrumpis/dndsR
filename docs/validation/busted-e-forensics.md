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

**4. Likelihood ordering of BUSTED-E against BUSTED: mostly not comparable; 41 candidate optimiser shortfalls (corrected 2026-11-23).** The first version of this report said BUSTED-E contains BUSTED as a special case, so its alternative fit could not be lower, and read 117 of 300 genes with a lower BUSTED-E fit as optimisation failures. **That premise is wrong for this implementation**, as the maintainer pointed out. In HyPhy 2.5.103's `BUSTED.bf`, with the error sink on, the positive-selection class of both test and background branches is limited to ω in [1, 100] (`range_1_100`, lines 205, 217, 532-533), the error class to ω ≥ 100 with weight in [0, 0.01] (`range_high`, `range_small_fraction`); without it, only the test branches' positive class is constrained (ω ≥ 1, line 537) and the background's third class keeps the model's default range. Because the error weight can be 0, BUSTED-E nests BUSTED exactly when BUSTED's estimates have test ω₃ ≤ 100 and background ω₃ in [1, 100] (classes with weight below 1e-6 treated as representable). Classification of the 117 genes whose BUSTED-E alternative log-likelihood is lower by more than 0.01 (same data, tree, SRV and multiple-hit settings; `Unconstrained model` log-likelihood in both):

| | Genes |
|---|---:|
| Not nested: BUSTED's estimate outside BUSTED-E's space (test ω₃ > 100: 38; background ω₃ outside [1, 100]: 21; both: 17) | 76 |
| Nested, BUSTED-E lower (median deficit 0.21; above 1 in 6) | **41** |

Only the 41 are candidate optimisation shortfalls; their deficits are mostly small (an alternative fit 0.21 short lowers the LRT by about 0.4), and whether more starts or another seed closes them is being checked with the queued reruns. The earlier conclusion that a fully optimised BUSTED-E "could give more false positives" is withdrawn as unsupported: across non-nested cases the two models are different hypotheses, and in BUSTED's fits the background third class has ω below 1 in 216 of 300 genes, which BUSTED-E cannot represent. Table: [data/busted_e_nesting.tsv](data/busted_e_nesting.tsv).

**5. Masking did not suppress true signal in study B's mixtures; the reported "power cost" is BH's adaptivity (reanalysis of existing results).** In the BUSTED mixture (116 selected genes): clean selected genes have identical LRTs masked and unmasked, and the same 36 of 72 have raw p ≤ 0.05 in both; selected genes with background-only corruption gain under masking (24 to 27). Yet BH-adjusted detection falls from 30 to 16, because the unmasked run's 32 invalid discoveries loosen BH's threshold for everyone. aBSREL's gene `p_value` behaves the same way (raw: 110 against 107 of 116; BH-adjusted 107 against 82, with 494 unmasked discoveries, 78% false). Study B's non-inferiority comparison of adjusted power therefore credited the unmasked arm with power that exists only because of invalid discoveries. The preregistered classifications stand as recorded; the documentation's wording ("masking costs 12 and 22 points of power") describes BH family composition, not lost signal, and should be corrected (proposal to the maintainer). Masking can still remove genuine signal when the selected sites form one block (20 points in stage 2B).

## Which explanation fits (maintainer's categories)

| Explanation | Assessment |
|---|---|
| dndsR implementation or integration error | Settings, foreground and parsing verified; direct-HyPhy equivalence pending (below). No evidence for it so far. |
| Wrong HyPhy configuration, model, foreground, interpretation | Configuration matches the documented BUSTED-E; foreground is the intended clade; p-values read from `test results`. No evidence. |
| Optimisation or convergence | **Not established.** Of 117 lower BUSTED-E fits, 76 are not comparable (different parameter spaces); 41 nested cases are candidate shortfalls, mostly small, being checked with a second seed (finding 4). |
| Failure under conditions outside the method's assumptions | **Largely this.** The corruption is a whole non-homologous segment of 10-121 codons in one sequence, which the authors' design (residual error, a class capped at 1%) does not claim to absorb; false positives are low for the shortest segments and near total for the longest. |
| Limitation within the intended scope | Possible for the shortest segments (5% at ≤ 10 codons, near the nominal rate); not established. |
| Unrealistically severe test | No: segment lengths come from windows flagged in real public alignments (stage 2B), so such errors occur in real data; but they are a different error class from the residual misalignment BUSTED-E targets. |

## Open questions (being checked, or proposed)

1. **Direct HyPhy equivalence and optimisation** (running, exploratory, about 8 core-hours): the 20 most significant BUSTED-E false positives and 10 nonsignificant controls rerun in HyPhy 2.5.103 directly on dndsR's files with (a) dndsR's exact settings and seed, (b) seed 2, (c) error-sink weight cap 5%, (d) 20%. (a) tests the wrapper; (b) the optimiser; (c-d) whether the cap is what forces the segment into the ω > 1 class. The 5% and 20% caps are sensitivity analyses only, not candidate settings.
2. **The 41 nested shortfalls** (queued reruns): whether seed 2 reaches or exceeds BUSTED's likelihood; no further study is proposed unless they persist.
3. **Masking and power on clean, block-selected genes** (existing stage 2B evidence; a confirmation would need new runs).

## What this means for dndsR now

- BUSTED-E is not a substitute for masking or for review of flagged alignments when the error is a non-homologous segment; dndsR's documentation should say so with these numbers.
- Nothing is reported upstream unless the nested shortfalls persist in the reruns.
- The documented "power cost of masking" should be reworded (finding 5).
