# Phylogenetic module: calibration, stage 2 (plan)

Approved by the maintainer on 2026-11-14, as a staged execution of the full study proposed in [phylo-calibration-pilot.md](../validation/phylo-calibration-pilot.md). This page amends nothing in [phylo-calibration-design.md](phylo-calibration-design.md). Its questions, simulator, parameters, criteria, and change log stand, and so do the pilot's results. This page adds model configurations, realistic alignment corruption, filtering, and start counts, and runs them in two stages. It was committed before any stage 2 job ran. Changes after that point go in the change log at the end, with their reason.

No package default, published interface, release, tag, or merge changes. Recommendations go to the maintainer for a separate decision.

## Questions

| | Question |
|---|---|
| S1 | Which BUSTED configuration keeps false positives at or below nominal under correctly specified nulls, resists multinucleotide substitutions and alignment errors, and keeps power, at what cost? |
| S2 | Does filtering locally dissimilar windows remove false positives caused by alignment errors, and how much genuine signal does it remove? |
| S3 | How often do a second and third optimizer start improve the fit or change the call? |
| S4 | Which configuration should dndsR recommend for exploratory analyses, and which for conservative genome-wide scans? |

## Configurations

All use dndsR's defaults otherwise, including `--srv Yes` with 3 synonymous rate classes, 3 ω classes, seed 1, 1 HyPhy thread, and `--kill-zero-lengths Yes`. Synonymous rate settings are the same in every primary comparison.

| | Options | Role |
|---|---|---|
| C0 | none (BUSTED, `--multiple-hits None`) | dndsR's default |
| C1 | `--multiple-hits Double` | |
| C2 | `--error-sink Yes` (BUSTED-E) | |
| C3 | `--error-sink Yes --multiple-hits Double` | |
| C4 | `--multiple-hits Double+Triple` | Targeted: null, multinucleotide, and two power cells only |

This is not a full factorial. C0-C3 are the 2 × 2 of error sink and double hits. Triple hits are tested only where they could matter.

## Scenarios

Simulated with `phylo-sim.R` on the P3A parameters (`params.json`, unchanged), five taxa, foreground the At+Al clade. Scenario list, gene counts, and seeds: `tools/dev/validation/phylo-stage2-scenarios.tsv`.

- **Clean nulls:** N0 at L50, L25, and L50 with divergence ×2; N1 at L50.
- **Clean alternatives:** A(3, 0.2), A(5, 0.1), A(10, 0.05), A(10, 0.1), A(10, 0.2), at L50.
- **Misspecification:** M1 (δ = 0.2), M3 (the design's random 12-codon window), M4 (indels, realigned by `phylo_prepare`).
- **Realistic corruption, from P3A.** `phylo-windows.py profile` applies dndsR's dissimilar-window rule to all 9,911 prepared P3A alignments. It agreed with dndsR's per-sequence window counts for all 48,765 sequences. Sequences with a dissimilar segment: At 3.2%, Al 2.0%, Cr 5.6%, Es 8.3%, Br 15.4%. Foreground segments have a median of 14 codons (10th-90th percentile 10-26). 18% touch a sequence end, and 48% are internal with a gap within 3 codons, the pattern of a frameshift or a misaligned exon boundary. Each corrupted scenario applies one segment per gene, with lengths drawn from the 684 empirical foreground segments, to a clean simulation:
  - R1, internal, in At: the same taxon's codons from another simulated gene, so the composition is realistic and the sequence unrelated;
  - R2, terminal, in At: as R1, at the start or end;
  - R3, frameshift, in At: the segment read one nucleotide out of frame, as compensating indels would leave it;
  - R4, internal, in Br, a background taxon;
  - R1 on A(10, 0.1): corruption on top of real selection.
- In P3A, flagged segments are evidence of possible error, not proof. Some may be real divergence. In simulation the truth is known, so the flag's rate on clean alternative genes measures how often genuine selection is flagged.

## Filtering (validation only, not added to the package)

`phylo-windows.py filter`, the same window rule:

- **F1, mask:** codons of a sequence inside its dissimilar windows become gaps.
- **F2, mask + gappy** (`mask_gappy` in unit names): F1, then columns gapped in more than half the sequences are dropped.

Both are applied to clean and corrupted alignments, nulls and alternatives, so the false positives removed can be weighed against the power lost.

## Stage 2A: bounded, stratified intermediate batch

Unit list (scenario × configuration × filter × starts): `tools/dev/validation/phylo-stage2-units.tsv`.

| Block | Scenarios | Configurations | Fits |
|---|---|---|---:|
| Configurations | All 17 (2,970 genes) | C0, C1, C2, C3 | 11,880 |
| Triple hits | N0-L50, M1, A(10, 0.1), A(10, 0.2) | C4 | 690 |
| Starts (seeds 1, 2, 3) | N0-L50, M1, A(10, 0.05), A(10, 0.1) | C0, C3 | 4,140 |
| F1 | N0-L50, A(10, 0.1), A(10, 0.2), M3, M4, R1-R4, R1 on A(10, 0.1) | C0, C3 | 3,120 |
| F2 | N0-L50, A(10, 0.1), M4, R1 | C0, C3 | 1,440 |
| **Total** | | | **21,270** |

Gene counts: 300 per clean null, 120 per alternative, and 150 per misspecified or corrupted null. At 300 genes, an observed rate of 0.5% has an exact 95% interval of 0.08-2.4%. At 150 genes and a rate near 90%, the interval is about ±5%.

For the start block, 1, 2, and 3 starts are read from the same run: each start is an independent HyPhy fit with its own seed, kept in the raw results, so "the best of seeds 1 and 2" is exact. Combining fits follows `.phylo_start_diagnostics`: the best alternative and the best null likelihoods, never the largest LRT.

## Measures, per unit and per paired comparison

- Denominators: genes simulated, tested, excluded (by reason), and failed (by reason). Failed fits are counted and classified, never dropped.
- False-positive rate or power at α = 0.05 and 0.01, with Clopper-Pearson 95% intervals. BH within each unit at q = 0.05.
- Paired comparisons on the same genes: calls changed, an exact McNemar test, and the difference in rate.
- Convergence: the share of genes where extra starts raise the best alternative or best null log-likelihood by more than 0.1 and 0.5 units; calls changed with 2 and 3 starts; HyPhy refits after a negative LRT; negative LRTs.
- Cost: core-seconds per fit by configuration, from the wall time and threads of each unit, all on one processor class.
- Filtering: false positives removed in corrupted nulls; power kept in clean and corrupted alternatives; codons and columns removed; the share of clean alternative genes flagged.

## Rules for stage 2B, fixed now

Stage 2A is evaluated and its report committed before stage 2B starts. Stage 2B scales only what 2A leaves open:

1. **Redundant configurations.** A configuration is dropped if, in every 2A scenario, its calls at α = 0.05 differ from those of a cheaper retained configuration in at most 1% of genes (and at most 2 genes).
2. **Unsafe configurations.** A configuration with a clean null cell whose 95% lower bound exceeds 0.05 is not recommended, and is scaled only to confirm that in one cell.
3. **Starts.** 2B uses the smallest number of starts k for which going from k to 3 changes at most 1% of calls in every 2A start scenario. The share of improved fits is reported either way.
4. **Triple hits.** C4 is scaled only if it changes at least 1% of calls against C1 in some 2A cell, or lowers the M1 false-positive rate.
5. **Filtering.** A filter is scaled only if it lowers false positives in some corrupted null cell. Its power cost is always reported.

Stage 2B then takes the retained configurations to the full study's scope: 1,000 genes per null cell over the length and divergence grid; the power grid (ω+ 2, 3, 5, 10 × p+ 0.02, 0.05, 0.1, 0.2); M1 at δ 0.1, 0.2, and 0.3; M2-M5 and R1-R4; 20 genome-wide sets of 500 genes, with corruption at the empirical per-taxon rates; and a 10-taxon tree. Its sizing, re-estimated from 2A's costs, goes in the change log before launch. If 2A shows that some comparisons are settled, 2B does not spend the rest of the approved budget on them.

## Reproducibility

- One commit, recorded in `freeze.json`. Stage 2 code is the code at that commit.
- dndsR installed from that commit (`git archive`) into a dedicated library, which is then made read-only. Nothing else installs into it.
- One image, by sha256 digest (the locally built development image with HyPhy 2.5.103, `sha256:e24b394df53597922308224fd2507956fb6f70bd146e470540b1c91179f1d17f`, not published). Each job checks the digest and the library before running.
- One processor class for every fit: AMD EPYC 9654 nodes (`agg`). The CPU model is recorded per unit.
- Seeds: simulation seeds in the scenario table, HyPhy seeds 1, 2, and 3.
- Resumable: each unit writes its outputs and then a `DONE` marker holding a hash of its specification. A rerun skips completed units and redoes incomplete ones.
- Raw results kept outside the repository: simulated and filtered alignments, results tables, provenance, logs, and each unit's HyPhy folder (JSON for every start) as a tarball.
- Results root: `dndsr-runtime/p3s2` on the HPC scratch space. Only summaries and the report come into the repository.

## Change log

| Change | Reason | Effect on results |
|---|---|---|
| F2 units and folders are named `mask_gappy` | A smoke test of the tooling (4 genes per scenario, in a separate root, discarded) found that `phylo_busted` rejects "+" in a run name | Before stage 2A; none |
| One start is the reported result; two and three combine likelihoods. Genes where HyPhy's reported LRT differs from twice its final likelihood difference are counted | The smoke test found HyPhy reporting LRT 0 after a refit while its final likelihoods differed by 0.017 (LRT 0.034). dndsR's `lrt_best_fits` uses likelihoods, so it inherits this | Before stage 2A; none |
| Stage 2A finished and reported in [phylo-calibration-stage2a.md](../validation/phylo-calibration-stage2a.md): 21,270 fits, 0 failed, 635 core-hours | As planned | None |
| Stage 2B not launched. Its proposed scope departs from rules 1 and 3 and from the full study's cell sizes, and is left to the maintainer | Applying the rules as written (no configuration redundant or unsafe, 3 starts, C4 and F2 allowed) costs very roughly 12,000 core-hours, against about 2,450 remaining; and the rules' 1% paired threshold is below the start-to-start noise that 2A measured (4-17.5% of calls in power cells) | Stage 2B's scope awaits a decision |
| A revised stage 2B design replaces the "Rules for stage 2B" above: [phylo-calibration-stage2b.md](phylo-calibration-stage2b.md), proposed and awaiting approval | The maintainer revised 2B's priorities and budget after 2A; 2A showed the 1% paired threshold is below optimizer noise | Stage 2B not launched |
