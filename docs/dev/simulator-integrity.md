# Validation simulator: integrity record

Requested by the maintainer on 2026-11-23 (`DECISIONS.md`) after two simulator defects were found while new cells were added. It records which simulator each study used, shows that every earlier study's inputs are what the current code generates and what the designs specify, and describes the regression test that now guards the simulator. Development tooling only; no dndsR code is involved.

## The two defects, and whether any result was affected

| Defect | Introduced | Effect | Reached a study? |
|---|---|---|---|
| `spec$error %in% ...` on a cell without an error type gave a zero-length condition | `49fd39f` (homologous indels) | the study L cells (`LS*`, `LG1`, `LA1`) stopped with an error; nothing was generated wrongly | No. Study L runs from its frozen copy (`760ca03`), which predates it. Found when RELAX cells were added; fixed at `9e7eb74` with `isTRUE()` |
| `spec$mix` partially matched a new field `mixv`, so the variable-lineage mixture also selected the fixed foreground | in the working tree only | selected genes in cell BMIXv got the fixed foreground as well as their random branches | No. Seen in a smoke test before any commit; the field was renamed (`random_lineages`) in `dbe95e5`, the first commit with that cell |

The first failed loudly; the second never reached a commit. Neither changed any reported result. The checks below show this from the data, not only from the history.

## Exact simulator versions used by each study

SHA-256 (first 16 hex digits) of the simulator files in each study's frozen source tree (`<root>/src`, made by `git archive` of the commit named).

| Root (study) | Commit | `phylo-lineage-sim2.R` | `phylo-lineage-sim.R` | `phylo-sim.R` | `phylo-windows.py` |
|---|---|---|---|---|---|
| study-L (L; pilot and full) | `760ca03` | 09845dbf97431f4b | 6eca0eaee2b554ce | 077f57f92e63ee67 | 3a6fec0e49f615a2 |
| study-R (R stage 1) | `760ca03` | 09845dbf97431f4b | 6eca0eaee2b554ce | 077f57f92e63ee67 | 3a6fec0e49f615a2 |
| study-B (B pilot) | `9bea53a` | 09845dbf97431f4b | 6eca0eaee2b554ce | 077f57f92e63ee67 | 3a6fec0e49f615a2 |
| study-B2 (B core arms) | `37aed18` | 29015362119f1eef | 6eca0eaee2b554ce | 077f57f92e63ee67 | 3a6fec0e49f615a2 |
| study-B3 (B BE3h and remaining arms) | `49fd39f` | 7571b01c5e13343f | 6eca0eaee2b554ce | 409859f94f89bd7f | 3a6fec0e49f615a2 |
| study-Q (multiple testing, stage 1) | `48aa75d` | 7a9d8ca66758b537 | 6eca0eaee2b554ce | 409859f94f89bd7f | 3a6fec0e49f615a2 |
| aln-pilot | `6962b4a` | a2f69aaacfa430f6 | 6eca0eaee2b554ce | 409859f94f89bd7f | 3a6fec0e49f615a2 |
| study-R2 (RELAX confirmation) | `9e7eb74` | 0ef9ff9007bcba4e | 6eca0eaee2b554ce | 409859f94f89bd7f | 3a6fec0e49f615a2 |
| aln-confirm, busted-e, study-Q2 | `dbe95e5` | 09c568e52bdf6535 | 6eca0eaee2b554ce | 409859f94f89bd7f | 3a6fec0e49f615a2 |

`phylo-sim.R` changed once (`49fd39f`: the optional homologous-insertion mode, default unchanged).

## Check 1: every study's cells reproduce under the current simulator

`tools/dev/validation/sim-reproduce.sh` regenerates, for each cell a study used (its own seed and `params.json`; 6-60 genes), the data with that study's frozen simulator and with the current one, and compares every output file (alignments, CDS, gene and branch truth tables) byte for byte. **37 of 37 cells identical**, across all nine frozen versions, including study R's legacy cells (C0-C2) and study L's (table: [data/simulator_reproduction.tsv](../validation/data/simulator_reproduction.tsv); cells in `sim-reproduce-cells.tsv`). The current simulator, with exact field access throughout, generates what every study used.

## Check 2: the inputs each study used came from the intended distributions

`sim-check.R` reads a simulated data set's truth tables and compares them with values derived from the design documents' numbers, not from the simulator's code: every unselected gene's branch ω mean (0.22 for the three-class cells, 0.2 for study L's, 0.2 × gene effect × branch effect where those vary, the gamma(0.5) class means by numerical integration for the RELAX cells, raised to k on the test branches); selected genes raised to 1.12 on exactly the intended branches (the fixed foreground for BMIX, BMIXh and BE3hS; one or two random branches for BMIXv); selection rates against 10% (binomial); L-A1's loss rule against its cutoff, exp(−0.245 + 0.7 Φ⁻¹(2/3)); four absent taxa per gene on T16m; no gaps in cells without indels and gaps in cells with them.

`sim-stored-check.sh` ran it on the inputs stored with every unit of every study: **all units of studies L (pilot and full), B (core, BE3h, remaining arms), Q, R2, the alignment pilot, the alignment confirmation and Q2 pass every check** ([data/simulator_stored_inputs.tsv](../validation/data/simulator_stored_inputs.tsv)). (The study R2 root also carries stage 1's units file, whose unit names the confirmation reuses with other cells; those mismatched pairings were dropped from the table, and each unit was checked against the cell it actually ran.)

What this establishes: the parameters every gene was simulated with are the designed ones. That the codon simulator turns parameters into sequences correctly was established earlier and independently: stage 2's simulator check against HyPhy's tree lengths ([phylo-calibration-pilot](../validation/phylo-calibration-pilot.md)), every error-free cell's calibrated false-positive rates (study B's E0, the lineage calibration), and aBSREL on true alignments matching E0 (study B's diagnostics). Study L's full run adds a direct check: on L-S0 (one ω per site, 0.2) the pooled estimate across 1,000 genes is reported against 0.2.

## Check 3: no ambiguous field access remains

Every `spec$...` and `par$...` access in `phylo-lineage-sim2.R` and `phylo-lineage-sim.R` is now `[["..."]]`, which never partially matches. No accessed name in any frozen version is a prefix of a different field of that version (the `mix`/`mixv` pair existed only in the working tree), so no frozen study was exposed.

## Regression test (in `tools/dev/check.sh`)

`tools/dev/test-simulator.sh`, run by `check.sh` in the image on every commit:

1. generates every cell the studies use (`sim-test-cells.tsv`, 25 cells and trees, test seeds 9001-9025) with `options(warn = 2, warnPartialMatchDollar = TRUE, warnPartialMatchArgs = TRUE, warnPartialMatchAttr = TRUE)`, so any warning, including a partial match, fails it;
2. runs `sim-check.R` on each (6-10 design checks per cell; for BMIXv a check that selected genes are not given the fixed foreground);
3. compares an MD5 of each cell's output with `sim-golden.tsv` (made from the current code after check 1 passed), so any change in what a seed generates fails the test until it is explained and the golden file updated deliberately.

**Mutation test:** with the `mixv` defect put back into a scratch copy, step 1 stops with "partial match of 'mix' to 'mixv'", and with warnings off, step 2 fails "each selected gene raised on 1 or 2 branches" and "not the fixed foreground set". Both layers catch it.

## Conclusion

No earlier conclusion is affected. Every study's inputs are reproduced byte for byte by the current simulator and match their designs' parameters.
