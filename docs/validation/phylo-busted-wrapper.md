# phylo_busted: wrapper validation

## Question

Phase 1 evidence, for `phylo_busted` on supplied alignments. Phase 2 (preparation from genomes, gene-tree mode, end to end) is in [phylo-prepare.md](phylo-prepare.md). Since phase 2, outcomes say `excluded` where phase 1 said `skipped`, and the example's check function also covers the genome example.

Does `phylo_busted` give HyPhy the intended inputs, and return exactly what HyPhy computes? This is **wrapper correctness**. Whether BUSTED, as run by this pipeline, controls false positives and detects selection on realistic data is **scientific calibration**: a separate question for phase 3, not answered here. Nothing on this page shows that a significant result means positive selection.

## Software

HyPhy 2.5.103 from Bioconda, `linux-64/hyphy-2.5.103-h526e2cb_0.conda`, sha256 `57993e1477c51b574a64f7cf7391ed3ab51696fd1d41482a4b58f6dd0a1bd548`, run in and outside the dndsR v0.4.0 image (`sha256:f1b42d9d338b406ec23daca1bb146ed0329db29d21a24a6871158996aa0572da`, R 4.6.2, ape 5.8.1), on a Linux amd64 node with Singularity CE 4.1.2. Also HyPhy 2.5.104 (`hyphy-2.5.104-h526e2cb_0.conda`, sha256 `034d5dd0de916f58c89e108c2137e66020699190247ff61577b729ec7ec770ee`). Code: `dev` after `6e21e46`, in the commit that adds this page.

## Inputs

The bundled example, `inst/extdata/phylo_example`, made by `tools/dev/validation/make-phylo-example.R` (base R, seed 20261113; rerunning it reproduces the files byte for byte). Seven alignments of 300 codons on an 8-taxon tree, simulated under a BUSTED-form model with ω 0.05, 0.5, and 1 (55%, 25%, 20% of sites), κ = 2, equal codon frequencies; in "episodic" orthogroups the third class has ω = 6 on the foreground (the clade of spA and spB, with its stem). Three orthogroups are designed failures. sha256:

| File | sha256 |
|---|---|
| `alignments/OG0001.fasta` | `0cd85a1e83c024d38d2e7678beccfed756c652c8a4d56b019b5aa32a1ae31ebc` |
| `alignments/OG0002.fasta` | `c1de60d97a330dd34f6dcc1b10125e419b30057b703d2143575a0fb8c8c9e142` |
| `alignments/OG0003.fasta` | `2c4b2a08b9a35e0f6fcc63ecb96e19099ae8dd186ce6915c08da6df1edef2ac3` |
| `alignments/OG0004.fasta` | `a8ea20290acbfd524d80a636c1e57a42a896c6657edb4987d3398dbeda3e4ca1` |
| `alignments/OG0005.fasta` | `a74277fb1b3a99d9c924dceb8226f0129f9e727d146afd9408447ddb455c93ff` |
| `alignments/OG0006.fasta` | `30ec1ac21cfd37e4bce89fd8fcc8810515772c9ac5db126b3290e2b50fe33b72` |
| `alignments/OG0007.fasta` | `7d157a74563386004ed8d65fe2dc053aaa92bffa0e11b74fd46bdfadb32d1b17` |
| `species_tree.nwk` | `c27e0de3959b52bfca3ebd760f6ef7f280131c97a955e9a4b486d5278e4b336e` |
| `foreground.tsv` | `92033150c351c5b4277b4576109b89ab7e0023b53e06ca3d039787396e927769` |

## 1. Agreement with HyPhy run directly

`tools/dev/validation/phylo-direct-hyphy.sh` runs HyPhy without any dndsR code: inputs prepared with `awk`, trees written by hand from the design (internal nodes unnamed, unlike dndsR's `Node<k>`), JSON read with Python. Its output is stored as `inst/extdata/phylo_example/truth/direct_hyphy.tsv`. `phylo_busted` was then run on the example, through the R function and through the CLI.

| Orthogroup | Design | LRT | p | log L, unconstrained | Labelled branches | dndsR against direct |
|---|---|---:|---:|---:|---:|---|
| OG0001 | episodic, terminal stops | 7.307789067929662 | 0.01294505114975508 | -2548.822260013832 | 3 | identical |
| OG0002 | null | 0 | 0.5 | -2168.442215707353 | 3 | identical |
| OG0003 | episodic, spH missing | 1.865712341835206 | 0.1967142026519154 | -2407.262274266562 | 3 | identical |
| OG0004 | null, spB missing | 0.3927652384500107 | 0.4108488885736624 | -2095.702913445463 | 1 | identical |

All four agree to every printed digit (largest difference 0), so the wrapper's input preparation, labels, settings, seed, and thread count match the hand-built direct runs. `verify_phylo_example` repeats this check for users, with a tolerance of 1e-6.

OG0003 was simulated with episodic selection and is not significant. This is one simulated gene and says nothing about power.

## 2. HyPhy behavior that the design depends on

Each was confirmed by running HyPhy 2.5.103 directly:

| Behavior | Evidence | Consequence |
|---|---|---|
| `{FG}` labels and `--branches FG` select exactly the labelled branches, internal ones included | JSON `tested` field | dndsR compares `tested` with its labels on every run |
| A rooted tree is unrooted: the two root branches become one | Input tree echoed in the JSON lost the root node; rooted and unrooted versions gave identical results | Root-adjacent stems are not identifiable |
| Given branch lengths are used as starting values | At a fixed seed, LRT 8.647 with no lengths and 7.224 with all lengths 1 | dndsR passes topology only |
| Any stop codon, or a length not divisible by 3, makes HyPhy exit 1, writing `errors.log` in its working folder | Run on an alignment with an internal TGA | Runs are isolated per orthogroup; stops are handled before HyPhy |
| An unknown branch set exits 1 with a message | `--branches XYZ` | |
| p = 0.5 × (1 − CChi2(LRT, 2)); when there is no evidence, LRT = 0, p = 0.5, and no constrained model is fitted | `BUSTED.bf` lines 1429-1431; OG0002 | Documented; parser handles the missing fit |
| Without `RANDOM_SEED`, repeated identical runs differ | Two runs: LRT 8.171 and 8.180 | dndsR always sets the seed |

## 3. Reproducibility of model fitting

`tools/dev/validation/phylo-seed-spread.sb` (Slurm job 21030235) ran the direct script on the example at 20 seeds with 1 thread, a repeat of seed 1, seed 1 with 2 threads, and seed 1 on HyPhy 2.5.104. Every run: [data/phylo_seed_spread_runs.tsv](data/phylo_seed_spread_runs.tsv).

| Orthogroup | LRT, min | median | max | p, min | max | Spread of best log L, unconstrained | Seeds with p < 0.05 |
|---|---:|---:|---:|---:|---:|---:|---:|
| OG0001 | 7.007 | 8.150 | 8.787 | 0.0062 | 0.0150 | 0.590 | 20 of 20 |
| OG0002 | 0 | 0 | 0 | 0.5 | 0.5 | 0.014 | 0 |
| OG0003 | 1.527 | 1.799 | 1.870 | 0.196 | 0.233 | 0.128 | 0 |
| OG0004 | 0.182 | 0.328 | 0.493 | 0.391 | 0.457 | 0.212 | 0 |

| Comparison, at seed 1 | Result |
|---|---|
| Same seed and threads, repeated | Bit-identical, all 4 |
| 1 thread against 2 | Differs, for example OG0001 LRT 7.308 against 7.359; OG0002 identical (0) |
| HyPhy 2.5.103 against 2.5.104 | Bit-identical, all 4 |

So results are deterministic for a given build, seed, and thread count, and not otherwise. The optimizer does not always reach the same maximum: the best unconstrained log-likelihood found varied by up to 0.59 across seeds. On this example no call at 0.05 changed, but a gene with p near a threshold could change call with the seed. How to handle this (report, several seeds, or more starting points) is open for phase 3.

## 4. Automated tests

| File | What | Independent expected values |
|---|---|---|
| `tests/testthat/test-phylo-tree.R` | Scopes, pruning, every identifiability rule, root, gene trees, invalid input | Worked out by hand from the topology |
| `tests/testthat/test-phylo-alignment.R` | Each exclusion, masking, metrics, flags, FASTA parsing | Hand-calculated p-distances and fractions |
| `tests/testthat/test-phylo-hyphy.R` | Parsing real HyPhy JSON, settings, arguments, failures with stand-in executables (error, truncated JSON, label mismatch), recorded seed and settings, missing HyPhy, CLI | Values read from the JSON with Python; HyPhy's p-value formula checked with `pchisq` |
| `tests/testthat/test-phylo-busted.R` | Live HyPhy on the example against `direct_hyphy.tsv`; `verify_phylo_example`; reuse of finished runs; rerun on a changed seed | Direct HyPhy runs (section 1). Skipped without HyPhy, and in fast mode |

The JSON fixtures in `tests/testthat/fixtures/hyphy/` are HyPhy 2.5.103 output with the input file path shortened to its name; nothing else was changed.

## Limits

- One synthetic example, 8 taxa, 300 codons, equal codon frequencies. Real alignments are larger and messier.
- Agreement with direct HyPhy shows the wrapper is faithful. It does not show that BUSTED's p-values are calibrated for these data, or that alignment QC thresholds are appropriate.
- Bit-identity across machines holds only within a processor class: see [phylo-prepare.md](phylo-prepare.md#processors). OG0001 differs by 2.5e-5 in its likelihood ratio between AVX-512 and AVX2 machines; the example check now runs HyPhy directly on the same machine.
- Gene-tree mode is tested by unit tests of labelling, not yet by a live HyPhy comparison.
