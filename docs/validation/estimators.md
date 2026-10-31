# dN/dS estimator validation

## Question

Does `calculate_dnds` report the dN and dS it should? dndsR's default is the Comeron (1995) method through `orthologr`, which calls a compiled estimator. The defaults were kept, as established behavior.

## Method

`tools/dev/validation/compare-estimators.R` sends the same aligned pairs to:

- orthologr's Comeron estimator, as used by dndsR
- PAML `yn00` 4.10.10 (YN00 and NG86)
- KaKs_Calculator 2.0.1 (LWL, MLWL, LPB, MYN)

The pairs are the designed example, and sequences simulated under a codon substitution model with equal codon frequencies, 600 codons, and 20 replicates for each combination of κ (2, 4), ω (0.1, 0.5, 1, 2), and branch length (0.1, 0.5, 1.5 substitutions per codon). Base R code simulates the sequences and computes true dN and dS from the same model. The true ratio equals ω by construction.

Separately, `tools/dev/validation/ng86.py` is an independent Nei-Gojobori calculation with no shared code, and it agrees with PAML's NG86 on the example's exact-count cases.

## Results

Selected results at κ = 2, mean estimate versus truth:

| ω | Length | True dN | Comeron dN | True dS | Comeron dS | YN00 dS | Median ratio, Comeron | Median ratio, YN00 | Comeron ratio > 1 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.5 | 0.048 | 0.045 | 0.475 | 0.509 | 0.483 | 0.09 | 0.10 | 0% |
| 0.5 | 0.5 | 0.130 | 0.130 | 0.261 | 0.273 | 0.252 | 0.48 | 0.54 | 0% |
| 1 | 0.1 | 0.033 | 0.034 | 0.033 | 0.039 | 0.034 | 0.86 | 1.04 | 30% |
| 1 | 0.5 | 0.167 | 0.159 | 0.167 | 0.181 | 0.161 | 0.91 | 1.07 | 20% |
| 1 | 1.5 | 0.500 | 0.475 | 0.500 | 0.573 | 0.502 | 0.84 | 1.01 | 0% |
| 2 | 0.5 | 0.194 | 0.189 | 0.097 | 0.115 | 0.101 | 1.73 | 2.01 | 100% |
| 2 | 1.5 | 0.581 | 0.544 | 0.290 | 0.376 | 0.309 | 1.45 | 1.90 | 100% |

Full results across κ = 2 and 4 are reproduced by the script. Across all neutral cells (ω = 1), Comeron gave a ratio above 1 in 0 to 35% of replicates. All ω = 2 cells gave 95 to 100%, and ω ≤ 0.5 gave 0%. Twenty replicates per cell makes these proportions approximate.

On the designed example:

- Identical sequences: all methods give dN = dS = 0.
- Synonymous changes only: Comeron and NG86 give dN = 0.
- Nonsynonymous changes only: NG86 and YN00 give dS = 0, but Comeron gives dS of about 0.008. Comeron's site-class approach counts some changes differently, and the MLWL and LWL methods in KaKs_Calculator also differ from NG86 in these extreme cases.

## Conclusions

- **No error found in dndsR's use of the estimator.** Comeron estimates track truth and fall among established methods.
- **The default is biased toward low ratios in this simulation.** It overestimates dS and underestimates dN relative to the model's truth, most at high divergence. YN00 was closest to truth on these simulations, which favor it because they share its assumptions about mutation. Real genomes have codon-usage bias, and these simulations do not address that.
- **`dNdS > 1` is a weak criterion for individual genes,** as explained in [interpretation](../interpretation.md).

## Decision for the maintainer

The default stays Comeron, because changing it would alter every result. KaKs_Calculator is not in the image, so its methods (such as `--dnds-method YN`) would need it added first. Whether to change the default or add it is an open scientific decision: see `docs/dev/DECISIONS.md`.

## Pinned orthologr

Results on the example were byte-identical between the orthologr in the published image and upstream commit `06c3098`, the commit the Dockerfile now pins.
