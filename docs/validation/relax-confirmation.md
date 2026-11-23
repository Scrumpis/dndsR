# RELAX calibration: confirmation on a held-out configuration (results)

Status: **confirmation run once to the design fixed before running ([relax-bootstrap-study.md](../dev/relax-bootstrap-study.md), "Confirmation", commit `9e7eb74`).** Outcome under the preregistered rule: **neither candidate p-value is confirmed.** RELAX stays flagged, its use limited to descriptive k; no column is proposed.

## Runs

T16m (16 taxa, 4 dropped at random per gene), site ω from a gamma distribution (shape 0.5, mean 0.3, 8 classes), k applied to the test branches; RELAX through `phylo_relax()` with HyPhy 2.5.103; one calibration set per unit (50 genes × 20 null replicates simulated by HyPhy from each gene's fitted RELAX null). Units R-N16 (k = 1; 552 of 600 genes tested: the others lacked both test taxa), R-K05 (k = 0.5; 183 tested), R-K2 (k = 2; 176). Root `study-R2`, Intel nodes. Summary: the stage 1 script unchanged; tables `data/relax_confirmation_*.tsv`.

**Cost: 95.7 core-hours allocated, above the 75 ceiling** (projection 63). The three calibration sets took 80 core-hours allocated against a projected 50: each ran on 16 allocated cores, and allocated rather than busy time is counted. Recorded, not reinterpreted.

## Results against the criteria

| Criterion | Result | Outcome |
|---|---|---|
| C-size: R-N16 at 0.05, interval inside 2.5-7.5% | HyPhy `p_value` 9.2% (7.0-12.0); `p_scaled` 6.3% (4.5-8.7); `p_pool` 6.5% (4.6-8.9) | **fail** for both (upper bounds above 7.5%; neither is anticonservative by its lower bound) |
| C-power: at least size-adjusted chi-square minus 5 points | k = 0.5: reference 32.2%; `p_scaled` 39.3%, `p_pool` 40.4%. k = 2: reference 14.8%; 18.8%, 21.0% | pass |
| C-FDR: 2,000 mixtures (552 nulls + 30 + 30), BH at 0.05, mean FDP ≤ 5% | `p_scaled` and `p_pool`: FDP 0, power 0; HyPhy `p_value`: FDP 65%, power 1.7% | passes as written, but with no discoveries it is uninformative |
| C-shared: heterogeneity of per-gene bootstrap means | p = 0.76, 0.89, 0.39 | no heterogeneity |
| C-upstream: HyPhy direct, LRT within 1e-3 in ≥ 95% of genes | 552 of 552 within 0.01, largest difference 3e-14 | pass |

Direction of significant k: 98.6-100% on the correct side of 1.

## Exploratory observation (not a criterion)

The null genes' mean LRT was 1.42 (se 0.08), against 1.30 in their own calibration set. The parametric bootstrap simulates from RELAX's fitted model, which cannot reproduce misspecification: here the data came from a gamma ω distribution that RELAX's three-class model does not match. In stage 1, where the simulator's classes matched RELAX's form more closely, the two agreed (1.29 against 1.25). This suggests the scaled correction under-corrects when RELAX's model is misspecified, which real genes usually are. It is a plausible explanation, not a tested one.

## What follows

- RELAX's `p_value` remains anticonservative (9.2% here, 8.2-9.3% before), and neither correction was confirmed on a second configuration. The documented status stays: RELAX is flagged, k is descriptive, and no p-value from it should be treated as calibrated.
- The ~930 core-hour study is not proposed. A further attempt would need a correction that accounts for misspecification (for example a bootstrap from a richer model than RELAX's own), which is a new method, not a confirmation, and would need the maintainer's approval.
