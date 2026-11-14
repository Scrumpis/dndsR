# dN/dS estimator validation

## Question

dndsR's default is the Comeron (1995) method, through `orthologr`. It stays the default for v0.3.0. Is it correct, and how does it behave?

An earlier version of this page said it reads low, based on simulations with uniform codon usage. Those simulations were misleading, as explained below. This page replaces that statement.

## Conclusions

1. **No implementation error was found in how dndsR uses the estimator, or in the estimator's behavior relative to related methods.** The code is libsequence's `Comeron95`. Independent KaKs_Calculator methods of the same family (MLWL, LPB) err in the same direction on the same data.
2. **It compresses ratios toward the middle, and more so at high divergence.** With realistic codon usage it overestimates very low ratios and underestimates high ones. This is a property of this family of site-class methods. YN00 was closer in these simulations, which favor it because they share its assumptions.
3. **The earlier low reading was mostly a simulation limitation.** Uniform codon usage is unrealistic, and it made Comeron about 10% low at every ω.
4. **A single gene's ratio is noisy, whichever method is used.** About half of truly neutral genes exceed 1. A ratio above 1 does not by itself distinguish positive selection from neutrality.

## Method

`tools/dev/validation/comeron-bias.R` simulates long pairs (50,000 codons, one per cell, so noise is small) under a Goldman-Yang codon model, with κ of 1, 2, 4, and 8; ω of 0.1, 0.5, 1, and 2; branch lengths of 0.05, 0.2, 0.6, and 1.5 substitutions per codon; and two codon usages: uniform, and empirical from *A. thaliana* CDS. True dN and dS use the model's mutational-opportunity definition of sites, so the true ratio equals ω.

`tools/dev/validation/comeron-genes.R` repeats this at realistic lengths (150 to 1,500 codons, 100 replicates per cell, κ of 2 and 4, *A. thaliana* codon usage), to measure per-gene behavior.

Compared with PAML yn00 4.10.10 (YN00, NG86) and KaKs_Calculator 2.0.1 (LWL, MLWL, LPB). Results: [data/](data/).

## Bias, with realistic codon usage

Estimated ratio divided by ω, averaged over κ (long pairs):

| ω | Length | Comeron | YN00 | MLWL | LPB | NG86 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.2 | 1.08 | 1.04 | 1.17 | 1.21 | 0.93 |
| 0.1 | 1.5 | 1.26 | 1.03 | 1.34 | 1.47 | 0.99 |
| 0.5 | 0.2 | 1.04 | 1.04 | 1.00 | 1.11 | 0.92 |
| 0.5 | 1.5 | 1.06 | 1.05 | 1.00 | 1.15 | 0.97 |
| 1 | 0.2 | 1.00 | 1.03 | 0.92 | 1.07 | 0.91 |
| 1 | 1.5 | 0.94 | 1.01 | 0.85 | 1.01 | 0.89 |
| 2 | 0.2 | 0.93 | 0.99 | 0.86 | 0.98 | 0.87 |
| 2 | 0.6 | 0.89 | 0.97 | 0.81 | 0.95 | 0.84 |
| 2 | 1.5 | 0.79 | 0.91 | 0.71 | 0.86 | 0.78 |

The cause is mainly dS. At ω = 2, Comeron's dS was 6 to 15% high from length 0.2 up, and at ω = 0.1 it was 5 to 19% low, while dN stayed within a few percent until very high divergence.

**Uniform against realistic codon usage,** mean ratio/ω across ω at lengths up to 0.6:

| κ | Comeron, uniform | Comeron, realistic | YN00, uniform | YN00, realistic |
|---:|---:|---:|---:|---:|
| 1 | 0.87 | 1.03 | 1.02 | 1.03 |
| 2 | 0.89 | 1.04 | 1.02 | 1.05 |
| 4 | 0.91 | 1.01 | 1.02 | 1.03 |
| 8 | 0.91 | 1.01 | 1.00 | 1.03 |

## Single genes, with realistic codon usage

κ = 2, 100 replicates per cell. Share of genes with an estimated ratio above 1:

| True ω | True dS | 150 codons | 300 | 600 | 1,500 |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.027 | 19% | 15% | 7% | 0% |
| 0.5 | 0.108 | 8% | 2% | 0% | 0% |
| 0.5 | 0.32 | 0% | 0% | 0% | 0% |
| 1 | 0.017 | 49% | 48% | 50% | 51% |
| 1 | 0.067 | 60% | 57% | 53% | 53% |
| 1 | 0.20 | 49% | 50% | 49% | 40% |
| 1 | 0.50 | 42% | 33% | 24% | 23% |
| 2 | 0.009 | 67% | 84% | 94% | 98% |
| 2 | 0.038 | 93% | 96% | 99% | 100% |
| 2 | 0.114 | 100% | 99% | 100% | 100% |

So:

- **Neutral genes cross 1 about half the time** at ordinary divergence. That is expected of any roughly unbiased estimate, and using 1 as a cutoff has a high false-positive rate among neutral genes. At very high divergence fewer cross, because the estimate reads low.
- **Genes under strong purifying selection rarely cross 1,** except short genes at low divergence (up to 19%, with ω = 0.5 and dS near 0.03).
- **Most genes with ω = 2 cross 1,** so the threshold has power. These simulations use one ω for the whole gene, which real genes rarely have: most positive selection acts on few sites, and would not raise a whole-gene ratio this much.
- **Short genes at low divergence often have undefined ratios** (up to 30% with 150 codons and dS below 0.01), because dS is 0.

## Benchmark against independent estimators

`tools/dev/validation/estimator-benchmark.R` simulates pairs at 540 settings: 100 to 1,600 codons; dS of 0.005 to 3; ω of 0.1 to 3; and either one ω per gene (`uniform`) or ω varying among sites, a gamma with shape 0.5 (`sites`). Each has 200 replicates, κ = 2, and *A. thaliana* codon usage. Estimators are Comeron (dndsR's default), PAML yn00 (YN00 and NG86), PAML codeml (pairwise ML), and KaKs_Calculator (MLWL and LPB). Results: [data/estimator_benchmark_cells.tsv](data/estimator_benchmark_cells.tsv). Summary: `estimator-benchmark-summary.py`.

Median estimated ratio over the true ratio, 400 codons, `uniform`:

| True ω | Estimator | dS 0.02 | 0.1 | 0.4 | 0.8 | 1.5 | 3 |
|---|---|---:|---:|---:|---:|---:|---:|
| 0.1 | Comeron | 0.98 | 1.03 | 1.13 | 1.13 | 1.22 | 1.48 |
| | YN00 | 1.01 | 1.05 | 1.05 | 0.99 | 1.00 | 1.47 |
| | ML | 0.96 | 1.00 | 1.07 | 1.04 | 1.03 | 1.07 |
| 1 | Comeron | 1.03 | 1.01 | 0.95 | 0.88 | 0.74 | 0.73 |
| | YN00 | 1.09 | 1.07 | 1.02 | 1.01 | 1.01 | 1.36 |
| | ML | 1.05 | 1.06 | 1.04 | 1.01 | 1.05 | 0.86 |
| 3 | Comeron | 0.91 | 0.78 | 0.55 | 0.39 | 0.28 | 0.26 |
| | YN00 | 1.05 | 0.91 | 0.74 | 0.65 | 0.59 | 0.48 |
| | ML | 1.06 | 1.05 | 0.97 | 0.95 | 0.88 | 0.33 |

Comeron is close at ordinary divergence for low ω, and falls increasingly short at high ω and dS. ML is closest, and yn00 next. **When ω varies among sites, all of them fall far short at high ω:** at a mean ω of 3, 400 codons, and dS 0.1, Comeron reads 0.55, yn00 0.63, and ML 0.68 of the truth. At a mean ω of 1 and dS 0.8, they read 0.50, 0.47, and 0.50. So a better pairwise estimator would not remove this.

**Failure conditions, Comeron.**

| Measure | Result |
|---|---|
| Undefined ratios, ω = 0.6 | At dS 0.005: 70% at 100 codons, 23% at 400, 0% at 1,600. At dS 0.02: 20.5% at 100 codons, none at 400. At dS 3: 58–70% |
| Replicates above 1, true ω = 1.5, 100 / 400 / 1,600 codons | dS 0.8: 62 / 70.5 / 94.5%. dS 1.5: 41 / 9.5 / 1%. Longer genes lose more, as the estimate tightens around a biased value |
| Replicates above 1, true ω = 1, 100 / 400 / 1,600 codons | dS 0.4: 45 / 37 / 25.5%. dS 0.8: 27 / 19.5 / 4%. About half at dS 0.2 or below |

So calls above 1 stop meaning the same thing around dS 0.4 to 0.8, and dndsR now reports pairs above 0.8. This threshold comes from these simulations, and from wrong pairings above: see [orthology](orthology.md).

**A bootstrap interval is not offered.** A 95% percentile bootstrap, with 100 resamples of codons on 50 replicates per setting, gives a median coverage of the true ratio, across lengths, of:

| Regime, true ω | dS 0.02 | 0.1 | 0.2 | 0.4 | 0.8 | 1.5 |
|---|---:|---:|---:|---:|---:|---:|
| `uniform`, 0.6 | 92 | 90 | 90 | 88 | 90 | 98 |
| `uniform`, 1 | 96 | 92 | 94 | 92 | 84 | 84 |
| `uniform`, 1.5 | 94 | 92 | 90 | 82 | 22 | 34 |
| `uniform`, 3 | 96 | 76 | 36 | 2 | 2 | 0 |
| `sites`, 1 | 94 | 90 | 56 | 16 | 0 | 0 |
| `sites`, 3 | 90 | 14 | 0 | 0 | 0 | 0 |

Resampling measures noise, not bias, so coverage fails in the cases that matter most for dN/dS above 1. Under the roadmap's rule, this is not established as useful, so no feature was added.

## Exact cases on the designed example

- Identical sequences: all methods give dN = dS = 0.
- Synonymous changes only: Comeron and NG86 give dN = 0.
- Nonsynonymous changes only: NG86 and YN00 give dS = 0, while Comeron gives a small nonzero dS. Comeron counts by site class, not by codon change. This is a known kind of difference between methods.

## Other methods

`calculate_dnds --dnds-method` also offers `NG`, `LWL`, `LPB`, `MLWL`, `YN`, and `MYN`, computed by KaKs_Calculator 2.0.1 on the codon alignments orthologr makes for every method. Comeron stays the default. KaKs_Calculator's `GY` is left out: on an example pair with no synonymous change it gave dS 0.00096 and dN/dS 50. `MLPB` is left out because no independent implementation was available to check it.

**What dndsR adds to KaKs_Calculator's output.**

- **Zeros.** KaKs_Calculator writes "NA" both for 0 and for an undefined value. dndsR reports 0 only where the alignment proves it. For NG, YN, and MYN, which count differences along mutational pathways, dN = 0 needs every step of every pathway between every compared codon pair to keep the amino acid, and dS = 0 needs every step to change it. A pathway through a stop codon decides neither. For LWL, LPB, and MLWL, which classify sites by degeneracy, 0 needs identical compared codons. Those methods count a synonymous transversion at a twofold site in dN, and their corrections can be negative. Identical proteins are not enough: on *A. thaliana* and *A. lyrata*, an earlier version of this rule gave dN = 0 for 27 YN pairs where all three NG86 implementations give dN > 0 (codon pairs such as Leu TTA and CTC, one of whose pathways passes through Phe TTC). It was fixed before release.
- **Cap.** KaKs_Calculator's YN writes 99 for a distance it cannot estimate (yn00, from which its code derives, caps at 99). dndsR reports NA and warns. This happened in 2 of 22,149 real pairs, where yn00 gives dS of 3.3 and 4.4.
- **Names.** Pairs are passed to KaKs_Calculator by number, so gene names containing "-" are kept. orthologr's own KaKs_Calculator path splits names at "-".

**Method.** `tools/dev/validation/dnds-methods.R` (Slurm: `dnds-methods.sb`) runs `calculate_dnds` with every method on the same pairs. It then checks:

1. that the codon alignments are byte-identical across methods;
2. that dndsR's values equal KaKs_Calculator run directly on those alignments, wherever KaKs_Calculator gives a number;
3. agreement with PAML yn00 4.10.10 (NG86, YN00, LWL85, LPB93) and with `ng86.py`, an independent NG86 in Python, on the codons KaKs_Calculator compares (codons with a gap, ambiguity, or stop in either sequence removed).

MLWL and MYN have no independent implementation here.

**Data.** The bundled example (32 and 33 pairs), and *A. thaliana* TAIR10.1 against *A. lyrata* v.1.0 (RefSeq, as in [public-genome.md](public-genome.md)): 22,149 ortholog pairs. CDS sha256 `5366bf48...` (*A. thaliana*) and `0678d08c...` (*A. lyrata*). dndsR at `ff2f758`, in a development image built locally from `tools/Dockerfile` at `8af8ca9` (`sha256:e24b394df53597922308224fd2507956fb6f70bd146e470540b1c91179f1d17f`, not published), on an Intel Xeon Gold 6148.

**Results, real pairs.**

- Alignments were byte-identical across all 7 methods. dndsR equalled KaKs_Calculator run directly, with a largest difference of 0, for every number KaKs_Calculator gave (about 44,000 per method).
- **No zero conflicts with yn00.** Every 0 that dndsR reported for NG, YN, LWL, or LPB was within yn00's rounding of 0: 252 values each for NG and YN, and 2 each for LWL and LPB.
- **NA left as NA.** NG had 47 values, all dS of saturated pairs: yn00 also fails on 40, and gives a median dS of 3.7 for the rest. YN had 27 (dN, the pairs above) and MYN 44. LWL and LPB had 158 each, and MLWL 251. For LWL, 149 are dN of pairs where NG gives dN = 0: identical proteins with synonymous differences, where the site-class estimate may be slightly negative, so 0 cannot be proved. Users who need a ratio for such genes should use NG or YN.

Ratio of dndsR's value to yn00's (or `ng86.py`'s), medians and shares within 5%, by yn00's dS:

| Method, against | Median dS ratio | dS within 5%, dS ≤ 0.5 | dS within 5%, 0.5-1 | dS within 5%, > 1 | dN within 5%, all | Spearman, dN/dS |
|---|---:|---:|---:|---:|---:|---:|
| NG, yn00 NG86 | 1.014 | 100% | 100% | 59% | 100% | 0.9999 |
| YN, yn00 YN00 | 1.004 | 90-97% | 43% | 23% | 28-94% | 0.997 |
| LWL, yn00 LWL85 | 1.005 | 94-98% | 96% | 68% | 92-96% | 0.998 |
| LPB, yn00 LPB93 | 1.003 | 95-99% | 94% | 64% | 93-98% | 0.998 |

- **NG reads dS about 1.4% higher than yn00's NG86, and dN about 0.4% lower,** consistently. On example pairs the substitution counts are identical, but KaKs_Calculator counts fewer synonymous sites: 226.5 against 230.1 on one gene, where `ng86.py` and yn00 both exclude changes to stop codons. KaKs_Calculator's convention for those sites is not documented. Against `ng86.py`, the median dS ratio is 1.004. Above dS 1, the Jukes-Cantor correction magnifies the difference.
- **YN agrees closely at low divergence, and less above dS 0.5.** YN00 estimates κ and codon frequencies from the pair, and the two programs' iterations differ there. On example pairs, dN and dS were within 0.4% of yn00.
- **LWL and LPB agree closely, with some dN up to half of yn00's.** Their dN rank-correlates at 0.998.
- **Methods rank genes alike, not identically.** Spearman correlations of dN/dS between methods are 0.95-0.99. Methods of the same family agree best: NG with LWL 0.991, YN with MYN 0.987.
- **Very short genes:** one 36-codon pair got dS of 9e-6 from YN and 1e-6 from MYN, so dN/dS of 4,300 and 29,000. yn00 gives dS 0.022, and NG, LWL, and LPB give 0.025-0.030. dndsR reports what KaKs_Calculator gives; treat a tiny positive dS with a large ratio as unreliable, whatever the method.

**Results, example pairs.** NG and YN: zeros agree with yn00 for all pairs. YN dN 0.996-1.004 and dS 0.999-1.001 of yn00. NG dS 1.010-1.023 and dN about 0.5% low. LWL and LPB dN 0.86-1.005 and dS up to 4.5% high. `tests/testthat/test-dnds-methods.R` checks NG and YN against stored yn00 values for three designed genes, and the zero rule by hand.

**Use.** Choose one method before looking at results, use it for the whole analysis, and say which one. The simulated benchmark above shows how YN00 and the site-class methods behave against a known truth. Methods differ by a few percent per gene, which matters only near a threshold such as dN/dS = 1.

## Limits of this validation

The model has one ω per gene, no recombination, no site variation in rate, and stationary codon usage from one species. Real data can differ in all of these. Different codon usage may shift the results.

## Saturation behavior in the source

`Comeron95` sets a term to 0 if its logarithm is undefined, and sets negative components to 0. Near saturation, or with very few differences, values are therefore truncated instead of flagged. `calculate_dnds` reports how many pairs have dS above 0.8, from the benchmark above.

## Pinned orthologr

Results on the example were byte-identical between the orthologr in the published image and upstream commit `06c3098`, the commit the Dockerfile pins.
