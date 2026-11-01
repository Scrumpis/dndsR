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

## Exact cases on the designed example

- Identical sequences: all methods give dN = dS = 0.
- Synonymous changes only: Comeron and NG86 give dN = 0.
- Nonsynonymous changes only: NG86 and YN00 give dS = 0, while Comeron gives a small nonzero dS. Comeron counts by site class, not by codon change. This is a known kind of difference between methods.

## Limits of this validation

The model has one ω per gene, no recombination, no site variation in rate, and stationary codon usage from one species. Real data can differ in all of these. Different codon usage may shift the results.

## Saturation behavior in the source

`Comeron95` sets a term to 0 if its logarithm is undefined, and sets negative components to 0. Near saturation, or with very few differences, values are therefore truncated instead of flagged. Treat high dS with caution. This is noted for a possible warning before v1.0.0.

## Pinned orthologr

Results on the example were byte-identical between the orthologr in the published image and upstream commit `06c3098`, the commit the Dockerfile pins.
