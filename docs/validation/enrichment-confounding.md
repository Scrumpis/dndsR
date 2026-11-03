# Enrichment: gene length, dS, and an adjusted test

Short genes and genes with small dS pass dN/dS > 1 far more often. Can that alone make a term look enriched, and can a simple adjustment remove it without losing real signals?

**Answer.** Yes to both, on this evidence. dndsR's existing tests report false terms in up to 70% of simulated genomes where omega has no link to annotation. One candidate met a rule fixed before the results, and is now `covariate_enrichment`. It held on new seeds and on a held-out pair from another family. It does not replace the defaults, which are the maintainer's decision.

## Design

The labels in [enrichment.md](enrichment.md) are shuffled, so they keep length and dS only as far as the shuffle does. Here, labels come from estimator noise:

- **Simulation (`sim-genome.R`).** Each gene keeps its real codons as an ancestor, and its length, dS, and terms. Omega is drawn independently of any term, and sequences then evolve under a codon model. **dndsR's real Comeron estimator gives each gene's dN/dS, and so its label.**
  - `long`: omega drawn from long genes' observed dN/dS, whatever a gene's length.
  - `bylen`: the same, scaled by the median dN/dS of the gene's length decile.
  - `spike_F_W`: `bylen`, then a share F of the genes in 5 chosen InterPro terms (20 to 300 genes) get omega W. This is a known signal, for power and false discovery.
- **Reassignments (`enrichment-methods.R`).** `permute`, `rotate` (around the genome, which keeps neighbors), and `strata5` and `strata10` (shuffles within bins of length, dS, and annotation richness). `strata10` keeps labels close to their real genes, so real associations can remain. It is not a clean null, and is reported anyway.
- **Pairs.** *A. thaliana* against *A. lyrata* (`AtAl`, median dS 0.14, 206 genes above 1), *Capsella rubella* (`AtCr`, 0.27, 17), and *Brassica rapa* (`AtBr`, 0.44, 3). GCF_000375325.1 and GCF_000309985.3 with the same *A. thaliana* terms as [enrichment.md](enrichment.md).

A null run counts as false if any term is significant at BH 0.05. For spikes, a term is false unless it is spiked, or overlaps spiked genes beyond chance (hypergeometric p < 0.001). Tables: [data/](data/), `enrichment_methods_*.tsv`, with 95% Wilson intervals.

The simulation reproduces the real gradient. In a pilot `bylen` run, positive genes fell from 2.22% to 0.14% across length quintiles, against 3.06% to 0.09% in real data. Across `bylen` runs, the median was 134.5 positive genes against 206 real ones, so if anything this understates the problem for *A. lyrata*. **Limits:** no indels, no recombination, one omega per gene, stationary codon usage.

## Candidates, and the rule

All run on the same labels.

| Name | Test |
|---|---|
| `fisher` | Fisher's test on all ancestors. Equals `go_enrichment(algorithm = "classic")` and `term_enrichment` |
| `go_default` | The real `go_enrichment()` defaults (`weight01`), on every fifth seed |
| `pb_len`, `pb_len_ds` | Poisson-binomial tail, with probabilities from logistic regression on splines of log length (4 df), or log length and log dS (4 df each) |
| `pb_add3` | The same, 3 df each |
| `pb_floor` | `pb_add3`, with no probability below a tenth of the positive rate |
| `pb_rich` | `pb_add3`, and log term count |
| `cmh_5x5`, `cmh_10x10` | Exact stratified test across bins of length and dS |
| `filter` | Fisher's test on genes of at least 150 amino acids and dS at least 0.05 |

An interaction term was tried first, in a pilot, and dropped before these runs: with about 100 positive genes the fit separates, giving probabilities near 1e-14 and p-values near 1e-263. That is why these use additive terms and low degrees of freedom. The candidates, including the floor of a tenth, were fixed in the script before any calibration run.

**Rule, stated before the results.** Eligible: at most 5% of runs with a false term in every null family, by point estimate, across pairs. Among those, take the highest power at spikes. The effect on real term lists is reported after, and does not count.

## Results

Largest false-term rate across 36 cells (3 pairs, GO and InterPro, 6 nulls), with `go_default` in 18:

| Method | Worst rate | Where | Cells above 5% |
|---|---:|---|---:|
| `fisher` | 69.8% | AtAl GO `bylen` | 8 |
| `go_default` | 60.0% | AtAl `strata10` (43.0% on the clean null, `bylen`) | 4 |
| `pb_len` | 11.6% | AtBr GO `permute` | 1 |
| `pb_len_ds` | 36.2% | AtBr GO `permute` | 12 |
| `pb_add3` | 27.6% | AtBr GO `rotate` | 7 |
| `pb_rich` | 37.2% | AtBr GO `permute` | 8 |
| `cmh_5x5` | 9.0% | AtAl GO `bylen` | 2 |
| `filter` | 9.7% | AtAl IPR `strata5` | 1 |
| **`pb_floor`** | **2.7%** (upper limit 5.2) | AtAl IPR `strata5` | **0** |
| **`cmh_10x10`** | **2.3%** (upper limit 4.7) | AtAl IPR `strata10` | **0** |

So `pb_floor` and `cmh_10x10` are eligible. **Without the floor the logistic tests fail badly when few genes are positive**, as in *B. rapa*, where most terms carry fits near zero.

Share of spiked InterPro terms found, mean across spike settings:

| Method | AtAl | AtCr | AtBr |
|---|---:|---:|---:|
| `pb_floor` | **77.3%** | **95.6%** | **92.8%** |
| `cmh_10x10` | 75.3% | 94.0% | 91.0% |
| `fisher`, for comparison | 70.9% | 95.6% | 94.5% |

So the rule selects `pb_floor`, which is what `covariate_enrichment` implements. It is not more powerful than Fisher's everywhere: in *B. rapa* it finds slightly fewer. In *A. lyrata* it finds more, and far fewer false terms. False discovery at spikes, AtAl:

| Spike (share, omega) | `go_default`, GO | `fisher`, GO | `pb_floor`, GO | `fisher`, IPR | `pb_floor`, IPR |
|---|---:|---:|---:|---:|---:|
| 0.05, 1.5 | 43.3% | 58.7% | 0.8% | 22.8% | 0.7% |
| 0.1, 3 | 15.4% | 35.9% | 2.2% | 8.2% | 0.7% |
| 0.2, 3 | 11.6% | 13.1% | 0.9% | 4.5% | 0.6% |
| 0.5, 3 | 6.0% | 5.7% | 1.9% | 2.7% | 0.6% |

**This revises the 9.5% in [enrichment.md](enrichment.md).** That figure comes from a shuffle within bins. With labels from estimator noise, `go_enrichment`'s defaults gave a false term in 43.0% of `bylen` runs at AtAl (95% interval 33.7–52.8, 100 runs), and in 3.3% of `long` runs (0.9–11.4, 60 runs). So length is the main driver. The excess falls with divergence: at AtCr, 5.0% (1.7–13.7, 60 runs) and 0%. At AtBr, with so few positive genes, it was 0% (0–6.4, 56 runs) and 11.1% (3.9–28.1, 27 runs), so these are rough.

## Real term lists

These did not affect the choice. *A. lyrata*, query side, at BH 0.05:

| Test | GO | InterPro |
|---|---:|---:|
| `go_enrichment` defaults | 4 | not applicable |
| Fisher's test | 9 | 24 |
| `covariate_enrichment` | **0** | **16** |

All 9 Fisher GO terms drop, including extracellular region, GO:0005576. Of the 24 InterPro terms, 15 stay, mainly F-box and Kelch families; 9 drop, among them defensin-like IPR010851; and IPR013089 joins. With *C. rubella* and *B. rapa*, `covariate_enrichment` finds no terms, with 17 and 3 genes above 1. That is low power, not a result about those genes.

**One outside comparison, chosen after seeing the result, so it supports rather than validates.** Mondragón-Palomino, Stam, John-Arputharaj, and Dresselhaus (2017, *BMC Evolutionary Biology*, doi:10.1186/s12862-017-1099-4) tested the defensin-like family in *Arabidopsis* species with codeml and FUBAR, and report that "neutral evolution prevails throughout the DEFL family". That fits defensin-like not surviving adjustment. It is one family, and we did not test the families that stay against any independent evidence.

## Confirmation

The method was selected on the runs above, so its rates there may flatter it. These runs use the shipped `covariate_enrichment`, seeds from 900,000 up, and the same pairs: 100 runs per null family, 150 for `bylen`, and 50 per spike setting. Tables: `enrichment_confirm_*.tsv`, in [data/](data/).

| Pair | Method | Worst false-term rate, of 12 null cells | Where | Cells above 5% | All null runs |
|---|---|---:|---|---:|---:|
| AtAl | `covariate_enrichment` | **5.0%** (2.2–11.2) | IPR `strata5` | 0 | 1.3% of 1,300 |
| AtAl | `fisher` | 65.3% (57.4–72.5) | GO `bylen` | 7 | 24.7% |
| AtCr | `covariate_enrichment` | 2.0% (0.6–7.0) | IPR `long` | 0 | 0.5% |
| AtCr | `fisher` | 6.0% (3.2–11.0) | IPR `bylen` | 1 | 2.1% |
| AtBr | `covariate_enrichment` | 1.0% (0.2–5.4) | IPR `permute` | 0 | 0.1% |
| AtBr | `fisher` | 3.3% (1.4–7.6) | GO `bylen` | 0 | 1.3% |

**It meets the rule again, with no margin.** The worst cell equals 5.0% with an upper limit of 11.2%, and is above the selection runs' 2.7% in the same cell. Upper limits above 5% occur in 9 of 36 cells, as expected with 100 runs per cell, so these runs do not rule out a true rate somewhat above 5% in a single setting. Across the 3,900 null runs, the rate is 0.6%.

Share of spiked InterPro terms found, and false discovery on GO:

| Spike (share, omega) | Pair | Power, `covariate_enrichment` | Power, `fisher` | GO false discovery, `covariate_enrichment` | GO false discovery, `fisher` |
|---|---|---:|---:|---:|---:|
| 0.1, 1.5 | AtAl | 82.8% | 67.2% | 2.0% | 49.0% |
| 0.1, 1.5 | AtCr | 89.2% | 85.6% | 0.5% | 1.8% |
| 0.1, 1.5 | AtBr | 87.2% | 90.0% | 1.4% | 1.9% |
| 0.2, 3 | AtAl | 100% | 100% | 0.4% | 15.8% |
| 0.2, 3 | AtCr | 100% | 100% | 0.7% | 1.6% |
| 0.2, 3 | AtBr | 99.2% | 100% | 1.4% | 1.5% |

## Held-out pair

Tomato against potato (`SlSt`, GCF_036512215.1 and GCF_000226075.1, median dS 0.09), with terms from UniProt proteome UP000004994. This is another family, another query genome, and another annotation source, and it played no part in the choice. Twice as many runs per cell: 200 per null family, 300 for `bylen`, 100 per spike setting.

| Method | Worst false-term rate, of 12 null cells | Where | Cells above 5% | All 2,600 null runs |
|---|---:|---|---:|---:|
| `covariate_enrichment` | **2.5%** (1.1–5.7) | IPR `long` | 0 | 0.8% |
| `fisher` | 7.5% (4.6–12.0) | GO `long` | 2 | 2.3% |

| Spike (share, omega) | Power, `covariate_enrichment` | Power, `fisher` | GO false discovery, `covariate_enrichment` | GO false discovery, `fisher` |
|---|---:|---:|---:|---:|
| 0.05, 3 | 17.2% | 9.6% | 3.7% | 5.3% |
| 0.1, 1.5 | 61.4% | 47.2% | 0.7% | 3.3% |
| 0.1, 3 | 67.4% | 50.6% | 0.6% | 2.7% |
| 0.2, 3 | 100% | 99.2% | 0.5% | 1.7% |

So the method holds on the held-out pair. Fisher's excess is much smaller here than in AtAl, which shows the size of that excess depends on the pair. Runs here have more positive genes than the real table: a median of 114 under `bylen` in the InterPro universe, against 43 real. Small signals are hard to find whatever the method, and half of the signals tested at low spike shares are missed.

**Real terms.** Of 19,811 pairs, 96 have dN/dS above 1. `covariate_enrichment` finds 0 of 1,734 GO and 0 of 745 InterPro terms. Fisher's test finds 4 GO terms and IPR034751 (Yippee domain, 3 of 11 genes, `p_adj` 0.0061 against 0.21 adjusted). `covariate-verify.py` recalculated both tables with no differences, and tables are identical at 1 and 2 threads. This is low power as well as adjustment, so it says nothing about those terms' biology.

**A finding on input.** The first tomato run stopped in `extract_cds`, because plastid *rps12* is trans-spliced across strands. `extract_cds` now says so and points to removing organelles. Numbers here come from nuclear annotations.

## The dS filter in the checklist

The [checklist](../publication-protocol.md) repeats the test without pairs above 3 times median dS. Is that repeat itself calibrated? `DNDSR_FILTER` runs `covariate-confirm.R` with that filter on the same labels and seeds as the confirmation, at AtAl and AtBr. Tables: `enrichment_filter_*.tsv`, in [data/](data/).

| Pair | Worst false-term rate with the filter, of 12 null cells | Where | Cells above 5% | All 1,300 null runs | Without the filter |
|---:|---:|---|---:|---:|---:|
| AtAl | 5.0% (2.2–11.2) | IPR `strata5` | 0 | 1.3% | 1.3% |
| AtBr | 1.0% (0.2–5.4) | IPR `permute` | 0 | 0.1% | 0.1% |

Spiked InterPro power at share 0.1, omega 1.5: 80.0% against 82.8% at AtAl, 87.2% against 87.2% at AtBr. Of the 700 runs compared directly, 2 changed their significant terms, so the filter does not add false terms in these simulations. The simulations do not contain mispaired genes, so they cannot show that it removes the false pairs it targets.

## How stable a term list is

On AtAl, the same test on the same pairs gives different lists when the input changes slightly. None of this was used to choose the method.

| Change, `covariate_enrichment` | GO | InterPro |
|---|---:|---:|
| Baseline, terms from a GAF file and InterProScan-style GFF attributes | 0 | 16 |
| Terms from the UniProt download in [the tutorial](../tutorials/public-genomes.md), with organelles removed: 22,145 of 22,146 pairs, one dN/dS changed | 1 | 19, with all 16 kept |
| Without pairs that OMA does not support (517) | not run | 13 |
| Without pairs above 3 times median dS (792, none of them positive) | 4 | 18, with all 16 kept |

Terms that come and go have `p_adj` just under or over 0.05, and many rest on 2 to 5 positive genes. The 4 GO terms gained under the filter, for example, each hold 2 positive genes among 11 or 12, with `p_adj` 0.049, against 0.056 without the filter. The filter is calibrated, as above, yet removing 792 genes changes the fit and the number of tests, which moves terms across the line with no new positive gene. **Stable terms are the ones to take seriously. Report the term database version, and check each term under the sensitivity runs in the [checklist](../publication-protocol.md).**

## Checks

- `covariate_enrichment` gives the same terms, and p-values within rounding, as the prototype, for 3 pairs.
- `covariate-verify.py` recalculated term counts, Fisher's p-values, and BH without dndsR on the AtAl GO table of 2,978 terms and InterPro table of 1,490 terms, with no differences. It does not refit the splines. Tests use `pbinom`, `phyper`, and hand arithmetic for the tail.
- Tables are byte-identical at 1 and 2 threads, for 3 pairs.

## Independent reproduction, at v0.4.0

Three questions, on 100 runs per null scenario, at AtAl and the held-out SlSt: 1,200 runs, 2 term families, and 2 thresholds. Table: [data/covariate_reproduction_cells.tsv](data/covariate_reproduction_cells.tsv). Scripts: `covariate-repro-task.sh`, `covariate-repro-summary.py`.

| Question | Compared | Result |
|---:|---|---|
| Do the simulations regenerate? | `sim-genome.R`, from the seed, on the newer image, against the stored file | 12 of 12 byte-identical |
| Does the result hold on the newer stack? | v0.4.0 candidate on the v0.3.0 image, against the original runs on the v0.2.0 image, same labels | The same significant terms in 2,399 of 2,400 comparisons. The other is one GO run, AtAl `permute` |
| Does it hold without dndsR's code? | `covariate-refit.py`: splines, logistic fit, tail, and BH written in Python from the method's description, with no R, on the same labeled tables | The same significant terms in all 4,800 comparisons, at both thresholds |

On one run compared in full, 2,942 GO and 1,490 InterPro terms: fitted probabilities agree within 4e-14, deviance agrees, no count differs, and `p_adj` agrees exactly. This closes the gap left by `covariate-verify.py`, which does not refit the splines.

**A failure on the way, in the check.** The first refit used plain Newton steps, and diverged when positives were rare and dS had a long tail: deviance 143,681 against R's 834.9. It reported false terms in up to 100% of runs. R's fit was right, and a damped step fixed the refit. Those outputs are not evidence, and `permute` and `rotate` had agreed only because a flat fit has nothing to diverge on.

**A reproduction is not a broader claim.** These rerun the same evidence, so its limits below stay as they are.

### Quantile thresholds

`--pos-threshold q0.95` was not part of the method's selection. On the same labels, worst false-term rate across null scenarios, except `strata10`, which keeps real associations:

| Threshold | AtAl | SlSt |
|---|---:|---:|
| `dNdS > 1` | 5.0% (2.2-11.2), InterPro `strata5` | 2.0% (0.6-7.0), InterPro `permute` |
| `q0.95` | **11.0% (6.3-18.6)**, InterPro `strata5` | 5.0% (2.2-11.2), InterPro `bylen` |

And 57% at AtAl InterPro `strata10`, against 3% at `dNdS > 1`. Shuffles within bins keep real associations in part, and a threshold passing many more genes is better able to detect them, so this may not be purely false terms. It cannot be told apart here, so **`q0.95` does not meet the rule, and calibration is not established for it.** Use `dNdS > 1` for a primary list.

## Limits

- **Confirmed at 4 plant pairs in 2 families, with 2 sources of terms, at `dNdS > 1`.** This does not establish false discovery control in general, and it has not been validated for detecting real selection in real genomes. Other clades, polyploids, other annotation pipelines, and quantile thresholds are untested or not established. For those, run the simulation on your own table.
- Calibration ran on the v0.2.0 image. The v0.3.0 image gives slightly different pairs: [stack comparison](stack-comparison.md).
- Only 100 to 300 runs per cell, so the interval on any one rate is wide.
- It treats genes as independent, and the fitted probabilities as known. Tandem families share history, and `rotate` tests that only in part.
- The floor is one setting, with no sensitivity analysis. It is why the test stays valid with few positive genes, and it can make tests conservative for terms of genes the model rules out.
- It does not adjust for annotation completeness.
- Power falls fast with divergence, whatever the method, since few genes then pass 1.
- A term here remains a lead. This is about false associations with term labels, not evidence of selection.

## Rebuild

```
Rscript tools/dev/validation/sim-genome.R <dnds_annot.tsv> <query CDS FASTA> <scenario> <seed> <sim.tsv.gz>
Rscript tools/dev/validation/enrichment-methods.R <dnds_annot.tsv> <go_ancestors.tsv> <sim.tsv.gz or permute:S ...> <terms_S.tsv>
python3 tools/dev/validation/enrichment-methods-summary.py <runs dir> <dnds_annot.tsv> <go_ancestors.tsv> <prefix>
Rscript tools/dev/validation/covariate-confirm.R <dnds_annot.tsv> <sim.tsv.gz or permute:S ...> <terms_S.tsv>
python3 tools/dev/validation/covariate-verify.py <dnds_annot.tsv> <table.tsv> <go_ancestors.tsv> q
```

Image `dndsr_v0.2.0_sha256-726657e9`, library at v0.3.0's code, seeds in the scripts' manifests: 100,000 and up for selection, 900,000 and up for confirmation. Input hashes for AtCr and AtBr: Crubella `947fda8f49900a8a` (FASTA), `322af545d18743bf` (GFF); Brapa `afd6ca07cea49c74` (FASTA), `1520fb8661ae3dbd` (GFF). For SlSt: tomato `2f7700e8342b2505` (FASTA), `98bedc4bfcb95c9d` (nuclear GFF), potato `578bbfb7e0362d8e` (FASTA), `21b1814de8924c0f` (nuclear GFF), UniProt table `2ffc67471b743fb1`.
