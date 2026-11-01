# Gene-level `dnds_state_contrast`: problems and a possible design

Status: disabled in v0.3.0, and stops with an error. No replacement has been implemented. A redesign needs the maintainer's scientific decision. Pairwise mode is separate, and its odds ratio matches a hand calculation in `tests/testthat/test-statistics.R`.

## Confirmed by running (commit `a1ef194`)

| Input | Expected | Actual |
|---|---|---|
| One table, 200 pairs, 30% with dN/dS above 1, `min_pairs = 1`, subgenomes B and C | An effect only if the groups differ | `log_or` about 0 and `p = 1`, by construction: each pair's value enters both groups |
| The same genes in two comparison files, default `min_pairs = 2` | Each gene combines two values | No results and no files: `gene_key` includes the comparison, so values never combine |
| Sequence names without a subgenome letter | Usable with `group_mode = "custom"` | Stops with a clear message by default, which is reasonable |

`gene_out_dir` also appears not to control where output is written: files landed in the working directory.

## From code review, not yet run

- One regions file is applied to both sides, so a subject gene can match a query region on a same-named sequence.
- BED starts are not converted, unlike `dnds_summary`.
- Zero-cell contrasts give no row, and mixed-model warnings are suppressed.

## Why the design matters

The question the mode seems meant to answer is whether genes in one group (subgenome, region) are more often in a selection state than genes in another. A pairwise dN/dS has one value shared by two genes, so it cannot be independent evidence for both, and comparing two ends of the same table measures nothing. Valid group comparisons need each gene's value to come from comparisons that distinguish the groups.

## Defensible options

1. **Paired, within-triad.** For homeolog triads, compare each gene's dN/dS against the same outgroup, then test within-triad differences. Each triad is the unit, so shared history is controlled. This needs a triad definition. It is the recommended starting point.
2. **Mixed model with a pair or gene random effect,** on rows that each use one side of a pair, chosen by an explicit rule. It is flexible, but the rule decides the answer.
3. **Remove gene level,** and document how to do this with pairwise mode plus a user-defined `group` column. This is simplest, and keeps the package honest until a design is chosen.

Whichever is chosen needs simulation under a known null to confirm false-positive rates before use.

For v1.0.0: choose option 1 or 3.
