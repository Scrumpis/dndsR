# Phylogenetic selection module: design

Approved by the maintainer on 2026-11-13: HyPhy as an external engine, BUSTED first, phases P0 and P1 now, stop before P2. This records the design and the P0 decisions. User documentation: [phylo-methods.md](../phylo-methods.md). Evidence: [phylo-busted-wrapper.md](../validation/phylo-busted-wrapper.md).

## Scope and separation

- dndsR wraps HyPhy and does not reimplement its models. The module is optional: `ape` is in `Suggests`, and HyPhy is a system requirement checked at run time.
- Pairwise code is untouched. No pairwise function, default, or output changed. Phylogenetic results are written to their own folder, `<output_dir>/<run_name>/`, have no `dN`, `dS`, or `dNdS` columns, and say which method made them.
- A significant BUSTED test and a pairwise dN/dS above 1 answer different questions and are never combined or compared in one table.
- The v1.0.0 roadmap is unchanged. This module is a separate line of work, and is not in its priorities.

## Layers

Everything below the method registry is shared, so aBSREL and RELAX add a registry entry and an exported function, not new infrastructure.

| Layer | File | Method-specific? |
|---|---|---|
| Alignment input and QC | `R/phylo_alignment.R` | No |
| Species tree, foreground hypotheses, pruning, identifiability, gene trees, Newick | `R/phylo_tree.R` | No |
| HyPhy discovery, version, isolated runs with seed and threads, JSON reading | `R/phylo_hyphy.R` | No |
| Method registry: settings and defaults, arguments, required JSON keys, parser, unit of test | `.hyphy_methods` in `R/phylo_hyphy.R` | Yes, one entry per method |
| Runner: per-orthogroup preparation, reuse, execution, label check, correction, provenance | `R/phylo_run.R` | No |
| Exported analysis | `R/phylo_busted.R` | Yes |
| Example and check | `R/phylo_example.R` | Example is BUSTED-specific |

## Data structures

- **Sequences.** One row per sequence: `seq_name`, `taxon`, `gene_id`, `sequence`. Names are `taxon` or `taxon|gene_id`. Several sequences per taxon are representable, so duplication-aware analyses can use the same structure later; P1 skips them as `multi_copy`.
- **Hypothesis.** The species tree, one row per reference branch with its descendant taxa, whether it is foreground, and whether it must exist in a gene tree (terminal and stem branches); clade sets for clade scope. Resolved once per run.
- **Branch map.** Per orthogroup, one row per branch given to HyPhy: name, terminal or not, the smaller side of its split, foreground or not, the species-tree branches it combines, and whether the topology came from the species tree or a gene tree.
- **Results.** One row per orthogroup considered, with `status` (`tested`, `skipped`, `failed`) and `reason`. Per-sequence QC is a separate table. A future branch-level method (aBSREL) adds a branch-level results table keyed by orthogroup and branch, with its own correction.

## Branch identity and identifiability

A branch is identified by the bipartition of the taxa present, never by a node name. This is what keeps a hypothesis's meaning after pruning.

1. The foreground file is resolved on the full, rooted species tree. `terminal`: each listed taxon's branch. `stem`: the branch above the most recent common ancestor of the listed taxa, which must be exactly that clade. `clade`: that stem and every branch below it. Stem and clade need a rooted tree; the root has no stem.
2. For an orthogroup, each species-tree branch maps to the split it induces on the taxa present. Branches with the same split become one branch once taxa are pruned. This includes the two branches at the root, which HyPhy joins, as it treats trees as unrooted (confirmed by running).
3. If a joined branch combines foreground and background branches, the hypothesis is **not identifiable** for that orthogroup, and it is skipped (`foreground_not_identifiable`), with the branches named. It is never relabelled. Consequences, all tested:
   - A stem hypothesis needs every clade member present, and the clade's sister lineage.
   - A clade hypothesis tolerates missing members, but not a missing sister lineage, because the stem would then include time before the clade.
   - A terminal hypothesis needs the taxon's sister lineage.
   - The stem of a clade beside the root is never identifiable, because it joins the other root branch.
4. No foreground branch left: `foreground_absent`.
5. With a gene tree, the foreground is still resolved on the species tree, so missing taxa are judged the same way. Every terminal and stem foreground branch, including each clade's stem, must exist in the gene tree, or the orthogroup is `foreground_discordant`. Branches inside a clade-scope foreground are foreground, even where the gene tree disagrees with the species tree inside the clade.

## Alignment policy

| Kind | Cases | Action |
|---|---|---|
| Exclusion | empty, unknown taxon, more than one sequence per taxon, unequal lengths, length not a multiple of 3, characters other than nucleotides, IUPAC codes, `-`, or `?`, a codon with one or two gap characters, a stop codon before the end of a sequence, a sequence with no determined codon, fewer than `min_taxa` sequences | Orthogroup skipped with the reason |
| Documented change | A stop codon that is the last non-gap codon of a sequence | Replaced by `---` and counted (`terminal_stops = "mask"`, default), or the orthogroup excluded (`"exclude"`). HyPhy refuses stop codons, so most CDS-derived alignments would otherwise fail |
| Reported only | Ambiguous codons, gappy sequences and columns, a sequence far more divergent than the rest, a high pairwise p-distance, a window of 10 codons mostly unlike every other sequence | Metrics and flags. Thresholds are heuristics, not validated, and never exclude or change data |

Nothing else is repaired, masked, trimmed, or realigned.

## HyPhy

- **Version.** 2.5.103: the newest GitHub release (2026-10-06), on Bioconda since 2026-10-09. 2.5.104 was tagged 2026-11-06 with no release notes; it gave bit-identical results on the example. HyPhy 2.5 is now in maintenance mode, with development moved to HyPhy 3 (first release v0.0.1, 2026-11-09, not on Bioconda). HyPhy 3 is not suitable yet; revisit when it has a stable, published release.
- **Defaults.** Every BUSTED setting defaults to HyPhy's and is passed explicitly: `--srv Yes`, `--rates 3`, `--syn-rates 3`, `--multiple-hits None`, `--error-sink No`, `--grid-size 250`, `--starting-points 1`, read from `BUSTED.bf` in 2.5.103. Alternatives are exposed. No dndsR default differs from HyPhy's.
- **Labels.** `{FG}` on branches and `--branches FG`. HyPhy's `tested` field is compared with dndsR's labels on every run; a difference fails the orthogroup.
- **Topology only.** Branch lengths are removed: HyPhy uses given lengths as starting values (results changed with them at a fixed seed), and species-tree lengths are not in its units.
- **Reproducibility.** HyPhy seeds itself randomly unless `RANDOM_SEED` is set, and results depend on the seed and on the thread count. dndsR always passes `ENV=RANDOM_SEED=<seed>` (default 1) and `CPU=<hyphy_cpu>` (default 1), and records both. Same build, seed, and threads gave bit-identical results.
- **Isolation.** Each orthogroup runs in its own folder, because HyPhy writes `errors.log` to its working folder.
- **Reuse.** A finished orthogroup is reused only if its prepared inputs, settings, seed, threads, and HyPhy version match its `record.json`.

## Statistics

- p-values are HyPhy's: `0.5 * (1 - CChi2(LRT, 2))` in `BUSTED.bf`, a 50:50 mixture of chi-square with 0 and 2 degrees of freedom. A likelihood ratio of 0 gives exactly 0.5, so p-values have a point mass at 0.5; QQ plots must allow for it.
- Benjamini-Hochberg across one family: every orthogroup with status `tested` in one run, under one foreground hypothesis and one set of settings. Raw p-values are kept. Skipped and failed orthogroups are outside the family, counted, and listed.
- Ascertainment: only single-copy orthogroups that pass QC and keep the foreground identifiable are tested. The provenance records counts by reason and, for each taxon, the number of tested orthogroups containing it.

## Phases

| Phase | Content | Status |
|---|---|---|
| P0 | Decisions, design | Done, this file and `DECISIONS.md` |
| P1 | `phylo_busted` on user-supplied codon alignments, with optional gene trees; QC; provenance; example; tests; container and conda changes, not published | Done; see the validation page |
| P2 | Orthogroup preparation from `extract_cds` output: orthogroups, MAFFT protein alignment, back-translation, QC | Not started |
| P3 | Scientific calibration: simulation, published controls, optimizer variability, multiple-hit settings | Not started |
| P4 | aBSREL: branch-level results and their correction | Not started |
| P5 | RELAX | Not started |
| P6 | Enrichment over the tested set, and plots | Not started |

## Open decisions for the maintainer

1. **Optimizer variability.** Over 20 seeds, one gene's LRT ranged 7.0 to 8.8, and its best unconstrained log-likelihood by 0.59. A gene near the threshold can change call with the seed. Options for P3: accept and report; several seeds per gene with the best likelihood per model; or more `--starting-points`. Each changes results, so it needs evidence and approval.
2. **Multiple hits.** HyPhy's default is `None`. Whether dndsR should recommend `Double` is for P3 simulation, not now.
3. **BUSTED-E (`--error-sink`)** as a standard sensitivity analysis, given HyPhy calls it experimental.
4. **Strict identifiability.** A missing sister lineage makes stem and clade hypotheses untestable. This is correct, but may exclude many orthogroups in patchy data. The alternative, accepting stems that include earlier time, changes the hypothesis.
5. **Gene trees.** Whether discordance near the foreground should exclude (as now), be flagged, or be handled otherwise.
6. **Development version.** `DESCRIPTION` stays 0.4.0, so provenance from this code says 0.4.0. A development version such as 0.4.0.9000 would distinguish it, but also changes the skill's `software-version`.
7. **Release line.** Whether this ships as 0.5.0 before 1.0.0, and the image that carries HyPhy.
