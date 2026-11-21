# Tutorial: comparing selection among lineages

From species genomes, orthogroups, and a species tree to branch-level dN/dS, branch tests of episodic selection, tests of relaxed or intensified selection, and one table per lineage that you can join to lineage traits such as LTR activity or genome-size change. Read [lineage-selection.md](../lineage-selection.md) for what each method tests before interpreting anything.

**Status: experimental, development code.** dndsR reproduces codeml and HyPhy run directly; the methods' error rates were measured only in simulation ([results](../validation/phylo-lineage-calibration.md)), where RELAX on 16 taxa and the free-ratio test with missing taxa gave more false positives than their level. The example below is synthetic: its numbers show how to read the outputs, not anything about real genomes.

## 1. Get dndsR, HyPhy, and PAML

You need a dndsR clone with these commands (`dndsr --help` lists `phylo_branch_omega`), HyPhy 2.5.103, and PAML 4.10.10. Images built from this version's `tools/Dockerfile` include both. Otherwise install them with micromamba (conda or mamba take the same arguments, untested):

```
micromamba create -y -p "$HOME/.dndsr/phylo" -c conda-forge -c bioconda hyphy=2.5.103 paml=4.10.10
```

and give `--hyphy "$HOME/.dndsr/phylo/bin/hyphy"` and `--codeml "$HOME/.dndsr/phylo/bin/codeml"` to the commands below, or put that `bin` folder on `PATH`. With the container route, bind it at install time (`tools/dndsr-install install ... --bind "$HOME/.dndsr"`), as in the [phylogenetic tutorial](phylogenetic-selection.md).

## 2. Run the example

Six small simulated genomes, 14 orthogroups, a species tree with `spF` as the outgroup, and a foreground: the clade of `spA` and `spB`.

```
dndsr phylo_example -o lineage_demo
cd lineage_demo
dndsr phylo_prepare --genome-file genomes.tsv --orthogroup-file orthogroup_table.tsv \
  --species-tree-file species_tree.nwk --foreground-file foreground.tsv -t 4
dndsr phylo_branch_omega --prepared-dir prepared -t 4
dndsr phylo_absrel --prepared-dir prepared -t 4
dndsr phylo_relax --prepared-dir prepared -t 4 --models Minimal
dndsr phylo_summarize --branch-omega-dir branch_omega --absrel-dir absrel --relax-dir relax
```

`phylo_prepare` builds codon alignments once; the three analyses reuse them. Eight orthogroups are analysed and six are excluded by design, with reasons, as in the phylogenetic tutorial.

## 3. Read the results

Values are for this example with PAML 4.10.10 and HyPhy 2.5.103; another machine can differ slightly in HyPhy's results (see the [phylogenetic tutorial](phylogenetic-selection.md)).

**Branch models**, `branch_omega/branch_omega_results.tsv`. OG01 to OG03 were simulated with episodic positive selection on the foreground clade:

| orthogroup | omega_bg | omega_fg | p_adj_two_vs_m0 | p_adj_fg_omega_gt1 |
|---|---:|---:|---:|---:|
| OG01 | 0.335 | 0.812 | 0.025 | 1.000 |
| OG02 | 0.369 | 1.666 | 0.002 | 0.409 |
| OG03 | 0.285 | 0.893 | 0.007 | 1.000 |
| OG04 | 0.508 | 0.243 | 0.207 | 1.000 |

The foreground's gene-wide omega is higher in the three selected genes, significantly after correction, yet only in OG02 does it exceed 1, and not significantly. That is typical: selection at a few sites raises the gene-wide average without lifting it above 1. Read branch-model omega as an average, never as a test of positive selection.

**Lineages**, `summary/summary_lineages.tsv` (selected columns):

| branch_id | n_estimated | pooled_omega (interval) | median_omega_reliable | absrel_n_sig_lineage |
|---|---:|---|---:|---:|
| spA | 7 | 0.90 (0.42-1.53) | 1.16 | 3 of 7 |
| spB | 7 | 0.80 (0.52-1.32) | 0.69 | 0 of 7 |
| mrca(spA,spB) | 7 | 0.48 (0.26-0.69) | 0.58 | 0 of 7 |
| spC | 8 | 0.30 (0.17-0.50) | 0.25 | 0 of 8 |
| spD | 8 | 0.44 (0.32-0.64) | 0.37 | 0 of 8 |
| spF | 0 | NA | NA | 0 |

- `pooled_omega` sums substitutions over genes; its interval resamples genes. The median of per-gene omegas answers a different question (the typical gene) and can differ, as for `spA`.
- `absrel_n_sig_lineage` counts genes with aBSREL evidence of episodic selection on that branch after the per-lineage correction (`p_bh_lineage`).
- `spF` and `mrca(spA,spD)` are the two branches at the root. Every analysis uses an unrooted tree, so they always join, and no estimate can be attributed to either; both are left out (`n_estimated` 0). This is why your species tree needs an outgroup.
- `summary_lineage_pairs.tsv` compares lineages on the genes both share; `summary_foreground.tsv` counts 3 of 8 genes with a significantly higher foreground omega, none significantly lower.

**RELAX**, `relax/relax_results.tsv`: 3 genes with significant intensification and 1 with relaxation after correction. The relaxed one, OG02, has `k = 0` and `k_at_bound` TRUE: an extreme fit at HyPhy's bound, not a measured strength. Intensification here fits the simulation, since selection pushed some omega classes far above 1; RELAX is not a test of positive selection, and intensified purifying selection would look the same.

## 4. Join lineage traits and fit a comparative model

`summary/summary_branch_ids.tsv` lists every species-tree branch with `ape_node` (for ancestral states from `ape::ace` or `phytools`) and `tree_label`; `summary/summary_species_tree_ids.nwk` is the species tree labelled with it. The example includes a synthetic trait table, `lineage_traits.tsv` (made-up LTR fractions and genome sizes). In R:

```r
lin <- read.delim("summary/summary_lineages.tsv")
traits <- read.delim("lineage_traits.tsv")
species <- merge(lin[lin$terminal, ], traits, by.x = "branch_id", by.y = "taxon")
species <- species[is.finite(species$pooled_omega), ]
species[, c("branch_id", "pooled_omega", "n_pooled", "ltr_fraction")]
tree <- ape::read.tree("summary/summary_species_tree_ids.nwk")
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, species$branch_id))
tree$edge.length <- rep(1, nrow(tree$edge)) # the example's species tree has no branch lengths
fit <- nlme::gls(log(pooled_omega) ~ ltr_fraction, data = species,
                 correlation = ape::corBrownian(phy = tree, form = ~branch_id))
round(summary(fit)$tTable, 3)
```

This joins five species (spF has no estimate) and fits a phylogenetic generalised least squares model: `ltr_fraction` 4.09, standard error 1.53, p = 0.075. **These numbers mean nothing**: the traits are invented, and five species cannot support such a model. The block shows the mechanics. In a real analysis use dated or substitution-scaled branch lengths, many more species, weights or sensitivity checks for lineages with few genes (`n_pooled`) or saturated dS (`frac_high_ds`), and internal branches (joined by `branch_id`) when the trait is a change along a branch, such as genome-size change.

## 5. Your own data

1. **Genomes or CDS, and orthogroups**, as in the [phylogenetic tutorial](phylogenetic-selection.md#5-your-own-data).
2. **A rooted species tree with an outgroup and labelled internal nodes**, for example `((((spA,spB)ab,spC)abc,(spD,spE)de)ingroup,spF);`. The labels become branch identifiers everywhere. The two branches at the root join in every analysis and cannot be attributed to one lineage, so keep lineages of interest away from the root.
3. **A foreground file** for the two-ratio model and RELAX, decided before you look at results. For RELAX, add a column `set`; rows with `reference` give the reference set, otherwise every other branch is the reference:

   ```
   scope	taxa	set
   clade	spA,spB	test
   clade	spD,spE	reference
   ```

Then the same commands, on your prepared folder. `-t` sets how many orthogroups run at once. codeml fits each model from three starting values (`--omega-starts`), and aBSREL on every branch costs the most; on the public Brassicaceae run below, 300 orthogroups of five taxa took 9 minutes (branch models), 5 (aBSREL), and 5 (RELAX, `Minimal`) with `-t 16`, about 5 core-hours in all. Time grows quickly with the number of taxa: in simulation, free-ratio fits on 16 taxa took about 15 core-minutes per gene.

## What a real run looks like

To show the methods on real data, all four ran on 300 orthogroups sampled at random (seed 20261120) from single-copy orthologs of five public RefSeq genomes: *Arabidopsis thaliana* (GCF_000001735.4), *A. lyrata* (GCF_000004255.2), *Capsella rubella* (GCF_000375325.1), *Eutrema salsugineum* (GCF_000478725.1), and *Brassica rapa* (GCF_000309985.2), with the species tree `(((At,Al),Cr),(Es,Br))`, foreground the *Arabidopsis* clade, codon alignments from `phylo_prepare`. Results are outside the repository; this is a demonstration, not a study.

| Lineage | Genes estimated | pooled omega (interval) | median omega | aBSREL genes (per-lineage BH) |
|---|---:|---|---:|---:|
| At | 299 | 0.216 (0.194-0.240) | 0.186 | 20 of 298 |
| Al | 299 | 0.201 (0.175-0.232) | 0.165 | 13 of 298 |
| *Arabidopsis* stem | 299 | 0.249 (0.212-0.299) | 0.211 | 10 of 298 |
| Cr | 300 | 0.214 (0.183-0.258) | 0.177 | 18 of 299 |
| Es | 277 | 0.224 (0.199-0.257) | 0.199 | 12 of 276 |
| Br | 277 | 0.175 (0.156-0.194) | 0.169 | 22 of 276 |

- All 300 orthogroups ran under codeml; aBSREL and RELAX excluded one, where HyPhy deleted a foreground branch as zero-length.
- The root's two branches (the stem of *A. thaliana*, *A. lyrata*, and *C. rubella*, and that of *E. salsugineum* and *B. rapa*) join and have no estimate. Es and Br have fewer estimates because, where one of them is missing, the other's branch joins the root too.
- One gene in 300 had a significantly different foreground omega (two-ratio), none had a foreground omega above 1, and RELAX found 5 relaxed and 3 intensified after correction. Pooled omegas lie close together, with overlapping intervals except Br's, which is lower.
- codeml reported "check convergence" for at least one model in 56 genes (`codeml_warning:<model>` in `qc_flags`); the three starting values disagreed in 2. Look at these before relying on their estimates.

With five species, nothing here supports a lineage-level conclusion: the table shows what the outputs look like on real genes, and how many genes per lineage remain after pruning.
