# Tutorial: testing for episodic selection across species

`phylo_prepare` and `phylo_busted` go from genomes and orthogroups to HyPhy BUSTED tests of whether some sites evolved with dN/dS above 1 on branches you choose. You do not build alignments or HyPhy inputs by hand. This is separate from dndsR's pairwise dN/dS. Read [phylo-methods.md](../phylo-methods.md) before interpreting results.

**Status: experimental, new in v0.5.0.** dndsR prepares inputs as designed and reproduces HyPhy run directly ([evidence](../validation/phylo-prepare.md)). A simulation study on one five-taxon tree recommends `--multiple-hits Double` and three starts, and treats masking as a sensitivity analysis; the defaults are unchanged. See [calibration](../phylo-methods.md#calibration-in-simulation) before choosing settings.

## 1. Get dndsR and HyPhy

You need a dndsR clone that has the commands, and HyPhy 2.5.103. Check the clone:

```
grep -q phylo_prepare NAMESPACE && echo "this clone has phylo_prepare"
```

Published dndsR images up to v0.4.0 include MAFFT but not HyPhy, so install HyPhy separately. With micromamba, which this was tested with (conda or mamba take the same arguments, untested):

```
micromamba create -y -p "$HOME/.dndsr/hyphy" -c conda-forge -c bioconda hyphy=2.5.103
"$HOME/.dndsr/hyphy/bin/hyphy" --version
```

The version line starts `HYPHY 2.5.103`. Then install dndsR from the clone, as in [getting started](getting-started.md), making HyPhy's folder visible in the container:

```
tools/dndsr-install install --engine singularity --sif dndsr.sif --bind "$HOME/.dndsr/hyphy"
export PATH="$HOME/.local/bin:$PATH"
dndsr --help
```

The command list includes `phylo_prepare` and `phylo_busted`. In an image or environment with `hyphy` on its `PATH`, leave out `--hyphy` below.

## 2. Run the example: genomes to results

The example has six small simulated genomes (FASTA and GFF3), 14 orthogroups, a species tree, and a foreground: the clade of spA and spB, with its stem. Eight orthogroups can be tested; six are designed to be excluded.

```
dndsr phylo_example
cd dndsr_phylo_example
dndsr phylo_prepare --genome-file genomes.tsv --orthogroup-file orthogroup_table.tsv \
  --species-tree-file species_tree.nwk --foreground-file foreground.tsv -t 4
dndsr phylo_busted --prepared-dir prepared --hyphy "$HOME/.dndsr/hyphy/bin/hyphy" -t 4
dndsr verify_phylo_example
```

Preparation takes about 20 seconds and reports `prepared 8, excluded 6, failed 0`; BUSTED takes about 40 seconds with 4 threads. The check ends with `All 11 phylogenetic example checks passed.` That confirms dndsR built the alignments as designed and reproduces HyPhy, which the check runs directly on this machine on inputs prepared without dndsR. It says nothing about your data.

`phylo_prepare` extracted each genome's CDS (in `prepared/cds/`), checked them, aligned the proteins with MAFFT, placed the codons under them, pruned the species tree for each orthogroup, and labelled its foreground. Its outputs are in `prepared/`: `alignments/`, `trees/` (labelled, as HyPhy will see them), `prepared_summary.tsv`, `prepared_members.tsv` (every gene of every orthogroup, used or not), `prepared_sequences.tsv`, and `prepared_provenance.json`.

## 3. Read the results

`busted/busted_results.tsv` has one row per orthogroup:

| orthogroup | status | stage | reason | exclusion_basis | p_value | p_adj |
|---|---|---|---|---|---:|---:|
| OG01 | tested | analysis | | | 0.0200 | 0.0799 |
| OG02 | tested | analysis | | | 0.105 | 0.281 |
| OG03 | tested | analysis | | | 0.00033 | 0.0027 |
| OG04, OG05, OG06, OG09, OG10 | tested | analysis | | | 0.309 or 0.5 | 0.5 |
| OG07 | excluded | prepare | multi_copy | policy | | |
| OG08 | excluded | prepare | cds_invalid | data_quality | | |
| OG11 | excluded | prepare | foreground_absent | no_data | | |
| OG12 | excluded | prepare | foreground_partly_unidentifiable | policy | | |
| OG13 | excluded | prepare | cds_invalid | data_quality | | |
| OG14 | excluded | prepare | foreground_unidentifiable | mathematical | | |

- **Tested.** `p_adj` is corrected across the 8 tested orthogroups (`family_size`). OG01 to OG03 were simulated with selection on the foreground; only OG03 is significant after correction. The others were simulated without; a p-value of exactly 0.5 means no evidence at all. One small example says nothing about power or false positives.
- **Excluded, and why.** OG07 has two spD genes (a paralog). OG08 and OG13 have a gene model that `extract_cds` rejected (an in-frame stop; a frameshift). OG11 has neither foreground taxon. OG12 lacks spC, so the clade's stem joins the branch above it: spA and spB remain separable, but testing only them would narrow the hypothesis. OG14 lacks spB and spC, so nothing in the foreground can be separated from background history. `detail` explains each.
- OG09 lacks spB: under a `clade` foreground, spA's branch, joined with the stem, still carries the hypothesis, and the flag `foreground_branches_absent:1` records it.
- `busted/busted_branches.tsv` lists the branches HyPhy tested. `prepared/trees/` has the labelled trees.

Results are for HyPhy 2.5.103, seed 1, and 1 HyPhy thread, on a processor with AVX-512; without it, OG01's p-value is 0.0199 rather than 0.0200. Another seed, thread count, or processor can change the values slightly; `--extra-seeds 2,3` reports how much, without changing `p_value`.

## 4. Gene trees

The example has gene trees for OG04 to OG06, with tips named as OrthoFinder names them. With `--gene-tree-dir`, each orthogroup needs a gene tree; the others are excluded as `gene_tree_missing`.

```
dndsr phylo_prepare --genome-file genomes.tsv --orthogroup-file orthogroup_table.tsv \
  --species-tree-file species_tree.nwk --foreground-file foreground.tsv \
  --gene-tree-dir gene_trees --run-name prepared_gt -t 4
dndsr phylo_busted --prepared-dir prepared_gt --run-name busted_gt --hyphy "$HOME/.dndsr/hyphy/bin/hyphy" -t 4
```

OG04's gene tree differs from the species tree away from the foreground, and OG05's beside it; both are tested, with `gt_*` columns and flags saying so. OG06's gene tree has no branch for spA and spB together, so the hypothesis cannot be placed on it: `foreground_not_in_gene_tree`. Add `--gene-tree-mode report` to test on the species tree and only report discordance.

## 5. Your own data

1. **Genomes** in a table with columns `taxon`, `fasta`, `gff`, or **CDS** in a table with `taxon`, `cds` (use `--cds-file`). Taxon names: letters, digits, and underscores, starting with a letter.
2. **Orthogroups**: a table with `orthogroup`, `taxon`, `gene_id`, or OrthoFinder's `Orthogroups.tsv` or `N0.tsv`. Gene IDs must match the CDS names: run OrthoFinder on the proteins from `dndsr extract_cds --export-proteins`, and add an `orthofinder_name` column if its species columns are not your taxon names.
3. **A rooted species tree** and **a foreground file**, tab-separated, decided before you look at results:

   ```
   scope	taxa
   stem	spA,spB
   ```

   Use `terminal` for a taxon's own branch, `stem` for the branch leading to a clade, or `clade` for the stem and everything in it.

```
dndsr phylo_prepare --genome-file genomes.tsv --orthogroup-file Orthogroups.tsv \
  --species-tree-file species.nwk --foreground-file foreground.tsv -t 16 -o results
dndsr phylo_busted --prepared-dir results/prepared -o results -t 16 \
  --multiple-hits Double --extra-seeds 2,3
```

These are the settings the [calibration](../phylo-methods.md#calibration-in-simulation) recommends, not HyPhy's defaults, which the example above used: `--multiple-hits Double` keeps multinucleotide changes from posing as selection (97% false positives without it in simulation, 0-2% with it), and `--extra-seeds 2,3` runs three starts. Report `p_best_fits_adj`, the Benjamini-Hochberg correction of `p_best_fits` (best alternative against best null fit over the three starts), not `p_adj`, which uses the first start only. A gene near your threshold remains uncertain: about 6% of calls near it differed between independent runs. None of this protects against alignment errors, the main source of false positives in the calibration; look at alignments flagged `dissimilar_window` before believing a result.

`-t` sets how many orthogroups run at once; each HyPhy run uses `--hyphy-cpu` threads (default 1). Results do not depend on `-t`, but do depend on `--seed` and `--hyphy-cpu`, which are recorded. A finished orthogroup is not rerun unless its inputs, settings, or HyPhy change; `--overwrite` forces it.

Before interpreting, count orthogroups by status and reason (both commands print them), and read the [limits](../phylo-methods.md#limits-and-interpretation). A significant result is evidence under the model for episodic positive selection somewhere on the foreground, not proof, and not a statement about particular sites.

Codon alignments made elsewhere can be tested directly: `dndsr phylo_busted --alignment-dir <folder> --species-tree-file <tree> --foreground-file <file>`. `dndsr phylo_example --type alignments` has an example.

## From R

The same steps as functions. MAFFT and HyPhy must be on the `PATH`, or named with `mafft =` and `hyphy =`.

```r
library(dndsR)
folder <- phylo_example(file.path(tempdir(), "dndsr_phylo_example"))
setwd(folder)
```

```r
prepared <- phylo_prepare("orthogroup_table.tsv", "species_tree.nwk", "foreground.tsv",
                          genome_file = "genomes.tsv", threads = 2)
res <- phylo_busted(prepared_dir = "prepared", threads = 2)
verify_phylo_example()
```

```r
table(res$status, res$exclusion_basis, useNA = "ifany")
res[res$status == "tested", c("orthogroup", "p_value", "p_adj", "qc_flags")]
```
