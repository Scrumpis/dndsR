# Tutorial: testing for episodic selection across species

`phylo_busted` runs HyPhy's BUSTED on orthogroups of several species, to test whether some sites evolved with dN/dS above 1 on branches you choose. It is separate from dndsR's pairwise dN/dS. Read [phylo-methods.md](../phylo-methods.md) before interpreting results.

**Status.** This is development code, not in any release: v0.4.0 does not have it. It reproduces HyPhy run directly, but has not been calibrated as a selection pipeline ([evidence](../validation/phylo-busted-wrapper.md)). In this phase you supply the codon alignments; dndsR does not yet build them.

## 1. Get dndsR and HyPhy

You need a dndsR clone that has the command, and HyPhy 2.5.103. Check the clone:

```
grep -q phylo_busted NAMESPACE && echo "this clone has phylo_busted"
```

Published dndsR images up to v0.4.0 do not include HyPhy, so install it separately. With micromamba, which this was tested with (conda or mamba take the same arguments, untested):

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

The command list includes `phylo_busted`. In an image or environment that already has `hyphy` on its `PATH`, leave out `--hyphy` below.

## 2. Run the example

The example has seven simulated orthogroups of up to eight taxa. Four can be tested; three are designed to be skipped.

```
dndsr phylo_example
cd dndsr_phylo_example
dndsr phylo_busted --alignment-dir alignments --species-tree-file species_tree.nwk \
  --foreground-file foreground.tsv --hyphy "$HOME/.dndsr/hyphy/bin/hyphy" -t 4
dndsr verify_phylo_example
```

This takes about a minute with 4 threads. It reports `tested 4, skipped 3, failed 0`, and the check ends with `All 9 phylogenetic example checks passed.` That confirms dndsR reproduces HyPhy run directly on this example. It says nothing about your data.

## 3. Read the results

`busted/busted_results.tsv` has one row per orthogroup:

| orthogroup | status | reason | p_value | p_adj |
|---|---|---|---:|---:|
| OG0001 | tested | | 0.0129 | 0.0518 |
| OG0002 | tested | | 0.5 | 0.5 |
| OG0003 | tested | | 0.197 | 0.393 |
| OG0004 | tested | | 0.411 | 0.5 |
| OG0005 | skipped | foreground_absent | | |
| OG0006 | skipped | internal_stop | | |
| OG0007 | skipped | multi_copy | | |

- `p_adj` is corrected across the 4 tested orthogroups (`family_size`). OG0001, simulated with selection on the foreground, has p = 0.013, and p_adj = 0.052 after correction.
- OG0002 has p = 0.5 exactly: its likelihood ratio is 0, so there is no evidence at all.
- OG0003 was also simulated with selection, and is not significant. One gene says nothing about power.
- OG0004 lacks spB, so spA's branch includes the clade's stem. Under the `clade` foreground that is still foreground. Under `stem`, it would be skipped as `foreground_not_identifiable`.
- Check `qc_flags`, and `busted/busted_branches.tsv` to see which branches HyPhy tested.

Results with HyPhy 2.5.103, seed 1, and 1 HyPhy thread. Another seed or thread count can change the values slightly.

## 4. Your own data

You need:

1. **Codon alignments**, one in-frame FASTA per single-copy orthogroup, named `<orthogroup>.fasta`, with sequences named `taxon` or `taxon|gene_id`.
2. **A rooted species tree** of the taxa, in Newick. Names: letters, digits, and underscores, starting with a letter.
3. **A foreground file**, tab-separated, decided before you look at results:

   ```
   scope	taxa
   stem	spA,spB
   ```

   Use `terminal` for a taxon's own branch, `stem` for the branch leading to a clade, or `clade` for the stem and everything in it.

```
dndsr phylo_busted --alignment-dir my_alignments --species-tree-file species.nwk \
  --foreground-file foreground.tsv -t 16 -o results --run-name busted_stemAB
```

`-t` sets how many orthogroups run at once; each HyPhy run uses `--hyphy-cpu` threads (default 1). Results do not depend on `-t`. They do depend on `--seed` and `--hyphy-cpu`, which are recorded. A finished orthogroup is not rerun unless its inputs, settings, or HyPhy change; add `--overwrite` to force it.

Before interpreting:

- Count skipped orthogroups by reason (the run prints them). Many `foreground_not_identifiable` means the foreground taxa or their sister lineages are often missing.
- Read the [limits](../phylo-methods.md#limits-and-interpretation). A significant result is evidence under the model for episodic positive selection somewhere on the foreground, not proof, and not a statement about particular sites.

## From R

The same steps as functions. HyPhy must be on the `PATH`, or named with `hyphy =` or `options(dndsR.hyphy = ...)`.

```r
library(dndsR)
folder <- phylo_example(file.path(tempdir(), "dndsr_phylo_example"))
setwd(folder)
```

```r
res <- phylo_busted("alignments", "species_tree.nwk", "foreground.tsv", threads = 2)
verify_phylo_example()
```

```r
res[, c("orthogroup", "status", "reason", "p_value", "p_adj", "qc_flags")]
branches <- read.delim("busted/busted_branches.tsv")
subset(branches, foreground)
```
