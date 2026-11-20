# Lineage-level analyses: does dndsR run the engines correctly?

Engineering correctness of `phylo_branch_omega`, `phylo_absrel`, `phylo_relax`, and `phylo_summarize`: dndsR gives codeml and HyPhy the intended inputs, reads their outputs correctly, and maps branches correctly after pruning. Whether the methods' results mean what they say is a separate question: [phylo-lineage-calibration.md](phylo-lineage-calibration.md).

## Software

PAML 4.10.10 and HyPhy 2.5.103 from Bioconda, with the v0.4.0 image (R 4.6.2) and development code of the commit that adds this page. Runs on Slurm, HPC Linux nodes.

## Against the engines run directly

`tests/testthat/test-phylo-lineage-engines.R`, on the alignment example's OG0001 (all eight taxa) and OG0004 (no spB, so the foreground clade collapses onto spA's branch, merged with the clade's stem). Inputs for the direct runs were made without dndsR by `tools/dev/validation/phylo-direct-lineage.sh` (PHYLIP files by awk, trees written by hand) and are stored in `inst/extdata/phylo_example/truth/direct/`; the test writes codeml's control files and HyPhy's commands out in full.

| Engine | Compared | Result |
|---|---|---|
| codeml M0, free-ratio, two-ratio, best of three starting omegas | Log-likelihoods; M0 omega; background and foreground omega; every terminal branch's free-ratio omega, read from codeml's own labelled tree | Equal within 1e-3 log-likelihood units and 0.1% (the tolerance allows for the order of sequences and branches, which differs between the two routes) |
| aBSREL, every branch tested | Adaptive model's log-likelihood; branches tested; per terminal branch: uncorrected and Holm p-values, LRT, baseline omega | Identical (JSON against JSON); dndsR's table within 1e-12 (written to 15 digits) |
| RELAX, test = spA, spB and their stem, reference = the rest | `k`, LRT, p-value; direction | Identical |

For OG0004, the foreground of the two-ratio model is spA's branch, named `spA` with `species_tree_branches` joining spA and the clade's stem, as the hand-written direct tree has it.

Result at the commit that adds this page: 12 tests, 166 expectations, 0 failed, none skipped.

## Without the engines

`tests/testthat/test-phylo-lineage.R`: codeml's output is parsed from stored files of codeml 4.10.10 (values checked against the files by hand: likelihoods, parameters, omegas, and branches matched by the taxa below them); the PAML tree carries `#1` marks and no internal names; likelihood-ratio tests, negative ratios, and Benjamini-Hochberg corrections against hand-computed values; the pooled omega against a hand calculation; RELAX's test and reference sets after pruning (identifiable, absent, overlapping); aBSREL's all-branch labelling; per-lineage Benjamini-Hochberg families; branch identifiers and `ape` node numbers; matched lineage pairs on shared genes only.

## The bundled example, end to end

Every command of [the tutorial](../tutorials/lineage-selection.md) was run through the command line on the synthetic genome example (`phylo_example`), with HyPhy and PAML given by path: 8 orthogroups tested by each method, 6 excluded by design, none failed. The designed effects show: the three genes simulated with episodic selection on the foreground have a significantly higher foreground omega in the two-ratio model; aBSREL finds episodic selection on spA in 3 of 7 genes after the per-lineage correction; the root's two branches have no lineage estimate, as they must. The R block of section 4 runs.

The image checks (`tools/dev/container-checks.sh`) run the same four commands in the freshly built image, which has PAML and HyPhy inside, on every container build.

## Not covered here

Statistical calibration, real biological data beyond a demonstration, and other engine versions. A different PAML or HyPhy version is announced with a note when run.
