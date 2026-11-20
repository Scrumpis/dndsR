#' Relaxation or intensification of selection with HyPhy RELAX (experimental)
#'
#' Compares the strength of selection on test branches with that on reference branches,
#' with HyPhy's RELAX (Wertheim et al. 2015), in every orthogroup of a [phylo_prepare()]
#' folder or a folder of codon alignments. Experimental: dndsR reproduces RELAX run
#' directly, but its error rates with dndsR's inputs have been measured only in simulation
#' (`docs/phylo-methods.md`).
#'
#' @section What it tests:
#' RELAX fits omega classes on reference branches and the same classes raised to a power
#' `k` on test branches. `k < 1` moves every class towards 1: **relaxation** (weaker
#' purifying and weaker positive selection). `k > 1` moves them away from 1:
#' **intensification**. The test compares `k = 1` with a free `k`. It is not a test of
#' positive selection: intensification can mean stronger purifying selection, and a
#' gene under positive selection can show either. Read `k` with the omega classes.
#'
#' @section Branch sets:
#' The foreground file's rows are the test set. Rows with `set` = `reference` give the
#' reference set; without them, every branch outside the test set is the reference. Both
#' sets must remain identifiable after pruning, and may not share a branch; otherwise the
#' orthogroup is excluded with a reason (`reference_*`, `test_reference_overlap`). Branches
#' in neither set are unclassified, with their own parameters.
#'
#' @section Outputs:
#' In `file.path(output_dir, run_name)`: `_results.tsv` (one row per orthogroup: `k`,
#' `direction`, `lrt`, `p_value`, `p_adj`, likelihoods, and the omega classes of each set),
#' `_branches.tsv` (with each branch's `set`), `_sequences.tsv`, `_provenance.json`, and
#' `hyphy/<orthogroup>/`. `p_adj` is Benjamini-Hochberg across tested orthogroups.
#'
#' @inheritParams phylo_busted
#' @param models "All" (HyPhy's default, also fits descriptive models) or "Minimal" (only
#'   the null and alternative, faster).
#' @param srv Synonymous rate variation: "No" (HyPhy's default for RELAX) or "Yes".
#' @return Invisibly, the results table.
#' @references Wertheim JO, et al. (2015) RELAX: detecting relaxed selection in a
#'   phylogenetic framework. Mol Biol Evol 32:820-832.
#' @seealso [phylo_prepare()], [phylo_busted()], [phylo_absrel()], [phylo_summarize()]
#' @export
phylo_relax <- function(prepared_dir = NULL,
                        alignment_dir = NULL,
                        species_tree_file = NULL,
                        foreground_file = NULL,
                        gene_tree_dir = NULL,
                        gene_tree_mode = NULL,
                        output_dir = ".",
                        run_name = "relax",
                        models = "All",
                        rates = 3L,
                        srv = "No",
                        multiple_hits = "None",
                        grid_size = 250L,
                        starting_points = 1L,
                        genetic_code = "Universal",
                        seed = 1L,
                        hyphy_cpu = 1L,
                        threads = 1L,
                        min_taxa = 4L,
                        terminal_stops = c("mask", "exclude"),
                        overwrite = FALSE,
                        hyphy = NULL) {
  terminal_stops <- match.arg(terminal_stops)
  settings <- .hyphy_settings("relax", list(models = models, rates = rates, srv = srv, multiple_hits = multiple_hits,
                                            grid_size = grid_size, starting_points = starting_points, genetic_code = genetic_code))
  .phylo_run("relax", alignment_dir, species_tree_file, foreground_file, gene_tree_dir, output_dir, run_name,
             settings, seed, hyphy_cpu, threads, min_taxa, terminal_stops, overwrite, hyphy,
             gene_tree_mode = gene_tree_mode, prepared_dir = prepared_dir)
}
