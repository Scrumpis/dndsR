#' Branch-specific episodic positive selection with HyPhy aBSREL (experimental)
#'
#' Tests each branch for episodic positive selection with HyPhy's aBSREL (Smith et al.
#' 2015), in every orthogroup of a [phylo_prepare()] folder or a folder of codon
#' alignments. Experimental: dndsR reproduces aBSREL run directly, but its error rates
#' with dndsR's inputs have been measured only in simulation (`docs/phylo-methods.md`).
#'
#' @section What it tests:
#' For each tested branch, whether some sites evolved with omega above 1 on that branch.
#' aBSREL fits a different number of omega classes per branch, chosen by AICc. A
#' significant branch is evidence under the model for episodic positive selection on some
#' sites of that branch, not across the gene, and not proof of adaptation.
#'
#' @section Branches tested:
#' `branches = "all"` (default) tests every branch: the exploratory analysis for comparing
#' lineages, which spends power on many tests. `"foreground"` tests only the foreground
#' hypothesis, as in [phylo_busted()], which is more powerful when the hypothesis was fixed
#' in advance. Branch identity follows the species tree: a branch of the pruned tree is
#' named by the species-tree branches it joins (`species_tree_branches`).
#'
#' @section Multiple testing:
#' Three levels, all reported. `p_raw` is HyPhy's uncorrected p-value for a branch.
#' `p_holm` is HyPhy's correction across the branches tested in that orthogroup: the k-th
#' smallest of m p-values times (m - k + 1), without Holm's running maximum, so it lies
#' between Holm's and Hochberg's procedures. Per orthogroup, `p_value` is the smallest
#' `p_holm`, in effect Hochberg's test that some tested branch was under selection (valid
#' when the branches' p-values are independent or positively dependent), and `p_adj`
#' applies Benjamini-Hochberg to it across
#' orthogroups. Per lineage, `p_bh_lineage` applies Benjamini-Hochberg to `p_raw` across the
#' orthogroups in which that species-tree branch was tested: use it to ask in how many genes
#' a given lineage shows episodic selection.
#'
#' @section Outputs:
#' In `file.path(output_dir, run_name)`: `_results.tsv` (one row per orthogroup),
#' `_branches.tsv` (one row per orthogroup and branch, with the baseline omega of every
#' branch and the test results of tested branches), `_sequences.tsv`, `_provenance.json`,
#' and `hyphy/<orthogroup>/` with HyPhy's inputs, log, and JSON.
#'
#' @inheritParams phylo_busted
#' @param branches "all" (default) or "foreground"; see Branches tested.
#' @param srv Synonymous rate variation: "No" (HyPhy's default for aBSREL) or "Yes".
#' @param multiple_hits "None" (HyPhy's default), "Double", or "Double+Triple".
#' @param foreground_file Needed only with `branches = "foreground"` and alignment input.
#' @return Invisibly, the results table.
#' @references Smith MD, et al. (2015) Less is more: an adaptive branch-site random effects
#'   model for efficient detection of episodic diversifying selection. Mol Biol Evol
#'   32:1342-1353.
#' @seealso [phylo_prepare()], [phylo_busted()], [phylo_relax()], [phylo_summarize()]
#' @export
phylo_absrel <- function(prepared_dir = NULL,
                         alignment_dir = NULL,
                         species_tree_file = NULL,
                         foreground_file = NULL,
                         gene_tree_dir = NULL,
                         gene_tree_mode = NULL,
                         output_dir = ".",
                         run_name = "absrel",
                         branches = c("all", "foreground"),
                         srv = "No",
                         multiple_hits = "None",
                         genetic_code = "Universal",
                         seed = 1L,
                         hyphy_cpu = 1L,
                         threads = 1L,
                         min_taxa = 4L,
                         terminal_stops = c("mask", "exclude"),
                         overwrite = FALSE,
                         hyphy = NULL) {
  branches <- match.arg(branches)
  terminal_stops <- match.arg(terminal_stops)
  settings <- .hyphy_settings("absrel", list(srv = srv, multiple_hits = multiple_hits, genetic_code = genetic_code))
  # Testing every branch needs no foreground, so a prepared folder's foreground is ignored.
  .phylo_run("absrel", alignment_dir, species_tree_file, if (branches == "all" && is.null(prepared_dir)) NULL else foreground_file,
             gene_tree_dir, output_dir, run_name, settings, seed, hyphy_cpu, threads, min_taxa, terminal_stops, overwrite, hyphy,
             gene_tree_mode = gene_tree_mode, prepared_dir = prepared_dir, all_branches = branches == "all")
}
