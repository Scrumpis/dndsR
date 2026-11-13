#' Test for episodic positive selection on foreground branches with HyPhy BUSTED
#'
#' Runs HyPhy's BUSTED (Murrell et al. 2015) on each orthogroup codon alignment, with the
#' foreground branches labelled from a species tree, and corrects across orthogroups.
#' BUSTED asks whether at least one site on at least one foreground branch evolved with
#' dN/dS above 1. It is a likelihood-ratio test across several sequences, and is a separate
#' analysis from pairwise dN/dS: a significant BUSTED result is not the same as, and cannot
#' be compared with, a pairwise dN/dS above 1.
#'
#' This is phase 1 of the phylogenetic module: it needs codon alignments you have already
#' made. It wraps HyPhy and has been checked to reproduce direct HyPhy runs; it has not yet
#' been calibrated as a selection pipeline. See `docs/phylo-methods.md`.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{alignment_dir}: one in-frame codon alignment per orthogroup, in FASTA,
#'     named \code{<orthogroup>.fasta} (or \code{.fa}, \code{.fas}, \code{.fna}).
#'     Sequence names are \code{taxon} or \code{taxon|gene_id}. Only single-copy
#'     orthogroups are tested: one sequence per taxon.
#'   \item \code{species_tree_file}: a rooted Newick tree of the taxa. Branch lengths and
#'     node labels are ignored. HyPhy receives the topology only, because it uses given
#'     lengths as starting values, and species-tree lengths are not in its units.
#'   \item \code{foreground_file}: tab-separated, with columns \code{scope} and
#'     \code{taxa}. Scope \code{terminal} labels the terminal branch of each listed taxon;
#'     \code{stem} labels the branch leading to the clade of the listed taxa; \code{clade}
#'     labels that stem and every branch within the clade. Rows are combined.
#'   \item \code{gene_tree_dir} (optional): \code{<orthogroup>.nwk} topologies to use
#'     instead of the pruned species tree, with tips named by taxon or sequence name.
#' }
#'
#' @section Missing taxa and identifiability:
#' The species tree is pruned to the taxa in each alignment. Branches that join when taxa
#' are missing become one branch. If such a branch would combine foreground and background
#' branches, the hypothesis cannot be tested for that orthogroup, and it is skipped with
#' reason \code{foreground_not_identifiable}; it is never relabelled. For example, a stem
#' hypothesis is not identifiable when a clade member, or the clade's sister lineage, is
#' missing. With a gene tree, every foreground branch must exist in it, or the orthogroup
#' is skipped as \code{foreground_discordant}.
#'
#' @section Statistics:
#' The p-value is HyPhy's, from a 50:50 mixture of chi-square distributions with 0 and 2
#' degrees of freedom, so it is exactly 0.5 when the likelihood ratio is 0. \code{p_adj}
#' applies Benjamini-Hochberg across one family: every orthogroup with status
#' \code{tested} in this run. Skipped and failed orthogroups are not in the family and are
#' reported with reasons. Raw p-values are kept.
#'
#' @section Reproducibility:
#' HyPhy draws random starting values, and its result depends on them and on the thread
#' count. dndsR always sets the seed (\code{seed}) and threads per HyPhy run
#' (\code{hyphy_cpu}), and records both. With the same HyPhy build, these give identical
#' results. A different seed can move the likelihood ratio of a gene; see
#' `docs/validation/phylo-busted-wrapper.md`.
#'
#' @section Outputs:
#' In \code{file.path(output_dir, run_name)}:
#' \itemize{
#'   \item \code{<run_name>_results.tsv}: one row per orthogroup, tested or not, with its status and reason.
#'   \item \code{<run_name>_branches.tsv}: every branch given to HyPhy, its split, and
#'     whether it is foreground.
#'   \item \code{<run_name>_sequences.tsv}: per-sequence QC.
#'   \item \code{<run_name>_provenance.json}: versions, settings, seed, inputs, counts.
#'   \item \code{hyphy/<orthogroup>/}: the exact alignment, tree, log, and JSON of each run.
#' }
#'
#' @param alignment_dir Folder of codon alignments, one per orthogroup.
#' @param species_tree_file Rooted Newick species tree.
#' @param foreground_file Foreground hypothesis (see Inputs).
#' @param gene_tree_dir Optional folder of gene trees, one per orthogroup.
#' @param output_dir Folder for results.
#' @param run_name Name of this analysis, used for its folder and files.
#' @param srv Synonymous rate variation: "Yes" (HyPhy's default), "No", or "Branch-site".
#' @param multiple_hits "None" (HyPhy's default), "Double", or "Double+Triple".
#' @param rates Number of omega rate classes (HyPhy's default 3).
#' @param syn_rates Number of synonymous rate classes when \code{srv} is not "No" (default 3).
#' @param error_sink "No" (HyPhy's default) or "Yes". HyPhy describes the error sink as an
#'   advanced experimental setting.
#' @param grid_size,starting_points HyPhy's optimisation settings (defaults 250 and 1).
#' @param genetic_code Only "Universal" is supported.
#' @param seed Random seed for HyPhy (default 1).
#' @param hyphy_cpu Threads for each HyPhy run (default 1). Changing it changes results
#'   slightly, as a different seed does.
#' @param threads Orthogroups to run at once.
#' @param min_taxa Fewest sequences an orthogroup may have (default 4).
#' @param terminal_stops "mask" (default) replaces a stop codon ending a sequence with
#'   gaps, and records it; "exclude" skips the orthogroup instead. Internal stop codons
#'   always exclude the orthogroup.
#' @param overwrite If \code{FALSE} (default), an orthogroup whose HyPhy inputs, settings,
#'   and version match a finished earlier run is not rerun.
#' @param hyphy Path to the HyPhy executable, if not on \code{PATH}.
#' @return Invisibly, the results table.
#' @references Murrell B, et al. (2015) Gene-wide identification of episodic selection.
#'   Mol Biol Evol 32:1365-1371. Wisotsky SR, et al. (2020) Synonymous site-to-site
#'   substitution rate variation dramatically inflates false positive rates of selection
#'   analyses. Mol Biol Evol 37:2430-2439.
#' @seealso [phylo_example()], [verify_phylo_example()]
#' @export
phylo_busted <- function(alignment_dir,
                         species_tree_file,
                         foreground_file,
                         gene_tree_dir = NULL,
                         output_dir = ".",
                         run_name = "busted",
                         srv = "Yes",
                         multiple_hits = "None",
                         rates = 3L,
                         syn_rates = 3L,
                         error_sink = "No",
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
  settings <- .hyphy_settings("busted", list(
    srv = srv, multiple_hits = multiple_hits, rates = rates, syn_rates = syn_rates, error_sink = error_sink,
    grid_size = grid_size, starting_points = starting_points, genetic_code = genetic_code
  ))
  .phylo_run("busted", alignment_dir, species_tree_file, foreground_file, gene_tree_dir, output_dir, run_name,
             settings, seed, hyphy_cpu, threads, min_taxa, terminal_stops, overwrite, hyphy)
}
