#' Test for episodic positive selection on foreground branches with HyPhy BUSTED
#'
#' Runs HyPhy's BUSTED (Murrell et al. 2015) on each orthogroup codon alignment, with the
#' foreground branches labelled from a species tree, and corrects across orthogroups.
#' BUSTED asks whether at least one site on at least one foreground branch evolved with
#' dN/dS above 1. It is a likelihood-ratio test across several sequences, and is a separate
#' analysis from pairwise dN/dS: a significant BUSTED result is not the same as, and cannot
#' be compared with, a pairwise dN/dS above 1.
#'
#' Experimental (new in v0.5.0). It reproduces HyPhy run directly. Its false-positive rate
#' and power were calibrated only in simulation, on one five-taxon tree with one
#' foreground, which recommends \code{multiple_hits = "Double"} and three starts
#' (\code{extra_seeds = 2:3}). It is not calibrated for arbitrary phylogenies or real data.
#' See `docs/phylo-methods.md`.
#'
#' @section Inputs:
#' Either \code{prepared_dir}, the output folder of [phylo_prepare()], or all of:
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
#'   \item \code{gene_tree_dir} (optional): \code{<orthogroup>.nwk} topologies, with tips
#'     named by taxon, sequence name, or gene ID.
#' }
#' With \code{prepared_dir}, every orthogroup phylo_prepare excluded is also listed, and the
#' labelled tree of each orthogroup must equal the one phylo_prepare wrote.
#'
#' @section Outcomes:
#' Each orthogroup is \code{tested}, \code{excluded} (with \code{reason} and
#' \code{exclusion_basis}: \code{mathematical}, \code{no_data}, \code{policy},
#' \code{data_quality}, or \code{input}), or \code{failed} (a step errored; basis
#' \code{software}). \code{stage} says whether this happened in preparation or analysis.
#'
#' @section Missing taxa and identifiability:
#' The species tree is pruned to the taxa in each alignment. Branches that join when taxa
#' are missing become one branch. A requested branch joined to a background branch cannot
#' be separated from it: if none can be separated, the orthogroup is excluded as
#' \code{foreground_unidentifiable} (mathematical); if some can, testing only those would
#' narrow the hypothesis, so it is excluded as \code{foreground_partly_unidentifiable}
#' (policy). A hypothesis row with no taxon present, or a terminal row lacking a taxon, is
#' \code{foreground_incomplete} (policy). Nothing is ever relabelled.
#'
#' @section Gene trees:
#' With \code{gene_tree_mode = "topology"} the gene tree is given to HyPhy; every requested
#' terminal and stem branch must exist in it (\code{foreground_not_in_gene_tree}
#' otherwise). With \code{"report"} the species tree is used. In both modes discordance with
#' the pruned species tree is reported (\code{gt_*} columns, \code{gene_tree_discordant} and
#' \code{gene_tree_discordant_near_foreground} flags), and never excludes an orthogroup.
#'
#' @section Statistics:
#' The p-value is HyPhy's, from a 50:50 mixture of chi-square distributions with 0 and 2
#' degrees of freedom, so it is exactly 0.5 when the likelihood ratio is 0. \code{p_adj}
#' applies Benjamini-Hochberg across one family: every orthogroup with status
#' \code{tested} in this run. Raw p-values are kept.
#'
#' @section Reproducibility and multiple starts:
#' HyPhy draws random starting values, and its result depends on them and on the thread
#' count. dndsR always sets the seed (\code{seed}) and threads per HyPhy run
#' (\code{hyphy_cpu}), and records both. \code{extra_seeds} runs further independent starts
#' and reports their spread and the best likelihoods found (\code{lrt_best_fits},
#' \code{p_best_fits}), and \code{p_best_fits_adj}, Benjamini-Hochberg over
#' \code{p_best_fits} across the same tested orthogroups as \code{p_adj}. The stage 2
#' calibration recommends three starts and reading \code{p_best_fits} and
#' \code{p_best_fits_adj} (\code{docs/phylo-methods.md}). \code{p_value} and \code{p_adj}
#' always come from \code{seed} alone, as before. The largest likelihood ratio from a single
#' start is never chosen.
#'
#' @section Outputs:
#' In \code{file.path(output_dir, run_name)}: \code{<run_name>_results.tsv} (one row per
#' orthogroup), \code{_branches.tsv}, \code{_sequences.tsv}, \code{_provenance.json}, and
#' \code{hyphy/<orthogroup>/} with the exact HyPhy inputs, logs, and JSON (extra starts in
#' \code{start_seed<k>/}).
#'
#' @param prepared_dir Output folder of [phylo_prepare()]. Replaces the next three.
#' @param alignment_dir Folder of codon alignments, one per orthogroup.
#' @param species_tree_file Rooted Newick species tree.
#' @param foreground_file Foreground hypothesis (see Inputs).
#' @param gene_tree_dir Optional folder of gene trees, one per orthogroup.
#' @param gene_tree_mode "topology" or "report" (see Gene trees). Defaults to phylo_prepare's
#'   choice, or "topology".
#' @param output_dir Folder for results.
#' @param run_name Name of this analysis, used for its folder and files.
#' @param srv Synonymous rate variation: "Yes" (HyPhy's default), "No", or "Branch-site".
#' @param multiple_hits "None" (HyPhy's default), "Double", or "Double+Triple".
#' @param rates Number of omega rate classes (HyPhy's default 3).
#' @param syn_rates Number of synonymous rate classes when \code{srv} is not "No" (default 3).
#' @param error_sink "No" (HyPhy's default) or "Yes" (BUSTED-E). HyPhy describes the error
#'   sink as an advanced experimental setting, and so does dndsR.
#' @param grid_size,starting_points HyPhy's optimisation settings (defaults 250 and 1).
#' @param kill_zero_lengths "Yes" (HyPhy's default) deletes internal branches estimated at
#'   zero length before fitting codon models. If that deletes a requested foreground branch,
#'   the orthogroup is excluded as \code{foreground_branch_zero_length} (policy), because
#'   only a narrower hypothesis would be tested. "No" keeps every branch.
#' @param genetic_code Only "Universal" is supported.
#' @param seed Random seed for HyPhy (default 1). Determines the reported result.
#' @param extra_seeds Further seeds for independent starts, reported as diagnostics.
#' @param hyphy_cpu Threads for each HyPhy run (default 1). Changing it changes results
#'   slightly, as a different seed does.
#' @param threads Orthogroups to run at once.
#' @param min_taxa Fewest sequences an orthogroup may have (default 4).
#' @param terminal_stops "mask" (default) replaces a stop codon ending a sequence with
#'   gaps, and records it; "exclude" excludes the orthogroup instead. Internal stop codons
#'   always exclude the orthogroup.
#' @param overwrite If \code{FALSE} (default), an orthogroup whose HyPhy inputs, settings,
#'   seeds, and version match a finished earlier run is not rerun.
#' @param hyphy Path to the HyPhy executable, if not on \code{PATH}.
#' @return Invisibly, the results table.
#' @references Murrell B, et al. (2015) Gene-wide identification of episodic selection.
#'   Mol Biol Evol 32:1365-1371. Wisotsky SR, et al. (2020) Synonymous site-to-site
#'   substitution rate variation dramatically inflates false positive rates of selection
#'   analyses. Mol Biol Evol 37:2430-2439.
#' @seealso [phylo_prepare()], [phylo_example()], [verify_phylo_example()]
#' @export
phylo_busted <- function(prepared_dir = NULL,
                         alignment_dir = NULL,
                         species_tree_file = NULL,
                         foreground_file = NULL,
                         gene_tree_dir = NULL,
                         gene_tree_mode = NULL,
                         output_dir = ".",
                         run_name = "busted",
                         srv = "Yes",
                         multiple_hits = "None",
                         rates = 3L,
                         syn_rates = 3L,
                         error_sink = "No",
                         grid_size = 250L,
                         starting_points = 1L,
                         kill_zero_lengths = "Yes",
                         genetic_code = "Universal",
                         seed = 1L,
                         extra_seeds = NULL,
                         hyphy_cpu = 1L,
                         threads = 1L,
                         min_taxa = 4L,
                         terminal_stops = c("mask", "exclude"),
                         overwrite = FALSE,
                         hyphy = NULL) {
  terminal_stops <- match.arg(terminal_stops)
  settings <- .hyphy_settings("busted", list(
    srv = srv, multiple_hits = multiple_hits, rates = rates, syn_rates = syn_rates, error_sink = error_sink,
    grid_size = grid_size, starting_points = starting_points, kill_zero_lengths = kill_zero_lengths, genetic_code = genetic_code
  ))
  if (identical(settings$error_sink, "Yes")) message("Note: error_sink = 'Yes' (BUSTED-E) is experimental in HyPhy and not validated in dndsR.")
  .phylo_run("busted", alignment_dir, species_tree_file, foreground_file, gene_tree_dir, output_dir, run_name,
             settings, seed, hyphy_cpu, threads, min_taxa, terminal_stops, overwrite, hyphy,
             extra_seeds = extra_seeds, gene_tree_mode = gene_tree_mode, prepared_dir = prepared_dir)
}
