#' Lineage-level summaries of phylogenetic selection analyses (experimental)
#'
#' Aggregates per-orthogroup results of [phylo_branch_omega()], [phylo_absrel()], and
#' [phylo_relax()] into one row per species-tree branch, with coverage, diagnostics, and
#' matched comparisons between lineages, and writes stable branch identifiers for joining
#' lineage traits measured elsewhere (LTR activity, genome-size change). Experimental.
#'
#' @section What the summaries are, and are not:
#' Everything here is descriptive. A lineage's pooled omega summarizes the genes tested on
#' that branch; it is not a test of adaptation, and a difference between lineages is not
#' evidence for one without a model of what else differs between them (branch length,
#' which genes were sampled, saturation, alignment quality). Formal tests stay per gene,
#' in the analyses summarized (`p_adj` columns there; `p_bh_lineage` for aBSREL).
#' Lineage-level hypotheses about traits belong in a phylogenetic comparative model fitted
#' with established tools (for example `phylolm`, `caper`, `nlme::gls` with
#' `ape::corBrownian`), using the tables and tree written here.
#'
#' @section Pooled omega:
#' For each species-tree branch, the free-ratio fits of every orthogroup in which the
#' branch exists, unmerged by pruning and with dS at most 1, are pooled as
#' `(sum N*dN / sum N) / (sum S*dS / sum S)`: substitutions per site summed over genes.
#' This weights genes by their information, so genes with few synonymous changes, whose
#' ratios are unstable or undefined, do not dominate, as they do in a mean of per-gene
#' ratios. A percentile interval comes from resampling orthogroups (`bootstrap` times).
#' The median and quartiles of per-gene omega, over reliable estimates only
#' (`reliable_free`), are given beside it.
#'
#' @section Matched comparisons:
#' Two lineages estimated from different genes differ partly because the genes differ.
#' `_lineage_pairs.tsv` compares pairs of branches (terminal branches by default) on the
#' orthogroups where both have reliable estimates: pooled omega of each on that shared set,
#' the median per-gene log2 ratio, and a sign test across genes, Benjamini-Hochberg
#' corrected across pairs. Genes are not independent, so treat these p-values as
#' descriptive.
#'
#' @section Outputs:
#' In `file.path(output_dir, run_name)`:
#' \itemize{
#'   \item `_branch_ids.tsv`: one row per species-tree branch: `branch_id` (taxon name
#'     for terminal branches; node label or `mrca(x,y)` otherwise), `terminal`,
#'     `ape_node` (the node below the branch, numbered as in `ape`, so ancestral states
#'     from `ape::ace` or `phytools` join directly), `parent_id`, descendant taxa, and the
#'     species tree's branch length.
#'   \item `_species_tree_ids.nwk`: the species tree with internal nodes labelled by
#'     `tree_label`, the Newick-safe form of `branch_id` (equal to it when the species
#'     tree's internal nodes carry plain labels; labelling them is recommended).
#'   \item `_lineages.tsv`: one row per branch: coverage, pooled omega with interval,
#'     per-gene omega quantiles, dS and branch-length diagnostics, flag rates, and aBSREL
#'     counts when given.
#'   \item `_lineage_pairs.tsv`: matched comparisons (see above).
#'   \item `_orthogroup_branches.tsv`: one row per orthogroup and species-tree branch,
#'     the long table for downstream models.
#'   \item `_foreground.tsv` and `_relax.tsv`: summaries of the two-ratio test and of
#'     RELAX, when available.
#'   \item `_provenance.json`.
#' }
#'
#' @param branch_omega_dir Output folder of [phylo_branch_omega()] (`output_dir/run_name`).
#' @param absrel_dir Output folder of [phylo_absrel()], optional.
#' @param relax_dir Output folder of [phylo_relax()], optional.
#' @param output_dir Folder for the summary.
#' @param run_name Name of the summary, used for its folder and files.
#' @param pairs Which branches to compare in matched pairs: "terminal" (default),
#'   "all", or "none".
#' @param alpha Significance level used for counts of significant tests.
#' @param bootstrap Resamples of orthogroups for the pooled-omega interval (0 for none).
#' @param seed Seed for the bootstrap.
#' @return Invisibly, the lineage table.
#' @seealso [phylo_branch_omega()], [phylo_absrel()], [phylo_relax()]
#' @export
phylo_summarize <- function(branch_omega_dir = NULL,
                            absrel_dir = NULL,
                            relax_dir = NULL,
                            output_dir = ".",
                            run_name = "summary",
                            pairs = c("terminal", "all", "none"),
                            alpha = 0.05,
                            bootstrap = 1000L,
                            seed = 1L) {
  .need_ape()
  pairs <- match.arg(pairs)
  if (is.null(branch_omega_dir) && is.null(absrel_dir) && is.null(relax_dir)) {
    stop("Give at least one of branch_omega_dir, absrel_dir, relax_dir.", call. = FALSE)
  }
  if (!grepl("^[A-Za-z0-9_.-]+$", run_name)) stop("run_name may contain only letters, digits, '.', '_' and '-'.", call. = FALSE)
  runs <- list(branch_omega = .phylo_read_run(branch_omega_dir, "branch_omega"),
               absrel = .phylo_read_run(absrel_dir, "absrel"), relax = .phylo_read_run(relax_dir, "relax"))
  runs <- Filter(Negate(is.null), runs)
  trees <- unique(vapply(runs, function(r) r$prov$species_tree$md5, character(1)))
  if (length(trees) > 1L) stop("The runs used different species trees; summarize runs of one species tree together.", call. = FALSE)
  tree <- .phylo_read_species_tree(runs[[1]]$prov$species_tree$path)
  if (!identical(unname(tools::md5sum(runs[[1]]$prov$species_tree$path)), trees)) {
    stop("The species tree changed since the analyses ran: ", runs[[1]]$prov$species_tree$path, call. = FALSE)
  }
  ids <- .phylo_branch_ids(tree, runs[[1]]$prov$species_tree$path)
  run_dir <- file.path(output_dir, run_name)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- function(what, ext = "tsv") file.path(run_dir, sprintf("%s_%s.%s", run_name, what, ext))
  write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_tsv(ids, out_file("branch_ids"))
  labelled <- tree
  labelled$node.label <- ids$tree_label[match(ape::Ntip(tree) + seq_len(tree$Nnode), ids$ape_node)]
  labelled$node.label[is.na(labelled$node.label)] <- "root"
  ape::write.tree(labelled, out_file("species_tree_ids", "nwk"))

  lineages <- ids[, c("branch_id", "terminal", "n_descendants", "species_tree_length")]
  long <- NULL
  if (!is.null(runs$branch_omega)) {
    bo <- .phylo_summarize_omega(runs$branch_omega, ids, bootstrap, seed)
    lineages <- merge(lineages, bo$lineages, by = "branch_id", all.x = TRUE, sort = FALSE)
    long <- bo$long
    if (pairs != "none") {
      keep <- if (pairs == "terminal") ids$branch_id[ids$terminal] else ids$branch_id
      write_tsv(.phylo_lineage_pairs(long, keep), out_file("lineage_pairs"))
    }
    write_tsv(long, out_file("orthogroup_branches"))
    fg <- .phylo_summarize_foreground(runs$branch_omega, alpha)
    if (!is.null(fg)) write_tsv(fg, out_file("foreground"))
  }
  if (!is.null(runs$absrel)) {
    lineages <- merge(lineages, .phylo_summarize_absrel(runs$absrel, ids, alpha), by = "branch_id", all.x = TRUE, sort = FALSE)
  }
  if (!is.null(runs$relax)) write_tsv(.phylo_summarize_relax(runs$relax, alpha), out_file("relax"))
  lineages <- lineages[match(ids$branch_id, lineages$branch_id), , drop = FALSE]
  write_tsv(lineages, out_file("lineages"))

  version_of <- function(pkg) if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  prov <- list(
    analysis = "summarize", dndsR = version_of("dndsR"), inputs = lapply(runs, function(r) list(dir = r$dir, results_md5 = r$results_md5)),
    species_tree = runs[[1]]$prov$species_tree, pairs = pairs, alpha = alpha, bootstrap = bootstrap, seed = seed,
    pooled_omega = "(sum N*dN / sum N) / (sum S*dS / sum S) over orthogroups where the branch exists unmerged and dS <= 1, free-ratio fits",
    interpretation = paste("Descriptive summaries of per-gene estimates. Not evidence of adaptation; lineage comparisons confound",
                           "branch length, gene sampling, and saturation. Fit lineage-level hypotheses in a phylogenetic comparative model.")
  )
  writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), out_file("provenance", "json"))
  message(sprintf("[summarize] %d species-tree branches. Results: %s", nrow(lineages), out_file("lineages")))
  invisible(lineages)
}

#' Read one analysis folder: results, branches, and provenance
#' @keywords internal
.phylo_read_run <- function(dir, analysis) {
  if (is.null(dir)) return(NULL)
  prov_file <- list.files(dir, pattern = "_provenance\\.json$", full.names = TRUE)
  if (length(prov_file) != 1L) stop("Expected one *_provenance.json in ", dir, call. = FALSE)
  prov <- jsonlite::read_json(prov_file)
  if (!identical(prov$analysis, analysis)) stop(dir, " holds a ", prov$analysis, " analysis, not ", analysis, ".", call. = FALSE)
  base <- sub("_provenance\\.json$", "", prov_file)
  read <- function(what) utils::read.delim(paste0(base, "_", what, ".tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  list(dir = dir, prov = prov, results = read("results"), branches = read("branches"),
       results_md5 = unname(tools::md5sum(paste0(base, "_results.tsv"))))
}

#' A branch name usable as a Newick label: anything but letters, digits, '.', '_' and '-' becomes '_'
#' @keywords internal
.phylo_tree_label <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9_.-]+", "_", x))

#' Stable identifiers for every branch of the species tree
#'
#' @param tree Species tree from [.phylo_read_species_tree()].
#' @param path Its file, read again for branch lengths, matched by descendant taxa.
#' @keywords internal
.phylo_branch_ids <- function(tree, path = NULL) {
  n <- ape::Ntip(tree)
  names_all <- c(tree$tip.label, .phylo_internal_names(tree))
  desc <- .phylo_descendants(tree)
  child <- tree$edge[, 2]; parent <- tree$edge[, 1]
  data.frame(
    branch_id = names_all[child],
    # Newick cannot hold brackets or commas in a label, so the exported tree uses this
    # form; it equals branch_id when the species tree's own node labels are plain.
    tree_label = .phylo_tree_label(names_all[child]),
    terminal = child <= n,
    ape_node = child,
    parent_id = ifelse(parent == n + 1L, NA_character_, names_all[parent]),
    descendant_taxa = vapply(desc[child], function(d) paste(sort(d, method = "radix"), collapse = ","), character(1)),
    n_descendants = lengths(desc[child]),
    species_tree_length = .phylo_edge_lengths(tree, path),
    stringsAsFactors = FALSE
  )
}

#' Branch lengths of the species tree file, in the edge order of `tree`
#' @keywords internal
.phylo_edge_lengths <- function(tree, path) {
  out <- rep(NA_real_, nrow(tree$edge))
  if (is.null(path)) return(out)
  raw <- ape::read.tree(path)
  if (is.null(raw$edge.length)) return(out)
  key <- function(tr) vapply(.phylo_descendants(tr)[tr$edge[, 2]], function(d) paste(sort(d, method = "radix"), collapse = ","), character(1))
  raw_keys <- key(raw)
  raw_keys[duplicated(raw_keys)] <- NA_character_
  raw$edge.length[match(key(tree), raw_keys)]
}

#' Pooled omega: substitutions per site, summed over orthogroups
#' @keywords internal
.phylo_pooled_omega <- function(NdN, N, SdS, S) {
  ok <- is.finite(NdN) & is.finite(N) & is.finite(SdS) & is.finite(S)
  if (!any(ok) || sum(SdS[ok]) <= 0) return(NA_real_)
  (sum(NdN[ok]) / sum(N[ok])) / (sum(SdS[ok]) / sum(S[ok]))
}

#' Lineage summaries of free-ratio branch omega
#' @keywords internal
.phylo_summarize_omega <- function(run, ids, bootstrap, seed) {
  res <- run$results; b <- run$branches
  tested <- res$orthogroup[res$status == "tested"]
  b <- b[b$orthogroup %in% tested, , drop = FALSE]
  if (!"omega_free" %in% names(b)) stop("The branch-omega run has no free-ratio fits; rerun with models including 'free'.", call. = FALSE)
  merged <- grepl("+", b$species_tree_branches, fixed = TRUE)
  b$branch_id <- ifelse(merged | !b$in_species_tree, NA_character_, b$species_tree_branches)
  long <- data.frame(orthogroup = b$orthogroup, branch_id = b$branch_id, species_tree_branches = b$species_tree_branches,
                     t = b$t_free, N = b$N_free, S = b$S_free, dN = b$dN_free, dS = b$dS_free, NdN = b$NdN_free, SdS = b$SdS_free,
                     omega = b$omega_free, flags = b$flags, reliable = b$reliable_free, stringsAsFactors = FALSE)
  starts <- res$qc_flags[match(long$orthogroup, res$orthogroup)]
  long$starts_disagree <- grepl("starts_disagree:free", starts, fixed = TRUE)
  usable <- long[!is.na(long$branch_id), , drop = FALSE]
  q <- function(x, p) if (length(x)) unname(stats::quantile(x, p, na.rm = TRUE)) else NA_real_
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old_seed)) rm(list = ".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)
  rows <- lapply(ids$branch_id, function(id) {
    x <- usable[usable$branch_id == id, , drop = FALSE]
    pool <- x[!grepl("high_ds", x$flags), , drop = FALSE]
    rel <- x[x$reliable %in% TRUE, , drop = FALSE]
    pooled <- .phylo_pooled_omega(pool$NdN, pool$N, pool$SdS, pool$S)
    ci <- c(NA_real_, NA_real_)
    if (bootstrap > 0 && nrow(pool) >= 5) {
      boots <- replicate(bootstrap, { i <- sample.int(nrow(pool), replace = TRUE); .phylo_pooled_omega(pool$NdN[i], pool$N[i], pool$SdS[i], pool$S[i]) })
      ci <- unname(stats::quantile(boots, c(0.025, 0.975), na.rm = TRUE))
    }
    has <- function(f) if (nrow(x)) mean(grepl(f, x$flags)) else NA_real_
    data.frame(
      branch_id = id,
      n_orthogroups_tested = length(tested),
      n_estimated = nrow(x),
      n_merged_or_absent = length(tested) - nrow(x),
      n_reliable = nrow(rel),
      n_pooled = nrow(pool),
      pooled_omega = pooled, pooled_omega_lo = ci[1], pooled_omega_hi = ci[2],
      median_omega_reliable = q(rel$omega, 0.5), q25_omega_reliable = q(rel$omega, 0.25), q75_omega_reliable = q(rel$omega, 0.75),
      frac_reliable_omega_gt1 = if (nrow(rel)) mean(rel$omega > 1) else NA_real_,
      median_ds = q(x$dS, 0.5), median_branch_length = q(x$t, 0.5), sum_syn_changes = sum(x$SdS, na.rm = TRUE),
      frac_high_ds = has("high_ds"), frac_few_synonymous = has("few_synonymous_changes"),
      frac_omega_upper_bound = has("omega_at_upper_bound"), frac_zero_length = has("zero_length"),
      frac_starts_disagree = if (nrow(x)) mean(x$starts_disagree) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  list(lineages = do.call(rbind, rows), long = long)
}

#' Matched comparisons between pairs of lineages, on shared orthogroups
#' @keywords internal
.phylo_lineage_pairs <- function(long, keep) {
  rel <- long[long$reliable %in% TRUE & long$branch_id %in% keep, , drop = FALSE]
  ids <- sort(unique(rel$branch_id), method = "radix")
  if (length(ids) < 2L) return(data.frame())
  out <- list()
  for (a in seq_along(ids)) for (b in seq_along(ids)) if (a < b) {
    x <- rel[rel$branch_id == ids[a], , drop = FALSE]; y <- rel[rel$branch_id == ids[b], , drop = FALSE]
    og <- intersect(x$orthogroup, y$orthogroup)
    x <- x[match(og, x$orthogroup), , drop = FALSE]; y <- y[match(og, y$orthogroup), , drop = FALSE]
    pos <- x$omega > 0 & y$omega > 0
    lr <- log2(x$omega[pos] / y$omega[pos])
    n_gt <- sum(lr > 0); n_lt <- sum(lr < 0)
    out[[length(out) + 1L]] <- data.frame(
      branch_a = ids[a], branch_b = ids[b], n_shared = length(og),
      pooled_omega_a = .phylo_pooled_omega(x$NdN, x$N, x$SdS, x$S),
      pooled_omega_b = .phylo_pooled_omega(y$NdN, y$N, y$SdS, y$S),
      median_log2_ratio = if (length(lr)) stats::median(lr) else NA_real_,
      n_a_higher = n_gt, n_b_higher = n_lt,
      sign_test_p = if (n_gt + n_lt) stats::binom.test(n_gt, n_gt + n_lt)$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, out)
  out$sign_test_p_adj <- stats::p.adjust(out$sign_test_p, method = "BH")
  out
}

#' Summary of the two-ratio test across orthogroups
#' @keywords internal
.phylo_summarize_foreground <- function(run, alpha) {
  r <- run$results
  if (is.null(r$two_ratio_status)) return(NULL)
  t <- r[r$status == "tested" & r$two_ratio_status %in% "tested" & !is.na(r$omega_fg), , drop = FALSE]
  if (!nrow(t)) return(NULL)
  b <- run$branches
  b <- b[b$orthogroup %in% t$orthogroup & !is.na(b$NdN_two), , drop = FALSE]
  fgb <- b$foreground %in% TRUE
  pos <- t$omega_fg > 0 & t$omega_bg > 0
  data.frame(
    n_tested = nrow(t),
    n_fg_higher_sig = sum(t$p_adj_two_vs_m0 <= alpha & t$omega_fg > t$omega_bg, na.rm = TRUE),
    n_fg_lower_sig = sum(t$p_adj_two_vs_m0 <= alpha & t$omega_fg < t$omega_bg, na.rm = TRUE),
    n_fg_omega_gt1_sig = sum(t$p_adj_fg_omega_gt1 <= alpha, na.rm = TRUE),
    median_log2_fg_over_bg = if (any(pos)) stats::median(log2(t$omega_fg[pos] / t$omega_bg[pos])) else NA_real_,
    pooled_omega_fg = .phylo_pooled_omega(b$NdN_two[fgb], b$N_two[fgb], b$SdS_two[fgb], b$S_two[fgb]),
    pooled_omega_bg = .phylo_pooled_omega(b$NdN_two[!fgb], b$N_two[!fgb], b$SdS_two[!fgb], b$S_two[!fgb]),
    alpha = alpha, stringsAsFactors = FALSE
  )
}

#' aBSREL counts per species-tree branch
#' @keywords internal
.phylo_summarize_absrel <- function(run, ids, alpha) {
  b <- run$branches
  tested <- run$results$orthogroup[run$results$status == "tested"]
  b <- b[b$orthogroup %in% tested & !is.na(b$p_raw) & !grepl("+", b$species_tree_branches, fixed = TRUE), , drop = FALSE]
  do.call(rbind, lapply(ids$branch_id, function(id) {
    x <- b[b$species_tree_branches %in% id, , drop = FALSE]
    data.frame(branch_id = id, absrel_n_tested = nrow(x),
               absrel_n_sig_lineage = sum(x$p_bh_lineage <= alpha, na.rm = TRUE),
               absrel_n_sig_holm = sum(x$p_holm <= alpha, na.rm = TRUE),
               absrel_median_baseline_omega = if (nrow(x)) stats::median(x$baseline_omega, na.rm = TRUE) else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

#' RELAX summary across orthogroups
#' @keywords internal
.phylo_summarize_relax <- function(run, alpha) {
  r <- run$results
  t <- r[r$status == "tested" & !is.na(r$k), , drop = FALSE]
  data.frame(
    n_tested = nrow(t),
    n_relaxation_sig = sum(t$p_adj <= alpha & t$k < 1, na.rm = TRUE),
    n_intensification_sig = sum(t$p_adj <= alpha & t$k > 1, na.rm = TRUE),
    median_log2_k = if (nrow(t)) stats::median(log2(pmax(t$k, 1e-8))) else NA_real_,
    alpha = alpha,
    note = "RELAX k scales omega classes of test relative to reference branches; it is not a test of positive selection.",
    stringsAsFactors = FALSE
  )
}
