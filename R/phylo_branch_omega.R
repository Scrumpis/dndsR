#' Branch-specific dN/dS with PAML codeml branch models (experimental)
#'
#' Fits codeml branch models to every orthogroup of a [phylo_prepare()] folder or a folder
#' of codon alignments: one omega for the whole tree (`M0`), one omega per branch
#' (free-ratio, `model = 1`), and one omega for foreground and one for background branches
#' (two-ratio, `model = 2`), with likelihood-ratio tests between them. Experimental: dndsR
#' reproduces codeml run directly, but its estimates and tests with dndsR's inputs have been
#' checked only in simulation (`docs/phylo-methods.md`).
#'
#' @section Models and tests:
#' Branch models give each branch, or each set of branches, one omega averaged over every
#' site of the gene. They are not branch-site models: a branch with episodic positive
#' selection at a few sites usually has omega well below 1, and an omega above 1 on a
#' branch is not by itself evidence of positive selection (short branches, few synonymous
#' changes, and alignment or orthology errors can produce it). Tests, each corrected by
#' Benjamini-Hochberg across tested orthogroups:
#' \itemize{
#'   \item `two_ratio` against `M0` (`p_two_vs_m0`, 1 df): does the foreground omega
#'     differ from the background omega?
#'   \item `free` against `M0` (`p_free_vs_m0`, branches minus 1 df): does omega vary among
#'     branches? Free-ratio estimates are exploratory: with many parameters, each branch's
#'     omega is noisy.
#'   \item `two_ratio` against `two_ratio_null`, where the foreground omega is fixed at 1
#'     (`p_fg_omega_gt1`): is the foreground omega, averaged over the gene, above 1? One
#'     sided, from the signed likelihood root.
#' }
#' The models assume one omega per branch (or set) for all sites, F3x4 codon frequencies by
#' default, one transition/transversion ratio, and the tree topology given.
#'
#' @section Optimization:
#' codeml can stop at a local optimum. Each model is fitted from each value in
#' `omega_starts` (default 0.4, 0.1, and 2), and the best likelihood is kept. The results
#' report the spread of likelihoods across starts and how many starts reached the best;
#' disagreement is flagged. A negative likelihood ratio (a nested model fitting better) is
#' set to 0 and flagged.
#'
#' @section Branches:
#' Branch identity follows the species tree (`species_tree_branches`); a branch of the
#' pruned, unrooted tree that joins several species-tree branches (for example the two
#' branches at the root) is named by all of them and flagged `merged_branch`. Per-branch
#' flags mark estimates to treat with care: `zero_length`, `few_synonymous_changes`
#' (fewer than one expected synonymous substitution, `S*dS < 1`), `omega_at_upper_bound`
#' (codeml's 999), and `high_ds` (dS above 1, where saturation makes estimates unreliable).
#' `reliable_free` combines them for the free-ratio estimates.
#'
#' @section Outputs:
#' In `file.path(output_dir, run_name)`: `_results.tsv` (one row per orthogroup: models'
#' likelihoods, parameters, tests, and optimizer diagnostics), `_branches.tsv` (one row per
#' orthogroup and branch: branch lengths, dN, dS, omega, and substitution counts under
#' each model, with flags), `_sequences.tsv`, `_provenance.json`, and `paml/<orthogroup>/`
#' with codeml's inputs and outputs for every model and start.
#'
#' @inheritParams phylo_busted
#' @param foreground_file Foreground for the two-ratio model; optional with alignment
#'   input (without it, only `M0` and `free` run).
#' @param models Any of "M0", "free", "two_ratio", "two_ratio_null". `M0` always runs.
#' @param codon_freq Codon frequency model: "F3x4" (default, codeml's `CodonFreq = 2`),
#'   "F1x4", or "F61".
#' @param omega_starts Initial omega values; each model is fitted from each.
#' @param kappa_start Initial transition/transversion ratio.
#' @param timeout Seconds allowed for one codeml run.
#' @param keep_runs If `TRUE` (default), `paml/<orthogroup>/` keeps codeml's control file
#'   and main output for every model and start, about 25 files per orthogroup, which lets
#'   an unchanged rerun reuse them. `FALSE` removes them once the results are in the
#'   tables, for genome-scale runs on file systems with file quotas.
#' @param codeml Path to PAML's codeml, if not on `PATH`.
#' @return Invisibly, the results table.
#' @references Yang Z (2007) PAML 4: phylogenetic analysis by maximum likelihood. Mol Biol
#'   Evol 24:1586-1591. Yang Z (1998) Likelihood ratio tests for detecting positive
#'   selection and application to primate lysozyme evolution. Mol Biol Evol 15:568-573.
#' @seealso [phylo_prepare()], [phylo_summarize()], [phylo_absrel()]
#' @export
phylo_branch_omega <- function(prepared_dir = NULL,
                               alignment_dir = NULL,
                               species_tree_file = NULL,
                               foreground_file = NULL,
                               gene_tree_dir = NULL,
                               gene_tree_mode = NULL,
                               output_dir = ".",
                               run_name = "branch_omega",
                               models = c("M0", "free", "two_ratio", "two_ratio_null"),
                               codon_freq = "F3x4",
                               omega_starts = c(0.4, 0.1, 2),
                               kappa_start = 2,
                               threads = 1L,
                               min_taxa = 4L,
                               terminal_stops = c("mask", "exclude"),
                               overwrite = FALSE,
                               timeout = 3600,
                               keep_runs = TRUE,
                               codeml = NULL) {
  .need_ape()
  terminal_stops <- match.arg(terminal_stops)
  models <- unique(c("M0", match.arg(models, c("M0", "free", "two_ratio", "two_ratio_null"), several.ok = TRUE)))
  if ("two_ratio_null" %in% models && !"two_ratio" %in% models) models <- c(models, "two_ratio")
  cf <- c(F1x4 = 1L, F3x4 = 2L, F61 = 3L)[tolower(codon_freq) == tolower(c("F1x4", "F3x4", "F61"))]
  if (length(cf) != 1L) stop("codon_freq must be F1x4, F3x4, or F61.", call. = FALSE)
  omega_starts <- suppressWarnings(as.numeric(omega_starts))
  if (!length(omega_starts) || anyNA(omega_starts) || any(omega_starts <= 0)) stop("omega_starts must be positive numbers.", call. = FALSE)
  kappa_start <- as.numeric(kappa_start); timeout <- as.integer(timeout)
  min_taxa <- as.integer(min_taxa); threads <- max(1L, as.integer(threads))
  if (is.na(min_taxa) || min_taxa < 3L) stop("min_taxa must be at least 3.", call. = FALSE)
  if (!grepl("^[A-Za-z0-9_.-]+$", run_name)) stop("run_name may contain only letters, digits, '.', '_' and '-'.", call. = FALSE)

  inp <- .phylo_resolve_inputs(prepared_dir, alignment_dir, species_tree_file, foreground_file, gene_tree_dir, gene_tree_mode,
                               require_foreground = FALSE)
  gene_tree_mode <- match.arg(inp$gene_tree_mode, c("topology", "report"))
  bin <- .codeml_path(codeml)
  version <- .codeml_version(bin)
  if (!identical(version, .paml_validated_version)) {
    message("Note: PAML ", version, " has not been checked with dndsR, which was validated with ", .paml_validated_version, ".")
  }
  tree <- .phylo_read_species_tree(inp$species_tree_file)
  hyp_fg <- if (is.null(inp$foreground_file)) NULL else .phylo_read_foreground(inp$foreground_file, tree, allow_reference = TRUE)
  if (is.null(hyp_fg)) models <- setdiff(models, c("two_ratio", "two_ratio_null"))
  files <- .phylo_list_alignments(inp$alignment_dir)
  if (!length(files) && is.null(inp$prepared)) stop("No alignments (.fa, .fas, .fasta, .fna) in ", inp$alignment_dir, call. = FALSE)
  run_dir <- file.path(output_dir, run_name)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  settings <- list(models = models, codon_freq = names(cf), codon_freq_code = unname(cf), omega_starts = omega_starts,
                   kappa_start = kappa_start, timeout = timeout, keep_runs = isTRUE(keep_runs))
  ctx <- list(method = "branch_omega", hyp = .phylo_all_branches_hypothesis(tree), hyp_fg = hyp_fg, settings = settings,
              bin = bin, version = version, run_dir = run_dir, min_taxa = min_taxa, terminal_stops = terminal_stops,
              overwrite = isTRUE(overwrite), gene_tree_dir = inp$gene_tree_dir, gene_tree_mode = gene_tree_mode,
              aliases = NULL, compare_labels = FALSE,
              prepared_tree_dir = if (is.null(inp$prepared)) NULL else inp$prepared$tree_dir)
  message(sprintf("[branch_omega] %d orthogroups, codeml %s, models %s, %d start(s) each, %d at a time",
                  length(files), version, paste(models, collapse = ","), length(omega_starts), threads))
  map_all <- function(hyp, present, gene_tree, mode, ctx) .phylo_map_foreground(hyp, present, gene_tree, mode, require_foreground = FALSE)
  res <- .parallel_lapply(names(files), function(og) .paml_one(og, files[[og]], ctx, map_all), threads = threads, what = "branch_omega")

  rows <- lapply(res, `[[`, "row")
  if (!is.null(inp$prepared)) {
    pre <- inp$prepared$summary[inp$prepared$summary$status != "prepared", , drop = FALSE]
    keep <- intersect(c("orthogroup", "status", "reason", "exclusion_basis", "detail", "n_taxa", "taxa_absent"), names(pre))
    rows <- c(rows, lapply(seq_len(nrow(pre)), function(i) c(as.list(pre[i, keep, drop = FALSE]), list(stage = "prepare"))))
  }
  results <- .phylo_bind_rows(rows)
  results <- results[order(results$orthogroup, method = "radix"), , drop = FALSE]
  results <- .paml_adjust(results)
  results$method <- "branch_omega"
  first <- c("orthogroup", "method", "status", "stage", "reason", "exclusion_basis", "detail", "n_taxa", "n_codons", "taxa_absent",
             "lnl_m0", "omega_m0", "kappa_m0", "lnl_free", "p_free_vs_m0", "p_adj_free_vs_m0",
             "two_ratio_status", "omega_bg", "omega_fg", "lnl_two_ratio", "p_two_vs_m0", "p_adj_two_vs_m0",
             "p_fg_omega_gt1", "p_adj_fg_omega_gt1")
  results <- results[, c(intersect(first, names(results)), setdiff(names(results), first)), drop = FALSE]
  results <- .calibration_column(results, "branch_omega")
  branches <- .phylo_bind_rows(lapply(res, `[[`, "branches"))
  sequences <- .phylo_bind_rows(lapply(res, `[[`, "sequences"))
  out_file <- function(what, ext = "tsv") file.path(run_dir, sprintf("%s_%s.%s", run_name, what, ext))
  write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_tsv(results, out_file("results"))
  write_tsv(branches, out_file("branches"))
  write_tsv(sequences, out_file("sequences"))

  tested <- results$status == "tested"
  label <- ifelse(tested, "tested", paste0(results$status, ":", results$reason))
  version_of <- function(pkg) if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  prov <- list(
    analysis = "branch_omega", dndsR = version_of("dndsR"), r = paste(R.version$major, R.version$minor, sep = "."),
    ape = version_of("ape"), paml_version = version, paml_validated_version = .paml_validated_version, codeml_path = bin,
    settings = settings, control_file_example = .paml_ctl(1L, omega_starts[1], 0L, kappa_start, unname(cf)),
    reproducibility = "codeml is deterministic for the same inputs, version, and starting values; the best of the starts is reported.",
    min_taxa = min_taxa, terminal_stops = terminal_stops,
    gene_tree_mode = if (is.null(inp$gene_tree_dir)) NA_character_ else gene_tree_mode,
    species_tree = list(path = inp$species_tree_file, md5 = unname(tools::md5sum(inp$species_tree_file))),
    foreground = if (is.null(inp$foreground_file)) NULL else list(path = inp$foreground_file, md5 = unname(tools::md5sum(inp$foreground_file))),
    gene_tree_dir = inp$gene_tree_dir, alignment_dir = inp$alignment_dir, prepared_dir = prepared_dir,
    multiple_testing = list(method = "BH", family = "orthogroups with status 'tested' in this run, separately for each test",
                            tests = c("p_two_vs_m0", "p_free_vs_m0", "p_fg_omega_gt1")),
    counts = as.list(table(label)),
    calibration = .calibration_notes("branch_omega", results),
    interpretation = "Branch-model omega averages over all sites of a gene. Omega above 1 on a branch is not by itself evidence of positive selection."
  )
  writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), out_file("provenance", "json"))
  message(sprintf("[branch_omega] tested %d, excluded %d, failed %d. Results: %s", sum(tested),
                  sum(results$status == "excluded"), sum(results$status == "failed"), out_file("results")))
  if (any(!tested)) {
    tab <- table(sprintf("%s:%s (%s, %s)", results$status[!tested], results$reason[!tested], results$exclusion_basis[!tested], results$stage[!tested]))
    message(paste0("  ", names(tab), ": ", as.integer(tab), collapse = "\n"))
  }
  if (any(results$status == "failed")) {
    warning(sum(results$status == "failed"), " orthogroup(s) failed. See ", out_file("results"), call. = FALSE)
  }
  if (any(tested)) .calibration_warn("branch_omega", results)
  invisible(results)
}

#' One model, every start: run (or reuse) codeml and keep the best fit
#' @keywords internal
.paml_fit <- function(ctx, dir, tree_lines, model, tips) {
  s <- ctx$settings
  code <- c(M0 = 0L, free = 1L, two_ratio = 2L, two_ratio_null = 2L)[[model]]
  fix <- as.integer(model == "two_ratio_null")
  fits <- list()
  errors <- character()
  # With the marked omega fixed at 1, codeml has no omega to start from, so one run.
  n_starts <- if (fix == 1L) 1L else length(s$omega_starts)
  for (k in seq_len(n_starts)) {
    d <- file.path(dir, model, sprintf("start%d", k))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    ctl <- .paml_ctl(code, if (fix == 1L) 1 else s$omega_starts[k], fix, s$kappa_start, s$codon_freq_code,
                     seqfile = "../../alignment.phy", treefile = "../tree.nwk")
    # The record of inputs, version, and settings goes in the control file as a codeml
    # comment, so an unchanged rerun can reuse the output without a separate file.
    record <- paste("* dndsR record:", ctx$version, unname(tools::md5sum(file.path(dir, "alignment.phy"))),
                    unname(tools::md5sum(file.path(dir, model, "tree.nwk"))))
    ctl <- c(record, ctl)
    ctl_path <- file.path(d, "codeml.ctl")
    mlc <- file.path(d, "codeml.mlc")
    reuse <- !ctx$overwrite && file.exists(ctl_path) && file.exists(mlc) && identical(readLines(ctl_path, warn = FALSE), ctl)
    if (!reuse) {
      unlink(mlc)
      writeLines(ctl, ctl_path)
      run <- .codeml_run(ctx$bin, d, s$timeout)
      if (run$timed_out) { errors <- c(errors, sprintf("start %d timed out after %d s", k, s$timeout)); next }
    }
    fit <- tryCatch(.codeml_parse(mlc, tips), error = function(e) conditionMessage(e))
    if (is.character(fit)) { unlink(mlc); errors <- c(errors, sprintf("start %d: %s", k, fit)); next }
    # codeml's other outputs are not read; its log is kept only when a run fails.
    unlink(file.path(d, .paml_aux_files))
    fit$start <- k
    fits[[length(fits) + 1L]] <- fit
  }
  if (!length(fits)) {
    return(list(ok = FALSE, reason = if (any(grepl("timed out", errors))) "codeml_timeout" else "codeml_error",
                detail = paste(model, paste(errors, collapse = "; "))))
  }
  lnl <- vapply(fits, `[[`, numeric(1), "lnl")
  best <- fits[[which.max(lnl)]]
  best$ok <- TRUE
  best$n_starts <- n_starts
  best$n_starts_ok <- length(fits)
  best$lnl_spread <- max(lnl) - min(lnl)
  best$n_starts_at_best <- sum(max(lnl) - lnl < 0.01)
  best
}

#' Branch table of one fit, matched to dndsR's branches by the taxa below each branch
#' @keywords internal
.paml_branch_values <- function(fit, branches, present, suffix) {
  side <- vapply(fit$branches$tips, .phylo_split_label, character(1), present = present)
  i <- match(branches$taxa_side, side)
  if (anyNA(i)) stop("codeml's tree does not match the branches dndsR wrote")
  b <- fit$branches[i, c("t", "N", "S", "omega", "dN", "dS", "NdN", "SdS"), drop = FALSE]
  names(b) <- paste0(names(b), "_", suffix)
  rownames(b) <- NULL
  b
}

#' Likelihood-ratio test between nested fits
#' @keywords internal
.paml_lrt <- function(alt, null) {
  if (is.null(alt) || is.null(null) || !isTRUE(alt$ok) || !isTRUE(null$ok)) return(list(lrt = NA_real_, df = NA_integer_, p = NA_real_, negative = FALSE))
  raw <- 2 * (alt$lnl - null$lnl)
  df <- alt$np - null$np
  lrt <- max(0, raw)
  list(lrt = lrt, df = df, p = stats::pchisq(lrt, df, lower.tail = FALSE), negative = raw < -1e-6)
}

#' Analyse one orthogroup with codeml
#' @keywords internal
.paml_one <- function(og, file, ctx, map_all) {
  st <- .phylo_og_setup(og, file, ctx, map_all)
  if (!is.null(st$done)) return(st$done)
  base <- st$base; info <- st$info; flags <- st$flags; map <- st$map; seq_qc <- st$seq_qc
  done <- function(status, reason, detail, extra = list(), branches = NULL) {
    list(row = .phylo_outcome(base, status, reason, detail, extra), branches = branches, sequences = seq_qc)
  }
  tp <- map$topology
  aln <- st$aln[match(tp$tree$tip.label, st$aln$taxon), , drop = FALSE]
  tips <- aln$taxon
  present <- tp$present
  branches <- cbind(orthogroup = og, map$branches[, c("branch", "terminal", "taxa_side", "species_tree_branches", "in_species_tree", "topology_source")],
                    stringsAsFactors = FALSE)
  dir <- file.path(ctx$run_dir, "paml", og)
  models <- ctx$settings$models

  # Two-ratio: the foreground must be identifiable on this orthogroup's tree.
  two <- list(status = "not_run", reason = NA_character_, marks = NULL)
  if (any(c("two_ratio", "two_ratio_null") %in% models)) {
    fgmap <- .phylo_map_foreground(ctx$hyp_fg, st$present, st$gene_tree, ctx$gene_tree_mode)
    if (!fgmap$ok) {
      two <- list(status = "excluded", reason = fgmap$reason, marks = NULL)
      models <- setdiff(models, c("two_ratio", "two_ratio_null"))
    } else {
      fg_sides <- fgmap$branches$taxa_side[fgmap$branches$foreground]
      branches$foreground <- branches$taxa_side %in% fg_sides
      marks <- logical(ape::Ntip(tp$tree) + tp$tree$Nnode)
      marks[tp$child[branches$foreground[match(tp$node_names[tp$child], branches$branch)]]] <- TRUE
      two <- list(status = "tested", reason = NA_character_, marks = marks)
    }
  }
  if (is.null(branches$foreground)) branches$foreground <- NA
  fits <- list()
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  .paml_write_phylip(aln, file.path(dir, "alignment.phy"))
  for (model in models) {
    tree_lines <- .paml_tree_lines(tp$tree, if (model %in% c("two_ratio", "two_ratio_null")) two$marks else NULL)
    dir.create(file.path(dir, model), showWarnings = FALSE)
    writeLines(tree_lines, file.path(dir, model, "tree.nwk"))
    fits[[model]] <- .paml_fit(ctx, dir, tree_lines, model, tips)
    if (!isTRUE(fits[[model]]$ok)) {
      if (model == "M0") return(done("failed", fits[[model]]$reason, fits[[model]]$detail, info, branches))
      if (model %in% c("two_ratio", "two_ratio_null")) two$status <- "failed"
      flags <- c(flags, sprintf("%s_failed", model))
    }
  }
  ok <- function(m) !is.null(fits[[m]]) && isTRUE(fits[[m]]$ok)
  for (m in names(fits)) if (ok(m)) {
    sfx <- c(M0 = "m0", free = "free", two_ratio = "two", two_ratio_null = "two_null")[[m]]
    vals <- tryCatch(.paml_branch_values(fits[[m]], branches, present, sfx), error = function(e) NULL)
    if (is.null(vals)) return(done("failed", "codeml_output_invalid", paste(m, ": codeml's tree does not match the branches dndsR wrote"), info, branches))
    if (m != "two_ratio_null") branches <- cbind(branches, vals, stringsAsFactors = FALSE)
    if (fits[[m]]$n_starts_at_best < fits[[m]]$n_starts_ok) flags <- c(flags, sprintf("starts_disagree:%s", m))
    if (length(fits[[m]]$warnings)) flags <- c(flags, sprintf("codeml_warning:%s", m))
  }
  # Per-branch reliability, from the free-ratio fit when present, otherwise M0.
  src <- if (ok("free")) "free" else "m0"
  get <- function(v) branches[[paste0(v, "_", src)]]
  bflags <- character(nrow(branches))
  add <- function(cond, flag) bflags <<- ifelse(cond & !is.na(cond), ifelse(nzchar(bflags), paste(bflags, flag, sep = ";"), flag), bflags)
  add(get("t") < 1e-5, "zero_length")
  add(get("SdS") < 1, "few_synonymous_changes")
  if (src == "free") add(get("omega") >= 998, "omega_at_upper_bound")
  add(get("dS") > 1, "high_ds")
  add(grepl("+", branches$species_tree_branches, fixed = TRUE), "merged_branch")
  add(!branches$in_species_tree, "not_in_species_tree")
  branches$flags <- bflags
  branches$reliable_free <- if (src == "free") !grepl("zero_length|few_synonymous_changes|omega_at_upper_bound|high_ds", bflags) else NA

  t2 <- .paml_lrt(fits$two_ratio, fits$M0)
  tf <- .paml_lrt(fits$free, fits$M0)
  tn <- .paml_lrt(fits$two_ratio, fits$two_ratio_null)
  for (x in list(list(t2, "two_vs_m0"), list(tf, "free_vs_m0"), list(tn, "fg_omega_eq1"))) if (isTRUE(x[[1]]$negative)) flags <- c(flags, paste0("lrt_negative:", x[[2]]))
  omega_fg <- if (ok("two_ratio")) fits$two_ratio$omega[2] else NA_real_
  p_gt1 <- if (is.na(tn$lrt) || is.na(omega_fg)) NA_real_ else stats::pnorm(sign(omega_fg - 1) * sqrt(tn$lrt), lower.tail = FALSE)
  fit_cols <- function(m, sfx) {
    f <- fits[[m]]
    if (!ok(m)) return(list())
    stats::setNames(list(f$lnl, f$np, f$kappa, f$tree_length, f$lnl_spread, f$n_starts_at_best, f$n_starts_ok),
                    paste0(c("lnl_", "np_", "kappa_", "tree_length_", "lnl_spread_", "starts_at_best_", "starts_ok_"), sfx))
  }
  values <- c(
    fit_cols("M0", "m0"), list(omega_m0 = fits$M0$omega[1]),
    fit_cols("free", "free"), list(lrt_free_vs_m0 = tf$lrt, df_free_vs_m0 = tf$df, p_free_vs_m0 = tf$p),
    list(two_ratio_status = two$status, two_ratio_reason = two$reason, n_fg_branches = sum(branches$foreground %in% TRUE)),
    fit_cols("two_ratio", "two_ratio"),
    list(omega_bg = if (ok("two_ratio")) fits$two_ratio$omega[1] else NA_real_, omega_fg = omega_fg,
         lrt_two_vs_m0 = t2$lrt, p_two_vs_m0 = t2$p),
    fit_cols("two_ratio_null", "two_ratio_null"),
    list(lrt_fg_omega_eq1 = tn$lrt, p_fg_omega_gt1 = p_gt1,
         n_branches = nrow(branches), n_branches_reliable = sum(branches$reliable_free %in% TRUE))
  )
  info$qc_flags <- paste(unique(flags), collapse = ";")
  if (!isTRUE(ctx$settings$keep_runs)) unlink(dir, recursive = TRUE)
  list(row = .phylo_outcome(base, "tested", NA_character_, NA_character_, c(info, values)), branches = branches, sequences = seq_qc)
}

#' Benjamini-Hochberg for each branch-model test, across tested orthogroups
#' @keywords internal
.paml_adjust <- function(results) {
  tested <- results$status == "tested"
  for (p in c("p_two_vs_m0", "p_free_vs_m0", "p_fg_omega_gt1")) {
    adj <- rep(NA_real_, nrow(results))
    if (!is.null(results[[p]])) {
      i <- which(tested & !is.na(results[[p]]))
      if (length(i)) adj[i] <- stats::p.adjust(results[[p]][i], method = "BH")
    }
    results[[sub("^p_", "p_adj_", p)]] <- adj
  }
  results
}
