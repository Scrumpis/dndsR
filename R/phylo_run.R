# Method-agnostic runner for the phylogenetic module: inputs, QC, foreground mapping,
# HyPhy execution, multiple starts, multiple testing, and provenance. Method specifics come
# from `.hyphy_methods`.
#
# Every orthogroup ends with one status: `tested` (the method ran and its output was read),
# `excluded` (a rule kept it out; `exclusion_basis` says which kind), or `failed` (a step
# errored). `stage` says whether this happened in phylo_prepare or in the analysis.

#' Bind rows with differing columns
#' @keywords internal
.phylo_bind_rows <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    r <- as.list(r)
    for (cl in setdiff(cols, names(r))) r[[cl]] <- NA
    as.data.frame(r[cols], stringsAsFactors = FALSE, check.names = FALSE)
  }))
}

#' List orthogroup alignments in a folder
#' @keywords internal
.phylo_list_alignments <- function(alignment_dir) {
  if (is.null(alignment_dir) || !dir.exists(alignment_dir)) {
    stop("Alignment folder not found: ", alignment_dir, call. = FALSE)
  }
  files <- sort(list.files(alignment_dir, pattern = "\\.(fa|fas|fasta|fna)$", ignore.case = TRUE, full.names = TRUE), method = "radix")
  og <- sub("\\.[^.]+$", "", basename(files))
  bad <- og[!grepl("^[A-Za-z0-9_.-]+$", og)]
  if (length(bad)) stop("Orthogroup file names may contain only letters, digits, '.', '_' and '-': ", paste(utils::head(bad, 5), collapse = ", "), call. = FALSE)
  if (anyDuplicated(og)) stop("Two alignments share an orthogroup name: ", paste(unique(og[duplicated(og)]), collapse = ", "), call. = FALSE)
  stats::setNames(files, og)
}

#' @keywords internal
.phylo_find_gene_tree <- function(gene_tree_dir, og) {
  if (is.null(gene_tree_dir)) return(NULL)
  hits <- file.path(gene_tree_dir, paste0(og, c(".nwk", ".newick", ".tre", ".tree", ".treefile", "_tree.txt")))
  hits <- hits[file.exists(hits)]
  if (length(hits)) hits[1] else NA_character_
}

#' One outcome row
#' @keywords internal
.phylo_outcome <- function(base, status, reason, detail, extra = list()) {
  c(base, list(status = status, reason = reason, exclusion_basis = .phylo_reason_basis(reason), detail = detail), extra)
}

#' Run one HyPhy start and read it
#' @keywords internal
.phylo_start <- function(ctx, m, dir, seed, aln, tree) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  json_path <- file.path(dir, paste0(ctx$method, ".json"))
  run <- .hyphy_run(ctx$bin, m$hyphy_name, m$args(ctx$settings, aln, tree, basename(json_path)), dir, ctx$hyphy_cpu, seed)
  if (run$status != 0L) {
    return(list(ok = FALSE, reason = "hyphy_error",
                detail = if (is.na(run$error)) paste("HyPhy exited with status", run$status) else run$error))
  }
  .phylo_read_start(ctx, m, dir)
}

#' @keywords internal
.phylo_read_start <- function(ctx, m, dir) {
  json_path <- file.path(dir, paste0(ctx$method, ".json"))
  parsed <- tryCatch({
    j <- .hyphy_read_json(json_path, m$required)
    list(values = m$parse(j), tested = .hyphy_tested_branches(j))
  }, error = function(e) conditionMessage(e))
  if (is.character(parsed)) return(list(ok = FALSE, reason = "hyphy_output_invalid", detail = parsed))
  log <- file.path(dir, "hyphy.log")
  log_text <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
  # With kill-zero-lengths (HyPhy's default), internal branches estimated at zero length
  # under the nucleotide model are deleted before the codon models are fitted.
  del <- grep("Deleted [0-9]+ zero-length internal branches:", log_text, value = TRUE)
  deleted <- if (length(del)) trimws(strsplit(gsub(".*: *`|`.*", "", del[1]), ",")[[1]]) else character()
  c(list(ok = TRUE, json = json_path, refit = any(grepl("Negative test LRT", log_text, fixed = TRUE)), deleted = deleted), parsed)
}

#' Diagnostics across independent starts. They never replace the primary result.
#'
#' The best alternative is the largest unconstrained log-likelihood found. The best null is
#' the largest constrained one; where HyPhy did not fit the constrained model because the
#' unconstrained fit showed no evidence, that fit is treated as a point of the null.
#' @keywords internal
.phylo_start_diagnostics <- function(starts, seeds) {
  unc <- vapply(starts, function(s) s$values$lnl_unconstrained, numeric(1))
  con <- vapply(starts, function(s) s$values$lnl_constrained, numeric(1))
  null_points <- ifelse(is.na(con), unc, con)
  lrt <- vapply(starts, function(s) s$values$lrt, numeric(1))
  best_alt <- max(unc)
  best_null <- max(null_points)
  lrt_best <- max(0, 2 * (best_alt - best_null))
  list(n_starts = length(starts), start_seeds = paste(seeds, collapse = ","),
       lrt_min_starts = min(lrt), lrt_max_starts = max(lrt),
       lnl_unconstrained_best = best_alt, lnl_null_best = best_null,
       lrt_best_fits = lrt_best, p_best_fits = 0.5 * stats::pchisq(lrt_best, 2, lower.tail = FALSE),
       primary_unconstrained_is_best = abs(unc[1] - best_alt) < 1e-6,
       primary_null_is_best = abs(null_points[1] - best_null) < 1e-6)
}

#' Benjamini-Hochberg correction across the orthogroups tested in a run
#'
#' `p_adj` corrects `p_value`, from the primary start, as always. With several starts,
#' `p_best_fits_adj` corrects `p_best_fits`, from the best alternative and best null fits
#' across starts (the stage 2 recommendation), over the same family: every orthogroup with
#' status `tested`. Any failed start makes an orthogroup `failed`, so every tested one has
#' `p_best_fits` when starts were added. Excluded and failed orthogroups get `NA`.
#'
#' @param results One row per orthogroup, with `status`, `p_value`, and, after several
#'   starts, `p_best_fits`.
#' @return `results` with `p_adj`, `p_adj_method`, `family_size`, and, when `p_best_fits`
#'   is present, `p_best_fits_adj`.
#' @keywords internal
.phylo_adjust <- function(results) {
  tested <- results$status == "tested"
  # With nothing tested (every orthogroup excluded or failed), there is no p_value column.
  bh <- function(p) {
    out <- rep(NA_real_, nrow(results))
    if (!is.null(p) && any(tested)) out[tested] <- stats::p.adjust(p[tested], method = "BH")
    out
  }
  results$p_adj <- bh(results$p_value)
  results$p_adj_method <- "BH"
  results$family_size <- sum(tested)
  if ("p_best_fits" %in% names(results)) results$p_best_fits_adj <- bh(results$p_best_fits)
  results
}

#' Analyse one orthogroup
#' @keywords internal
.phylo_one <- function(og, file, ctx) {
  base <- list(orthogroup = og, stage = "analysis", alignment_md5 = unname(tools::md5sum(file)))
  done <- function(status, reason, detail, extra = list(), branches = NULL, sequences = NULL) {
    list(row = .phylo_outcome(base, status, reason, detail, extra), branches = branches, sequences = sequences)
  }
  aln <- .phylo_read_alignment(file)
  qc <- .phylo_check_alignment(aln, ctx$hyp$taxa, ctx$min_taxa, ctx$terminal_stops)
  if (!qc$ok) return(done("excluded", qc$reason, qc$detail))
  aln <- qc$aln
  present <- aln$taxon
  info <- c(as.list(qc$metrics), list(taxa_absent = paste(setdiff(ctx$hyp$taxa, present), collapse = ",")))
  seq_qc <- cbind(orthogroup = og, qc$sequences, stringsAsFactors = FALSE)

  gene_tree <- NULL
  gt_path <- .phylo_find_gene_tree(ctx$gene_tree_dir, og)
  if (!is.null(gt_path)) {
    if (is.na(gt_path)) return(done("excluded", "gene_tree_missing", paste("no gene tree for", og, "in", ctx$gene_tree_dir), info, sequences = seq_qc))
    g <- .phylo_try_gene_tree(gt_path, aln, ctx$aliases)
    if (is.null(g$tree)) return(done(g$status, g$reason, g$detail, info, sequences = seq_qc))
    gene_tree <- g$tree
  }
  map <- .phylo_map_foreground(ctx$hyp, present, gene_tree, ctx$gene_tree_mode)
  info <- c(info, map$info)
  if (!map$ok) return(done("excluded", map$reason, map$detail, info, sequences = seq_qc))
  fg <- sort(map$branches$branch[map$branches$foreground], method = "radix")
  if (!is.null(ctx$prepared_tree_dir)) {
    expected <- file.path(ctx$prepared_tree_dir, paste0(og, ".nwk"))
    if (!file.exists(expected) || !identical(readLines(expected, warn = FALSE), map$newick)) {
      return(done("failed", "prepared_tree_mismatch", paste("the labelled tree differs from the one phylo_prepare wrote:", expected), info, sequences = seq_qc))
    }
  }
  flags <- qc$flags
  if (isTRUE(map$info$fg_absent_branches > 0)) flags <- c(flags, sprintf("foreground_branches_absent:%d", map$info$fg_absent_branches))
  if (isTRUE(map$info$gt_missing_species_splits > 0)) flags <- c(flags, "gene_tree_discordant")
  if (isTRUE(map$info$gt_discordant_near_foreground)) flags <- c(flags, "gene_tree_discordant_near_foreground")
  info <- c(info, list(n_fg_branches = length(fg), fg_branches = paste(fg, collapse = ",")))

  dir <- file.path(ctx$run_dir, "hyphy", og)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  aln_path <- file.path(dir, "alignment.fas")
  tree_path <- file.path(dir, "tree.nwk")
  record_path <- file.path(dir, "record.json")
  .phylo_write_alignment(aln, aln_path)
  writeLines(map$newick, tree_path)
  seeds <- c(ctx$seed, ctx$extra_seeds)
  record <- list(method = ctx$method, settings = ctx$settings, seed = ctx$seed, extra_seeds = ctx$extra_seeds,
                 hyphy_cpu = ctx$hyphy_cpu, hyphy_version = ctx$hyphy_version,
                 alignment_md5 = unname(tools::md5sum(aln_path)), tree_md5 = unname(tools::md5sum(tree_path)))
  record_json <- as.character(jsonlite::toJSON(record, auto_unbox = TRUE, pretty = TRUE, digits = NA))
  reuse <- !ctx$overwrite && file.exists(record_path) &&
    identical(paste(readLines(record_path, warn = FALSE), collapse = "\n"), record_json)
  branches <- cbind(orthogroup = og, map$branches, stringsAsFactors = FALSE)
  m <- .hyphy_methods[[ctx$method]]
  start_dirs <- c(dir, file.path(dir, paste0("start_seed", ctx$extra_seeds)))
  if (!reuse) {
    unlink(c(record_path, file.path(dir, paste0(ctx$method, ".json"))))
    unlink(list.files(dir, pattern = "^start_seed", full.names = TRUE), recursive = TRUE)
  }
  starts <- list()
  for (k in seq_along(seeds)) {
    rel <- if (k == 1L) "" else "../"
    s <- if (reuse) .phylo_read_start(ctx, m, start_dirs[k]) else
      .phylo_start(ctx, m, start_dirs[k], seeds[k], paste0(rel, "alignment.fas"), paste0(rel, "tree.nwk"))
    if (!isTRUE(s$ok)) {
      return(done("failed", s$reason, paste0(if (k > 1L) sprintf("start with seed %d: ", seeds[k]), s$detail),
                  info, branches, seq_qc))
    }
    if (!identical(s$tested, fg)) {
      gone <- intersect(s$deleted, fg)
      if (length(gone) && all(s$tested %in% fg) && setequal(c(s$tested, gone), fg)) {
        return(done("excluded", "foreground_branch_zero_length", paste0(
          "HyPhy deleted requested branch(es) ", paste(gone, collapse = ","), " as zero-length (kill_zero_lengths = Yes), ",
          "so it would test only ", paste(s$tested, collapse = ","), ", a narrower hypothesis"), info, branches, seq_qc))
      }
      return(done("failed", "label_mismatch", paste0("HyPhy tested ", paste(s$tested, collapse = ","), "; dndsR labelled ", paste(fg, collapse = ",")),
                  info, branches, seq_qc))
    }
    starts[[k]] <- s
  }
  if (!reuse) writeLines(record_json, record_path)
  if (starts[[1]]$refit) flags <- c(flags, "hyphy_refit_after_negative_lrt")
  if (length(starts[[1]]$deleted)) flags <- c(flags, sprintf("hyphy_deleted_zero_length_branches:%d", length(starts[[1]]$deleted)))
  info$qc_flags <- paste(flags, collapse = ";")
  diag <- if (length(seeds) > 1L) .phylo_start_diagnostics(starts, seeds) else list(n_starts = 1L)
  row <- .phylo_outcome(base, "tested", NA_character_, NA_character_,
                        c(info, starts[[1]]$values, diag,
                          list(json = file.path("hyphy", og, basename(starts[[1]]$json)),
                               json_md5 = unname(tools::md5sum(starts[[1]]$json)))))
  list(row = row, branches = branches, sequences = seq_qc)
}

#' Run a HyPhy method over a folder of orthogroup alignments
#' @keywords internal
.phylo_run <- function(method, alignment_dir, species_tree_file, foreground_file, gene_tree_dir,
                       output_dir, run_name, settings, seed, hyphy_cpu, threads, min_taxa,
                       terminal_stops, overwrite, hyphy, extra_seeds = NULL,
                       gene_tree_mode = NULL, prepared_dir = NULL) {
  .need_ape()
  prepared <- NULL
  if (!is.null(prepared_dir)) {
    if (!is.null(alignment_dir) || !is.null(species_tree_file) || !is.null(foreground_file)) {
      stop("With prepared_dir, the alignments, species tree, and foreground come from phylo_prepare; do not give them too.", call. = FALSE)
    }
    prepared <- .phylo_read_manifest(prepared_dir)
    alignment_dir <- prepared$alignment_dir
    species_tree_file <- prepared$species_tree_file
    foreground_file <- prepared$foreground_file
    # The hypothesis must be the one prepared: same species tree and foreground files.
    for (what in c("species_tree", "foreground")) {
      file <- if (what == "species_tree") species_tree_file else foreground_file
      if (!identical(unname(tools::md5sum(file)), prepared$inputs[[what]]$md5)) {
        stop("The ", what, " in ", prepared_dir, " differs from the one phylo_prepare used. Rerun phylo_prepare.", call. = FALSE)
      }
    }
    if (is.null(gene_tree_dir)) gene_tree_dir <- prepared$gene_tree_dir
    if (is.null(gene_tree_mode)) gene_tree_mode <- prepared$gene_tree_mode
  }
  if (is.null(gene_tree_mode)) gene_tree_mode <- "topology" else if (is.null(alignment_dir) || is.null(species_tree_file) || is.null(foreground_file)) {
    stop("Give prepared_dir (from phylo_prepare), or alignment_dir, species_tree_file, and foreground_file.", call. = FALSE)
  }
  if (!is.null(gene_tree_dir) && !dir.exists(gene_tree_dir)) stop("Gene tree folder not found: ", gene_tree_dir, call. = FALSE)
  gene_tree_mode <- match.arg(gene_tree_mode, c("topology", "report"))
  terminal_stops <- match.arg(terminal_stops, c("mask", "exclude"))
  seed <- as.integer(seed); hyphy_cpu <- as.integer(hyphy_cpu); threads <- max(1L, as.integer(threads))
  min_taxa <- as.integer(min_taxa)
  extra_seeds <- if (length(extra_seeds)) as.integer(extra_seeds) else integer()
  if (is.na(seed)) stop("seed must be an integer.", call. = FALSE)
  if (anyNA(extra_seeds) || seed %in% extra_seeds || anyDuplicated(extra_seeds)) {
    stop("extra_seeds must be distinct integers, different from seed.", call. = FALSE)
  }
  if (is.na(hyphy_cpu) || hyphy_cpu < 1L) stop("hyphy_cpu must be a positive integer.", call. = FALSE)
  if (is.na(min_taxa) || min_taxa < 3L) stop("min_taxa must be at least 3.", call. = FALSE)
  if (!grepl("^[A-Za-z0-9_.-]+$", run_name)) stop("run_name may contain only letters, digits, '.', '_' and '-'.", call. = FALSE)
  bin <- .hyphy_path(hyphy)
  version <- .hyphy_version(bin)
  if (!identical(version, .hyphy_validated_version)) {
    message("Note: HyPhy ", version, " has not been checked with dndsR, which was validated with ", .hyphy_validated_version, ".")
  }
  tree <- .phylo_read_species_tree(species_tree_file)
  hyp <- .phylo_read_foreground(foreground_file, tree)
  files <- .phylo_list_alignments(alignment_dir)
  if (!length(files) && is.null(prepared)) stop("No alignments (.fa, .fas, .fasta, .fna) in ", alignment_dir, call. = FALSE)
  run_dir <- file.path(output_dir, run_name)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  ctx <- list(method = method, hyp = hyp, settings = settings, seed = seed, extra_seeds = extra_seeds,
              hyphy_cpu = hyphy_cpu, hyphy_version = version, bin = bin, run_dir = run_dir, min_taxa = min_taxa,
              terminal_stops = terminal_stops, overwrite = isTRUE(overwrite), gene_tree_dir = gene_tree_dir,
              gene_tree_mode = gene_tree_mode, aliases = NULL,
              prepared_tree_dir = if (is.null(prepared)) NULL else prepared$tree_dir)
  message(sprintf("[%s] %d orthogroups, HyPhy %s, %d at a time with %d HyPhy thread(s) each%s", method, length(files), version,
                  threads, hyphy_cpu, if (length(extra_seeds)) sprintf(", %d starts each", 1L + length(extra_seeds)) else ""))
  res <- .parallel_lapply(names(files), function(og) .phylo_one(og, files[[og]], ctx), threads = threads, what = method)

  rows <- lapply(res, `[[`, "row")
  if (!is.null(prepared)) {
    pre <- prepared$summary[prepared$summary$status != "prepared", , drop = FALSE]
    keep <- intersect(c("orthogroup", "status", "reason", "exclusion_basis", "detail", "n_taxa", "taxa_absent"), names(pre))
    rows <- c(rows, lapply(seq_len(nrow(pre)), function(i) c(as.list(pre[i, keep, drop = FALSE]), list(stage = "prepare"))))
  }
  results <- .phylo_bind_rows(rows)
  results <- results[order(results$orthogroup, method = "radix"), , drop = FALSE]
  results <- .phylo_adjust(results)
  tested <- results$status == "tested"
  results$method <- method
  first <- c("orthogroup", "method", "status", "stage", "reason", "exclusion_basis", "detail", "n_taxa", "n_codons",
             "taxa_absent", "n_fg_branches", "fg_branches", "lrt", "p_value", "p_adj", "p_adj_method", "family_size")
  results <- results[, c(intersect(first, names(results)), setdiff(names(results), first)), drop = FALSE]
  if ("p_best_fits_adj" %in% names(results)) {   # beside p_best_fits, which it corrects
    rest <- setdiff(names(results), "p_best_fits_adj")
    at <- match("p_best_fits", rest)
    results <- results[, c(rest[seq_len(at)], "p_best_fits_adj", rest[-seq_len(at)]), drop = FALSE]
  }
  branches <- .phylo_bind_rows(lapply(res, `[[`, "branches"))
  sequences <- .phylo_bind_rows(lapply(res, `[[`, "sequences"))

  out_file <- function(what, ext = "tsv") file.path(run_dir, sprintf("%s_%s.%s", run_name, what, ext))
  write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_tsv(results, out_file("results"))
  write_tsv(branches, out_file("branches"))
  write_tsv(sequences, out_file("sequences"))

  label <- ifelse(results$status == "tested", "tested", paste0(results$status, ":", results$reason))
  absent <- strsplit(ifelse(is.na(results$taxa_absent), "", results$taxa_absent), ",")
  presence <- vapply(hyp$taxa, function(tx) sum(tested & !vapply(absent, function(a) tx %in% a, logical(1))), integer(1))
  version_of <- function(pkg) if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  prov <- list(
    analysis = method, dndsR = version_of("dndsR"), r = paste(R.version$major, R.version$minor, sep = "."),
    ape = version_of("ape"), hyphy_version = version, hyphy_validated_version = .hyphy_validated_version,
    hyphy_path = bin, settings = settings, random_seed = seed, extra_seeds = extra_seeds, hyphy_cpu = hyphy_cpu,
    reproducibility = paste("Results are reproducible for the same HyPhy build, settings, random_seed, and hyphy_cpu.",
                            "threads only sets how many orthogroups run at once. p_value and p_adj always come",
                            "from random_seed; with extra seeds, p_best_fits and p_best_fits_adj combine every start."),
    min_taxa = min_taxa, terminal_stops = terminal_stops,
    gene_tree_mode = if (is.null(gene_tree_dir)) NA_character_ else gene_tree_mode,
    species_tree = list(path = species_tree_file, md5 = unname(tools::md5sum(species_tree_file))),
    foreground = list(path = foreground_file, md5 = unname(tools::md5sum(foreground_file)),
                      hypothesis = hyp$spec, foreground_taxa = hyp$foreground_taxa),
    gene_tree_dir = gene_tree_dir, alignment_dir = alignment_dir, prepared_dir = prepared_dir,
    prepared_provenance_md5 = if (is.null(prepared)) NULL else prepared$provenance_md5,
    multiple_testing = list(method = "BH", family = "orthogroups with status 'tested' in this run, under one foreground hypothesis and one set of settings",
                            family_size = sum(tested)),
    counts = as.list(table(label)), exclusions_by_basis = as.list(table(results$exclusion_basis[!tested])),
    tested_orthogroups_per_taxon = as.list(presence),
    ascertainment = "Only single-copy orthogroups that pass QC and keep the foreground identifiable are tested. Excluded and failed orthogroups are listed with reasons; conclusions apply to the tested set only."
  )
  writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), out_file("provenance", "json"))

  message(sprintf("[%s] tested %d, excluded %d, failed %d. Results: %s", method, sum(tested),
                  sum(results$status == "excluded"), sum(results$status == "failed"), out_file("results")))
  if (any(!tested)) {
    tab <- table(sprintf("%s:%s (%s, %s)", results$status[!tested], results$reason[!tested], results$exclusion_basis[!tested], results$stage[!tested]))
    message(paste0("  ", names(tab), ": ", as.integer(tab), collapse = "\n"))
  }
  if (any(results$status == "failed")) {
    warning(sum(results$status == "failed"), " orthogroup(s) failed and are not in the multiple-testing family. See ",
            out_file("results"), call. = FALSE)
  }
  invisible(results)
}
