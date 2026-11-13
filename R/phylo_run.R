# Method-agnostic runner for the phylogenetic module: inputs, QC, foreground mapping,
# HyPhy execution, multiple testing, and provenance. Method specifics come from
# `.hyphy_methods`.

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
  if (!length(files)) stop("No alignments (.fa, .fas, .fasta, .fna) in ", alignment_dir, call. = FALSE)
  og <- sub("\\.[^.]+$", "", basename(files))
  bad <- og[!grepl("^[A-Za-z0-9_.-]+$", og)]
  if (length(bad)) stop("Orthogroup file names may contain only letters, digits, '.', '_' and '-': ", paste(utils::head(bad, 5), collapse = ", "), call. = FALSE)
  if (anyDuplicated(og)) stop("Two alignments share an orthogroup name: ", paste(unique(og[duplicated(og)]), collapse = ", "), call. = FALSE)
  stats::setNames(files, og)
}

#' @keywords internal
.phylo_find_gene_tree <- function(gene_tree_dir, og) {
  if (is.null(gene_tree_dir)) return(NULL)
  hits <- file.path(gene_tree_dir, paste0(og, c(".nwk", ".newick", ".tre", ".tree", ".treefile")))
  hits <- hits[file.exists(hits)]
  if (length(hits)) hits[1] else NA_character_
}

#' Analyse one orthogroup
#' @keywords internal
.phylo_one <- function(og, file, ctx) {
  base <- list(orthogroup = og, status = "skipped", reason = NA_character_, detail = NA_character_,
               alignment_md5 = unname(tools::md5sum(file)))
  skip <- function(reason, detail, extra = list()) {
    list(row = c(base[setdiff(names(base), c("reason", "detail"))], list(reason = reason, detail = detail), extra),
         branches = NULL, sequences = NULL)
  }
  aln <- .phylo_read_alignment(file)
  qc <- .phylo_check_alignment(aln, ctx$hyp$taxa, ctx$min_taxa, ctx$terminal_stops)
  if (!qc$ok) return(skip(qc$reason, qc$detail))
  aln <- qc$aln
  present <- aln$taxon
  info <- c(as.list(qc$metrics), list(taxa_absent = paste(setdiff(ctx$hyp$taxa, present), collapse = ",")))
  seq_qc <- cbind(orthogroup = og, qc$sequences, stringsAsFactors = FALSE)

  gene_tree <- NULL
  gt_path <- .phylo_find_gene_tree(ctx$gene_tree_dir, og)
  if (!is.null(gt_path)) {
    if (is.na(gt_path)) return(skip("gene_tree_missing", paste("no gene tree for", og, "in", ctx$gene_tree_dir), info))
    gene_tree <- tryCatch(.phylo_read_gene_tree(gt_path, stats::setNames(aln$taxon, aln$seq_name)),
                          error = function(e) conditionMessage(e))
    if (is.character(gene_tree)) return(skip(gene_tree, gt_path, info))
  }
  map <- .phylo_map_foreground(ctx$hyp, present, gene_tree)
  if (!map$ok) {
    r <- skip(map$reason, map$detail, info)
    r$sequences <- seq_qc
    return(r)
  }
  fg <- sort(map$branches$branch[map$branches$foreground], method = "radix")
  info <- c(info, list(n_fg_branches = length(fg), fg_branches = paste(fg, collapse = ","),
                       qc_flags = paste(qc$flags, collapse = ";")))

  dir <- file.path(ctx$run_dir, "hyphy", og)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  aln_path <- file.path(dir, "alignment.fas")
  tree_path <- file.path(dir, "tree.nwk")
  json_path <- file.path(dir, paste0(ctx$method, ".json"))
  record_path <- file.path(dir, "record.json")
  .phylo_write_alignment(aln, aln_path)
  writeLines(map$newick, tree_path)
  record <- list(method = ctx$method, settings = ctx$settings, seed = ctx$seed, hyphy_cpu = ctx$hyphy_cpu,
                 hyphy_version = ctx$hyphy_version,
                 alignment_md5 = unname(tools::md5sum(aln_path)), tree_md5 = unname(tools::md5sum(tree_path)))
  record_json <- jsonlite::toJSON(record, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  reuse <- !ctx$overwrite && file.exists(json_path) && file.exists(record_path) &&
    identical(paste(readLines(record_path, warn = FALSE), collapse = "\n"), as.character(record_json))
  branches <- cbind(orthogroup = og, map$branches, stringsAsFactors = FALSE)
  m <- .hyphy_methods[[ctx$method]]
  if (!reuse) {
    unlink(c(json_path, record_path))
    run <- .hyphy_run(ctx$bin, m$hyphy_name, m$args(ctx$settings, "alignment.fas", "tree.nwk", basename(json_path)),
                      dir, ctx$hyphy_cpu, ctx$seed)
    if (run$status != 0L) {
      r <- skip("hyphy_error", if (is.na(run$error)) paste("HyPhy exited with status", run$status) else run$error, info)
      r$row$status <- "failed"
      r$branches <- branches
      r$sequences <- seq_qc
      return(r)
    }
  }
  parsed <- tryCatch({
    j <- .hyphy_read_json(json_path, m$required)
    list(values = m$parse(j), tested = .hyphy_tested_branches(j))
  }, error = function(e) conditionMessage(e))
  fail <- function(reason, detail) {
    r <- skip(reason, detail, info)
    r$row$status <- "failed"
    r$branches <- branches
    r$sequences <- seq_qc
    r
  }
  if (is.character(parsed)) return(fail("hyphy_output_invalid", parsed))
  if (!identical(parsed$tested, fg)) {
    return(fail("label_mismatch", paste0("HyPhy tested ", paste(parsed$tested, collapse = ","), "; dndsR labelled ", paste(fg, collapse = ","))))
  }
  if (!reuse) writeLines(record_json, record_path)
  log_text <- readLines(file.path(dir, "hyphy.log"), warn = FALSE)
  notes <- character()
  if (any(grepl("Negative test LRT", log_text, fixed = TRUE))) notes <- c(notes, "hyphy_refit_after_negative_lrt")
  if (length(notes)) info$qc_flags <- paste(c(if (nzchar(info$qc_flags)) info$qc_flags, notes), collapse = ";")
  row <- c(base, list(status = "tested"), info, parsed$values,
           list(json = file.path("hyphy", og, basename(json_path)), json_md5 = unname(tools::md5sum(json_path))))
  row <- row[!duplicated(names(row), fromLast = TRUE)]
  list(row = row, branches = branches, sequences = seq_qc)
}

#' Run a HyPhy method over a folder of orthogroup alignments
#' @keywords internal
.phylo_run <- function(method, alignment_dir, species_tree_file, foreground_file, gene_tree_dir,
                       output_dir, run_name, settings, seed, hyphy_cpu, threads, min_taxa,
                       terminal_stops, overwrite, hyphy) {
  .need_ape()
  if (!is.null(gene_tree_dir) && !dir.exists(gene_tree_dir)) stop("Gene tree folder not found: ", gene_tree_dir, call. = FALSE)
  terminal_stops <- match.arg(terminal_stops, c("mask", "exclude"))
  seed <- as.integer(seed); hyphy_cpu <- as.integer(hyphy_cpu); threads <- max(1L, as.integer(threads))
  min_taxa <- as.integer(min_taxa)
  if (is.na(seed)) stop("seed must be an integer.", call. = FALSE)
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
  run_dir <- file.path(output_dir, run_name)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  ctx <- list(method = method, hyp = hyp, settings = settings, seed = seed, hyphy_cpu = hyphy_cpu,
              hyphy_version = version, bin = bin, run_dir = run_dir, min_taxa = min_taxa,
              terminal_stops = terminal_stops, overwrite = isTRUE(overwrite), gene_tree_dir = gene_tree_dir)
  message(sprintf("[%s] %d orthogroups, HyPhy %s, %d at a time with %d HyPhy thread(s) each", method, length(files), version, threads, hyphy_cpu))
  res <- .parallel_lapply(names(files), function(og) .phylo_one(og, files[[og]], ctx), threads = threads, what = method)

  results <- .phylo_bind_rows(lapply(res, `[[`, "row"))
  tested <- results$status == "tested"
  results$p_adj <- NA_real_
  results$p_adj[tested] <- stats::p.adjust(results$p_value[tested], method = "BH")
  results$p_adj_method <- "BH"
  results$family_size <- sum(tested)
  results$method <- method
  first <- c("orthogroup", "method", "status", "reason", "detail", "n_taxa", "n_codons", "taxa_absent", "n_fg_branches",
             "fg_branches", "lrt", "p_value", "p_adj", "p_adj_method", "family_size")
  results <- results[, c(intersect(first, names(results)), setdiff(names(results), first)), drop = FALSE]
  branches <- .phylo_bind_rows(lapply(res, `[[`, "branches"))
  sequences <- .phylo_bind_rows(lapply(res, `[[`, "sequences"))

  out_file <- function(what, ext = "tsv") file.path(run_dir, sprintf("%s_%s.%s", run_name, what, ext))
  write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_tsv(results, out_file("results"))
  write_tsv(branches, out_file("branches"))
  write_tsv(sequences, out_file("sequences"))

  counts <- as.list(table(ifelse(results$status == "tested", "tested", paste0(results$status, ":", results$reason))))
  presence <- vapply(hyp$taxa, function(tx) sum(tested & !vapply(strsplit(ifelse(is.na(results$taxa_absent), "", results$taxa_absent), ","), function(a) tx %in% a, logical(1))), integer(1))
  version_of <- function(pkg) if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  prov <- list(
    analysis = method, dndsR = version_of("dndsR"), r = paste(R.version$major, R.version$minor, sep = "."),
    ape = version_of("ape"), hyphy_version = version, hyphy_validated_version = .hyphy_validated_version,
    hyphy_path = bin, settings = settings, random_seed = seed, hyphy_cpu = hyphy_cpu,
    reproducibility = "Results are reproducible for the same HyPhy build, settings, random_seed, and hyphy_cpu. threads only sets how many orthogroups run at once.",
    min_taxa = min_taxa, terminal_stops = terminal_stops,
    species_tree = list(path = species_tree_file, md5 = unname(tools::md5sum(species_tree_file))),
    foreground = list(path = foreground_file, md5 = unname(tools::md5sum(foreground_file)),
                      hypothesis = hyp$spec, foreground_taxa = hyp$foreground_taxa),
    gene_tree_dir = gene_tree_dir, alignment_dir = alignment_dir,
    multiple_testing = list(method = "BH", family = "orthogroups with status 'tested' in this run, under one foreground hypothesis and one set of settings",
                            family_size = sum(tested)),
    counts = counts,
    tested_orthogroups_per_taxon = as.list(presence),
    ascertainment = "Only single-copy orthogroups that pass QC and keep the foreground identifiable are tested. Excluded orthogroups are listed with reasons; conclusions apply to the tested set only."
  )
  writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), out_file("provenance", "json"))

  message(sprintf("[%s] tested %d, skipped %d, failed %d. Results: %s", method, sum(tested),
                  sum(results$status == "skipped"), sum(results$status == "failed"), out_file("results")))
  other <- results$status != "tested"
  if (any(other)) {
    tab <- table(paste0(results$status[other], ":", results$reason[other]))
    message(paste0("  ", names(tab), ": ", as.integer(tab), collapse = "\n"))
  }
  if (any(results$status == "failed")) {
    warning(sum(results$status == "failed"), " orthogroup(s) failed in HyPhy and are not in the multiple-testing family. See ",
            out_file("results"), call. = FALSE)
  }
  invisible(results)
}
