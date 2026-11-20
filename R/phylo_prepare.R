#' Prepare codon alignments and trees for phylogenetic selection tests
#'
#' Builds everything [phylo_busted()] needs from genomes or coding sequences and an
#' orthogroup table: one codon alignment per single-copy orthogroup, aligned as protein with
#' MAFFT and translated back to codons; the species tree pruned to each orthogroup with its
#' foreground labelled; normalized gene trees, if given; QC; and a record of every
#' orthogroup and member that was not used, with the reason. Then run
#' `phylo_busted(prepared_dir = ...)`.
#'
#' Preparation does not change the hypothesis or repair data. Each orthogroup is kept or
#' excluded as a whole; a terminal stop codon is removed and counted (or, with
#' `terminal_stops = "exclude"`, excludes the orthogroup); nothing else is edited.
#'
#' @section Coding sequences:
#' Give one of:
#' \itemize{
#'   \item \code{genome_file}: tab-separated with columns \code{taxon}, \code{fasta},
#'     \code{gff}, and optionally \code{orthofinder_name}. CDS are extracted with
#'     [extract_cds()] (longest valid isoform per gene, named by gene ID) into
#'     \code{cds/<taxon>/}. Genes it excludes keep its reason (\code{cds_invalid}).
#'   \item \code{cds_file}: tab-separated with columns \code{taxon}, \code{cds} (a CDS
#'     FASTA named by gene ID), and optionally \code{orthofinder_name}.
#' }
#' Relative paths are read from the current folder, or failing that from the table's folder.
#'
#' @section Orthogroups:
#' \itemize{
#'   \item Long format (\code{orthogroup_format = "long"}): columns \code{orthogroup},
#'     \code{taxon}, \code{gene_id}, one row per gene.
#'   \item OrthoFinder (\code{"orthofinder"}): \code{Orthogroups.tsv}, or
#'     \code{Phylogenetic_Hierarchical_Orthogroups/N0.tsv} (orthogroups named by
#'     \code{HOG}). Species columns are named by proteome; they must equal a taxon or its
#'     \code{orthofinder_name}. Genes must be named as in the CDS, as when OrthoFinder is
#'     run on the proteins from \code{extract_cds --export-proteins}.
#' }
#' An orthogroup with more than one gene from a taxon is excluded (\code{multi_copy}). The
#' member table keeps every gene, so duplication-aware analyses can use it later.
#'
#' @section Gene trees:
#' \code{gene_tree_dir} holds \code{<orthogroup>.nwk} (or OrthoFinder's
#' \code{<orthogroup>_tree.txt}), with tips named by taxon, gene ID, \code{taxon_geneID},
#' or \code{<orthofinder_name>_<geneID>}. They are checked against the orthogroup's genes,
#' compared with the pruned species tree, and written with tips named by taxon.
#'
#' @param orthogroup_file Orthogroup table (see Orthogroups).
#' @param species_tree_file Rooted Newick species tree.
#' @param foreground_file Foreground hypothesis, as for [phylo_busted()].
#' @param genome_file,cds_file Genome or CDS table (see Coding sequences). Give one.
#' @param gene_tree_dir Optional folder of gene trees.
#' @param orthogroup_format "auto" (default), "long", or "orthofinder".
#' @param output_dir Folder for results.
#' @param run_name Name of this preparation, used for its folder.
#' @param min_taxa Fewest taxa an orthogroup may have (default 4).
#' @param terminal_stops "remove" (default) removes a final stop codon and counts it;
#'   "exclude" excludes the orthogroup instead.
#' @param mafft_strategy "linsi" (default; L-INS-i, accurate for few sequences) or "auto".
#' @param gene_tree_mode "topology" (default): gene trees will be given to HyPhy.
#'   "report": the species tree is used, and gene trees are only compared with it.
#' @param threads Orthogroups (and genomes) to process at once.
#' @param overwrite Re-extract CDS even if they exist.
#' @param mafft Path to the MAFFT executable, if not on \code{PATH}.
#' @return Invisibly, the summary table, one row per orthogroup.
#' @seealso [phylo_busted()], `docs/phylo-methods.md`
#' @export
phylo_prepare <- function(orthogroup_file,
                          species_tree_file,
                          foreground_file,
                          genome_file = NULL,
                          cds_file = NULL,
                          gene_tree_dir = NULL,
                          orthogroup_format = c("auto", "long", "orthofinder"),
                          output_dir = ".",
                          run_name = "prepared",
                          min_taxa = 4L,
                          terminal_stops = c("remove", "exclude"),
                          mafft_strategy = c("linsi", "auto"),
                          gene_tree_mode = c("topology", "report"),
                          threads = 1L,
                          overwrite = FALSE,
                          mafft = NULL) {
  .need_ape()
  orthogroup_format <- match.arg(orthogroup_format)
  terminal_stops <- match.arg(terminal_stops)
  mafft_strategy <- match.arg(mafft_strategy)
  gene_tree_mode <- match.arg(gene_tree_mode)
  threads <- max(1L, as.integer(threads))
  min_taxa <- as.integer(min_taxa)
  if (is.na(min_taxa) || min_taxa < 3L) stop("min_taxa must be at least 3.", call. = FALSE)
  if (!grepl("^[A-Za-z0-9_.-]+$", run_name)) stop("run_name may contain only letters, digits, '.', '_' and '-'.", call. = FALSE)
  if (is.null(genome_file) == is.null(cds_file)) stop("Give exactly one of genome_file or cds_file.", call. = FALSE)
  if (!is.null(gene_tree_dir) && !dir.exists(gene_tree_dir)) stop("Gene tree folder not found: ", gene_tree_dir, call. = FALSE)
  mafft_bin <- .mafft_path(mafft)
  mafft_version <- .mafft_version(mafft_bin)

  tree <- .phylo_read_species_tree(species_tree_file)
  hyp <- .phylo_read_foreground(foreground_file, tree, allow_reference = TRUE) # reference rows serve phylo_relax
  run_dir <- file.path(output_dir, run_name)
  for (d in c("alignments", "proteins", "trees")) dir.create(file.path(run_dir, d), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(gene_tree_dir)) dir.create(file.path(run_dir, "gene_trees"), showWarnings = FALSE)

  table_path <- if (is.null(genome_file)) cds_file else genome_file
  sources <- .phylo_read_source_table(table_path, if (is.null(genome_file)) c("taxon", "cds") else c("taxon", "fasta", "gff"))
  unknown <- setdiff(sources$taxon, hyp$taxa)
  if (length(unknown)) stop("Taxa in ", table_path, " are not in the species tree: ", paste(unknown, collapse = ", "), call. = FALSE)
  aliases <- stats::setNames(sources$taxon, sources$orthofinder_name)

  if (!is.null(genome_file)) {
    message(sprintf("[phylo_prepare] extracting CDS from %d genomes", nrow(sources)))
    cds_dir <- file.path(run_dir, "cds")
    .parallel_lapply(seq_len(nrow(sources)), function(i) {
      suppressMessages(extract_cds(comparison_name = sources$taxon[i], query_fasta = sources$fasta[i],
                                   query_gff = sources$gff[i], output_dir = cds_dir, overwrite = overwrite))
      NULL
    }, threads = threads, what = "phylo_prepare")
    stem <- tools::file_path_sans_ext(basename(sources$fasta))
    sources$cds <- file.path(cds_dir, sources$taxon, paste0(stem, "_CDS.fasta"))
    sources$qc <- file.path(cds_dir, sources$taxon, paste0(stem, "_CDS_qc.tsv"))
  } else {
    sources$qc <- NA_character_
  }
  cds <- lapply(seq_len(nrow(sources)), function(i) .phylo_read_cds(sources$cds[i], sources$qc[i]))
  names(cds) <- sources$taxon

  members <- .phylo_read_orthogroups(orthogroup_file, orthogroup_format, aliases)
  bad_taxa <- setdiff(unique(members$taxon), sources$taxon)
  if (length(bad_taxa)) stop("Orthogroups name taxa without genomes or CDS: ", paste(bad_taxa, collapse = ", "), call. = FALSE)
  ogs <- unique(members$orthogroup)
  bad <- ogs[!grepl("^[A-Za-z0-9_.-]+$", ogs)]
  if (length(bad)) stop("Orthogroup names may contain only letters, digits, '.', '_' and '-': ", paste(utils::head(bad, 5), collapse = ", "), call. = FALSE)

  ctx <- list(hyp = hyp, cds = cds, run_dir = run_dir, min_taxa = min_taxa, terminal_stops = terminal_stops,
              mafft = mafft_bin, mafft_args = .mafft_args(mafft_strategy), gene_tree_dir = gene_tree_dir,
              gene_tree_mode = gene_tree_mode, aliases = aliases)
  message(sprintf("[phylo_prepare] %d orthogroups, MAFFT %s (%s), %d at a time", length(ogs), mafft_version, mafft_strategy, threads))
  res <- .parallel_lapply(ogs, function(og) .phylo_prepare_one(og, members[members$orthogroup == og, , drop = FALSE], ctx),
                          threads = threads, what = "phylo_prepare")

  summary <- .phylo_bind_rows(lapply(res, `[[`, "row"))
  member_tab <- .phylo_bind_rows(lapply(res, `[[`, "members"))
  sequences <- .phylo_bind_rows(lapply(res, `[[`, "sequences"))
  first <- c("orthogroup", "status", "reason", "exclusion_basis", "detail", "n_taxa", "taxa", "taxa_absent", "n_codons")
  summary <- summary[, c(intersect(first, names(summary)), setdiff(names(summary), first)), drop = FALSE]
  file.copy(species_tree_file, file.path(run_dir, "species_tree.nwk"), overwrite = TRUE)
  file.copy(foreground_file, file.path(run_dir, "foreground.tsv"), overwrite = TRUE)
  out_file <- function(what, ext = "tsv") file.path(run_dir, sprintf("%s_%s.%s", run_name, what, ext))
  write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_tsv(summary, out_file("summary"))
  write_tsv(member_tab, out_file("members"))
  write_tsv(sequences, out_file("sequences"))

  md5 <- function(p) if (is.null(p) || is.na(p)) NA_character_ else unname(tools::md5sum(p))
  version_of <- function(pkg) if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  label <- ifelse(summary$status == "prepared", "prepared", paste0(summary$status, ":", summary$reason))
  prov <- list(
    step = "phylo_prepare", dndsR = version_of("dndsR"), r = paste(R.version$major, R.version$minor, sep = "."),
    ape = version_of("ape"), biostrings = version_of("Biostrings"),
    mafft = list(path = mafft_bin, version = mafft_version, strategy = mafft_strategy, arguments = ctx$mafft_args),
    settings = list(min_taxa = min_taxa, terminal_stops = terminal_stops, gene_tree_mode = gene_tree_mode,
                    orthogroup_format = attr(members, "format"), genetic_code = "Universal"),
    inputs = list(
      orthogroup_file = list(path = orthogroup_file, md5 = md5(orthogroup_file)),
      species_tree = list(path = species_tree_file, md5 = md5(species_tree_file)),
      foreground = list(path = foreground_file, md5 = md5(foreground_file), hypothesis = hyp$spec),
      source_table = list(path = table_path, md5 = md5(table_path), kind = if (is.null(genome_file)) "cds" else "genome"),
      sources = lapply(seq_len(nrow(sources)), function(i) list(
        taxon = sources$taxon[i],
        fasta_md5 = if (!is.null(sources$fasta)) md5(sources$fasta[i]) else NULL,
        gff_md5 = if (!is.null(sources$gff)) md5(sources$gff[i]) else NULL,
        cds_md5 = md5(sources$cds[i]))),
      gene_tree_dir = gene_tree_dir),
    counts = as.list(table(label)),
    exclusions_by_basis = as.list(table(summary$exclusion_basis[summary$status != "prepared"])),
    prepared_orthogroups_per_taxon = as.list(vapply(hyp$taxa, function(tx) {
      sum(summary$status == "prepared" & vapply(strsplit(ifelse(is.na(summary$taxa), "", summary$taxa), ","), function(t) tx %in% t, logical(1)))
    }, integer(1))),
    ascertainment = "Only single-copy orthogroups whose genes pass CDS and alignment checks, with at least min_taxa taxa and an identifiable foreground, are prepared. Every other orthogroup is in the summary with its reason; conclusions apply to the prepared set only."
  )
  writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), out_file("provenance", "json"))
  manifest <- list(format = 1L, run_name = run_name, alignment_dir = "alignments", species_tree_file = "species_tree.nwk",
                   foreground_file = "foreground.tsv", gene_tree_dir = if (is.null(gene_tree_dir)) NULL else "gene_trees",
                   tree_dir = "trees", gene_tree_mode = gene_tree_mode, summary = basename(out_file("summary")),
                   provenance = basename(out_file("provenance", "json")))
  writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"), file.path(run_dir, "prepare_manifest.json"))

  message(sprintf("[phylo_prepare] prepared %d, excluded %d, failed %d. Summary: %s", sum(summary$status == "prepared"),
                  sum(summary$status == "excluded"), sum(summary$status == "failed"), out_file("summary")))
  other <- summary$status != "prepared"
  if (any(other)) {
    tab <- table(sprintf("%s:%s (%s)", summary$status[other], summary$reason[other], summary$exclusion_basis[other]))
    message(paste0("  ", names(tab), ": ", as.integer(tab), collapse = "\n"))
  }
  if (any(summary$status == "failed")) warning(sum(summary$status == "failed"), " orthogroup(s) failed in preparation. See ", out_file("summary"), call. = FALSE)
  invisible(summary)
}

#' Read a genome or CDS table, resolving relative paths
#' @keywords internal
.phylo_read_source_table <- function(path, required) {
  if (!file.exists(path)) stop("Table not found: ", path, call. = FALSE)
  t <- utils::read.delim(path, colClasses = "character", strip.white = TRUE, comment.char = "#")
  missing <- setdiff(required, names(t))
  if (length(missing) || !nrow(t)) stop(path, " needs a header and columns: ", paste(required, collapse = ", "), call. = FALSE)
  if (anyDuplicated(t$taxon)) stop("Each taxon may appear once in ", path, call. = FALSE)
  .phylo_check_tip_names(t$taxon, path)
  resolve <- function(p) ifelse(file.exists(p) | grepl("^/", p), p, file.path(dirname(path), p))
  for (col in intersect(c("fasta", "gff", "cds"), names(t))) {
    t[[col]] <- resolve(t[[col]])
    gone <- !file.exists(t[[col]])
    if (any(gone)) stop("Files in ", path, " not found: ", paste(t[[col]][gone], collapse = ", "), call. = FALSE)
  }
  if (is.null(t$orthofinder_name)) t$orthofinder_name <- t$taxon
  t$orthofinder_name[is.na(t$orthofinder_name) | !nzchar(t$orthofinder_name)] <- t$taxon[is.na(t$orthofinder_name) | !nzchar(t$orthofinder_name)]
  t
}

#' Read a CDS FASTA, and the extract_cds QC table if there is one
#' @keywords internal
.phylo_read_cds <- function(path, qc_path) {
  s <- Biostrings::readDNAStringSet(path)
  ids <- sub("[[:space:]].*$", "", names(s))
  if (anyDuplicated(ids)) stop("Duplicate sequence names in ", path, call. = FALSE)
  qc <- if (!is.na(qc_path) && file.exists(qc_path)) utils::read.delim(qc_path, colClasses = c(id = "character")) else NULL
  list(seq = stats::setNames(toupper(as.character(s)), ids), qc = qc, path = path)
}

#' Read orthogroups into one row per gene
#' @keywords internal
.phylo_read_orthogroups <- function(path, format, aliases) {
  if (!file.exists(path)) stop("Orthogroup file not found: ", path, call. = FALSE)
  t <- utils::read.delim(path, colClasses = "character", check.names = FALSE, na.strings = character(), comment.char = "")
  if (format == "auto") {
    format <- if (all(c("orthogroup", "taxon", "gene_id") %in% names(t))) "long"
    else if (names(t)[1] %in% c("Orthogroup", "HOG")) "orthofinder"
    else stop("Cannot tell the format of ", path, ": use columns orthogroup, taxon, gene_id, or OrthoFinder's Orthogroups.tsv or N0.tsv.", call. = FALSE)
  }
  if (format == "long") {
    if (!all(c("orthogroup", "taxon", "gene_id") %in% names(t))) stop(path, " needs columns orthogroup, taxon, gene_id.", call. = FALSE)
    out <- data.frame(orthogroup = trimws(t$orthogroup), taxon = trimws(t$taxon), gene_id = trimws(t$gene_id), stringsAsFactors = FALSE)
  } else {
    id_col <- names(t)[1]
    if (!id_col %in% c("Orthogroup", "HOG")) stop(path, " is not an OrthoFinder Orthogroups.tsv or N0.tsv.", call. = FALSE)
    species_cols <- setdiff(names(t), c("Orthogroup", "HOG", "OG", "Gene Tree Parent Clade"))
    unknown <- setdiff(species_cols, names(aliases))
    if (length(unknown)) {
      stop("OrthoFinder columns not matched to a taxon: ", paste(unknown, collapse = ", "),
           ". Name taxa the same, or give them in the orthofinder_name column of the genome or CDS table.", call. = FALSE)
    }
    parts <- lapply(species_cols, function(sp) {
      g <- strsplit(t[[sp]], ",", fixed = TRUE)
      n <- lengths(g)
      data.frame(orthogroup = rep(t[[id_col]], n), taxon = rep(unname(aliases[[sp]]), sum(n)),
                 gene_id = trimws(unlist(g)), stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, parts)
  }
  out <- out[nzchar(out$gene_id) & nzchar(out$orthogroup), , drop = FALSE]
  out <- out[order(out$orthogroup, out$taxon, out$gene_id, method = "radix"), , drop = FALSE]
  attr(out, "format") <- format
  out
}

#' Prepare one orthogroup
#' @keywords internal
.phylo_prepare_one <- function(og, m, ctx) {
  copies <- table(m$taxon)
  m$copies <- as.integer(copies[m$taxon])
  m$cds_status <- NA_character_
  m$used <- FALSE
  present <- sort(unique(m$taxon), method = "radix")
  base <- list(orthogroup = og, n_taxa = length(present), taxa = paste(present, collapse = ","),
               taxa_absent = paste(setdiff(ctx$hyp$taxa, present), collapse = ","))
  out <- function(status, reason = NA_character_, detail = NA_character_, extra = list(), sequences = NULL) {
    list(row = c(base, list(status = status, reason = reason, exclusion_basis = .phylo_reason_basis(reason), detail = detail), extra),
         members = cbind(m, stringsAsFactors = FALSE), sequences = sequences)
  }
  if (any(copies > 1L)) {
    multi <- copies[copies > 1L]
    return(out("excluded", "multi_copy", paste("more than one gene from:", paste(sprintf("%s (%d)", names(multi), as.integer(multi)), collapse = ", "))))
  }

  # Coding sequences.
  seqs <- character(nrow(m))
  for (i in seq_len(nrow(m))) {
    src <- ctx$cds[[m$taxon[i]]]
    id <- m$gene_id[i]
    hit <- if (id %in% names(src$seq)) id else {
      alt <- names(src$seq)[.normalize_gff_id(names(src$seq)) == .normalize_gff_id(id)]
      if (length(alt) == 1L) alt else NA_character_
    }
    if (!is.na(hit)) { seqs[i] <- src$seq[[hit]]; m$cds_status[i] <- "found"; next }
    q <- if (!is.null(src$qc)) src$qc[src$qc$id == id | .normalize_gff_id(src$qc$id) == .normalize_gff_id(id), , drop = FALSE] else NULL
    if (!is.null(q) && nrow(q)) {
      why <- if (!isTRUE(as.logical(q$multiple_of_3[1]))) "length not a multiple of 3" else if (isTRUE(as.logical(q$internal_stop[1]))) "in-frame stop codon" else "not written by extract_cds"
      m$cds_status[i] <- paste("invalid:", why)
    } else {
      m$cds_status[i] <- "not found"
    }
  }
  if (any(grepl("^invalid", m$cds_status))) {
    k <- grepl("^invalid", m$cds_status)
    return(out("excluded", "cds_invalid", paste0("invalid CDS, excluded by extract_cds: ",
                                                 paste(sprintf("%s %s (%s)", m$taxon[k], m$gene_id[k], sub("^invalid: ", "", m$cds_status[k])), collapse = "; "))))
  }
  if (any(m$cds_status == "not found")) {
    k <- m$cds_status == "not found"
    return(out("excluded", "gene_not_found", paste("no CDS named:", paste(sprintf("%s %s", m$taxon[k], m$gene_id[k]), collapse = ", "))))
  }

  # Sequence checks. Nothing is changed except a removed final stop codon.
  table <- .phylo_codon_table()
  seq_rows <- list()
  for (i in seq_len(nrow(m))) {
    s <- seqs[i]
    if (grepl("[^ACGTRYKMSWBDHVN]", s)) return(out("excluded", "invalid_character", paste(m$taxon[i], m$gene_id[i], "has characters other than nucleotides and IUPAC codes")))
    if (nchar(s) %% 3L != 0L || nchar(s) == 0L) return(out("excluded", "length_not_multiple_of_3", sprintf("%s %s: %d nt", m$taxon[i], m$gene_id[i], nchar(s))))
    codons <- substring(s, seq.int(1L, nchar(s), 3L), seq.int(3L, nchar(s), 3L))
    aa <- unname(table[codons])
    stop_at <- which(!is.na(aa) & aa == "*")
    final_stop <- length(codons) %in% stop_at
    internal <- setdiff(stop_at, length(codons))
    if (length(internal)) return(out("excluded", "internal_stop", sprintf("%s %s: stop codon at codon %d", m$taxon[i], m$gene_id[i], internal[1])))
    if (final_stop && ctx$terminal_stops == "exclude") return(out("excluded", "terminal_stop", paste(m$taxon[i], m$gene_id[i], "ends in a stop codon")))
    if (final_stop) { codons <- codons[-length(codons)]; seqs[i] <- paste(codons, collapse = "") }
    seq_rows[[i]] <- data.frame(orthogroup = og, taxon = m$taxon[i], gene_id = m$gene_id[i], cds_codons = length(codons),
                                terminal_stop_removed = final_stop, starts_with_atg = identical(codons[1], "ATG"),
                                ambiguous_codons = sum(grepl("[^ACGT]", codons)), stringsAsFactors = FALSE)
  }
  sequences <- do.call(rbind, seq_rows)
  info <- list(terminal_stops_removed = sum(sequences$terminal_stop_removed))
  if (length(present) < ctx$min_taxa) {
    return(out("excluded", "too_few_taxa", sprintf("%d taxa, fewer than min_taxa = %d", length(present), ctx$min_taxa), info, sequences))
  }

  # Trees and the foreground, before aligning.
  seq_df <- data.frame(seq_name = paste0(m$taxon, "|", m$gene_id), taxon = m$taxon, gene_id = m$gene_id, stringsAsFactors = FALSE)
  gene_tree <- NULL
  gt_path <- .phylo_find_gene_tree(ctx$gene_tree_dir, og)
  if (!is.null(gt_path)) {
    if (is.na(gt_path)) return(out("excluded", "gene_tree_missing", paste("no gene tree for", og), info, sequences))
    g <- .phylo_try_gene_tree(gt_path, seq_df, ctx$aliases)
    if (is.null(g$tree)) return(out(g$status, g$reason, g$detail, info, sequences))
    gene_tree <- g$tree
  }
  map <- .phylo_map_foreground(ctx$hyp, present, gene_tree, ctx$gene_tree_mode)
  info <- c(info, map$info)
  if (!map$ok) return(out("excluded", map$reason, map$detail, info, sequences))

  # Protein alignment with MAFFT, then codons placed under each residue.
  aa <- vapply(seqs, function(s) {
    codons <- substring(s, seq.int(1L, nchar(s), 3L), seq.int(3L, nchar(s), 3L))
    a <- unname(table[codons]); a[is.na(a)] <- "X"
    paste(a, collapse = "")
  }, character(1), USE.NAMES = FALSE)
  input <- file.path(ctx$run_dir, "proteins", paste0(og, ".input.faa"))
  aligned_path <- file.path(ctx$run_dir, "proteins", paste0(og, ".faa"))
  writeLines(as.vector(rbind(paste0(">", m$taxon), aa)), input)
  status <- suppressWarnings(system2(ctx$mafft, c(ctx$mafft_args, shQuote(input)), stdout = aligned_path, stderr = FALSE))
  aligned <- if (identical(as.integer(status), 0L) && file.exists(aligned_path)) .phylo_read_fasta_simple(aligned_path) else NULL
  unlink(input)
  if (is.null(aligned) || !setequal(names(aligned), m$taxon)) {
    return(out("failed", "alignment_failed", sprintf("MAFFT exited with status %s or wrote an incomplete alignment", status), info, sequences))
  }
  aligned <- toupper(aligned[m$taxon])
  codon_aln <- character(length(aligned))
  for (i in seq_along(aligned)) {
    if (!identical(gsub("-", "", aligned[[i]], fixed = TRUE), aa[i])) {
      return(out("failed", "backtranslation_mismatch", paste("aligned protein differs from the translation of", m$taxon[i]), info, sequences))
    }
    chars <- strsplit(aligned[[i]], "", fixed = TRUE)[[1]]
    codons <- substring(seqs[i], seq.int(1L, nchar(seqs[i]), 3L), seq.int(3L, nchar(seqs[i]), 3L))
    out_codons <- rep("---", length(chars))
    out_codons[chars != "-"] <- codons
    codon_aln[i] <- paste(out_codons, collapse = "")
  }
  aln <- data.frame(seq_name = seq_df$seq_name, taxon = m$taxon, gene_id = m$gene_id, sequence = codon_aln, stringsAsFactors = FALSE)
  # The codon alignment, without gaps, must be exactly the coding sequences used.
  if (!identical(gsub("-", "", codon_aln, fixed = TRUE), seqs)) {
    return(out("failed", "backtranslation_mismatch", "codon alignment does not reproduce the coding sequences", info, sequences))
  }
  qc <- .phylo_check_alignment(aln, ctx$hyp$taxa, ctx$min_taxa, "exclude")
  if (!qc$ok) return(out("failed", "backtranslation_mismatch", paste("codon alignment failed its own checks:", qc$reason), info, sequences))

  aln_path <- file.path(ctx$run_dir, "alignments", paste0(og, ".fasta"))
  writeLines(as.vector(rbind(paste0(">", aln$seq_name), aln$sequence)), aln_path)
  tree_path <- file.path(ctx$run_dir, "trees", paste0(og, ".nwk"))
  writeLines(map$newick, tree_path)
  if (!is.null(gene_tree)) {
    gt <- gene_tree
    gt$edge.length <- NULL
    gt$node.label <- NULL
    ape::write.tree(gt, file.path(ctx$run_dir, "gene_trees", paste0(og, ".nwk")))
  }
  m$used <- TRUE
  flags <- qc$flags
  if (isTRUE(map$info$fg_absent_branches > 0)) flags <- c(flags, sprintf("foreground_branches_absent:%d", map$info$fg_absent_branches))
  if (isTRUE(map$info$gt_missing_species_splits > 0)) flags <- c(flags, "gene_tree_discordant")
  if (isTRUE(map$info$gt_discordant_near_foreground)) flags <- c(flags, "gene_tree_discordant_near_foreground")
  if (any(!sequences$starts_with_atg)) flags <- c(flags, paste0("no_start_codon:", sequences$taxon[!sequences$starts_with_atg]))
  extra <- c(info, as.list(qc$metrics[setdiff(names(qc$metrics), c("n_taxa", "terminal_stops_masked"))]),
             list(n_fg_branches = sum(map$branches$foreground), qc_flags = paste(flags, collapse = ";"),
                  alignment_md5 = unname(tools::md5sum(aln_path)), tree_md5 = unname(tools::md5sum(tree_path))))
  out("prepared", extra = extra, sequences = sequences)
}

#' Read FASTA into a named character vector (names are the first word of each header)
#' @keywords internal
.phylo_read_fasta_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^>", lines)
  if (!length(starts)) return(NULL)
  ends <- c(starts[-1] - 1L, length(lines))
  seqs <- mapply(function(s, e) if (e > s) paste(gsub("[[:space:]]", "", lines[(s + 1L):e]), collapse = "") else "", starts, ends)
  stats::setNames(unname(seqs), sub("^>[[:space:]]*([^[:space:]]*).*", "\\1", lines[starts]))
}

#' Find MAFFT
#' @keywords internal
.mafft_path <- function(mafft = NULL) {
  if (!is.null(mafft)) {
    if (!file.exists(mafft)) stop("MAFFT not found at ", mafft, call. = FALSE)
    return(normalizePath(mafft))
  }
  candidates <- c(mafft, getOption("dndsR.mafft"), Sys.getenv("DNDSR_MAFFT", ""), Sys.which("mafft"))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  stop("MAFFT was not found. It is in the dndsR image; otherwise install it (conda install -c conda-forge -c bioconda mafft) ",
       "or give its path with --mafft.", call. = FALSE)
}

#' MAFFT version, such as "7.526"
#' @keywords internal
.mafft_version <- function(bin) {
  out <- tryCatch(suppressWarnings(system2(bin, "--version", stdout = TRUE, stderr = TRUE)), error = function(e) character())
  v <- regmatches(out, regexpr("[0-9]+\\.[0-9]+", out))
  if (!length(v)) stop("Could not read the MAFFT version from '", bin, " --version'.", call. = FALSE)
  v[1]
}

#' MAFFT arguments: one thread, so results do not depend on threads
#' @keywords internal
.mafft_args <- function(strategy) {
  switch(strategy,
         linsi = c("--localpair", "--maxiterate", "1000", "--thread", "1", "--quiet"),
         auto = c("--auto", "--thread", "1", "--quiet"))
}

#' Read a phylo_prepare manifest
#' @keywords internal
.phylo_read_manifest <- function(prepared_dir) {
  path <- file.path(prepared_dir, "prepare_manifest.json")
  if (!file.exists(path)) stop("No prepare_manifest.json in ", prepared_dir, ". Run phylo_prepare first.", call. = FALSE)
  m <- jsonlite::read_json(path, simplifyVector = TRUE)
  at <- function(x) if (is.null(x)) NULL else file.path(prepared_dir, x)
  summary <- utils::read.delim(file.path(prepared_dir, m$summary), colClasses = c(orthogroup = "character"), na.strings = "NA")
  prov <- file.path(prepared_dir, m$provenance)
  list(alignment_dir = at(m$alignment_dir), species_tree_file = at(m$species_tree_file), foreground_file = at(m$foreground_file),
       gene_tree_dir = at(m$gene_tree_dir), tree_dir = at(m$tree_dir), gene_tree_mode = m$gene_tree_mode,
       summary = summary, provenance_md5 = unname(tools::md5sum(prov)),
       inputs = jsonlite::read_json(prov, simplifyVector = TRUE)$inputs)
}
