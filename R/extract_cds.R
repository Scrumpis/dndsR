#' Extract CDS (and optionally proteins) from genome+GFF using Bioconductor
#'
#' Extracts CDS sequences from a genome FASTA and GFF3 using TxDb/GenomicFeatures.
#' Optionally translates CDS to proteins and writes an amino-acid FASTA.
#'
#' In batch mode, reads a whitespace-delimited comparison file (tabs/spaces) with either:
#' \itemize{
#'   \item headered columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff
#'   \item or headerless with the above in the first 5 columns
#' }
#' Subject inputs may be left blank/NA to run query-only extraction.
#'
#' @param comparison_name Unique identifier for a comparison (e.g., "CheAl_v_CheFo").
#' @param query_fasta Path to the query genome FASTA (required unless using comparison_file).
#' @param subject_fasta Path to the subject genome FASTA (optional in single-genome mode).
#' @param query_gff Path to the query GFF3.
#' @param subject_gff Path to the subject GFF3 (optional in single-genome mode).
#' @param output_dir Output directory; one subdir per comparison.
#' @param overwrite Overwrite existing outputs if TRUE (default FALSE).
#' @param verbose Print progress messages if TRUE (default FALSE).
#' @param threads Integer; maximum number of parallel workers (forked processes; default 1).
#'   In batch mode, workers are applied across unique genome extraction tasks.
#'   On Windows, runs sequentially (no forking).
#' @param comparison_file Optional comparison table path or data.frame.
#' @param group_by "gene" or "tx" for how CDS are grouped/named.
#' @param export_proteins If TRUE, also write translated proteins as AA FASTA.
#' @param genetic_code Integer NCBI code for translation (1 = Standard).
#' @param keep_internal_stops If FALSE, drop sequences with internal stops (post-translation).
#' @param cds_suffix Filename suffix for CDS FASTA outputs. Default "_CDS.fasta".
#' @param protein_suffix Filename suffix for AA FASTA outputs. Default "_AA.fasta".
#'
#' @details
#' TxDb genome build metadata is not required for CDS extraction and is intentionally ignored.
#'
#' @return A list (or list-of-lists in batch) with paths to CDS/protein files and metadata.
#' @export
extract_cds <- function(comparison_name = NULL,
                        query_fasta = NULL,
                        subject_fasta = NULL,
                        query_gff = NULL,
                        subject_gff = NULL,
                        output_dir = ".",
                        overwrite = FALSE,
                        verbose = FALSE,
                        threads = 1,
                        comparison_file = NULL,
                        group_by = c("gene", "tx"),
                        export_proteins = FALSE,
                        genetic_code = 1L,
                        keep_internal_stops = FALSE,
                        cds_suffix = "_CDS.fasta",
                        protein_suffix = "_AA.fasta") {

  group_by <- match.arg(group_by)

  threads <- suppressWarnings(as.integer(threads))
  if (is.na(threads) || threads < 1L) threads <- 1L

  # ---- dependency guards with clear errors ----
  .need_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Install it (or run inside the dndsR container).",
           call. = FALSE)
    }
  }
  .need_pkg("GenomicFeatures")
  .need_pkg("txdbmaker")
  .need_pkg("Biostrings")
  .need_pkg("Rsamtools")
  .need_pkg("GenomeInfoDb")
  .need_pkg("tools")
  .need_pkg("parallel")

  vmsg <- function(...) if (isTRUE(verbose)) message(...)

  # ---- helper: ensure FASTA is indexed; return FaFile and available seqnames ----
  .open_indexed_fasta <- function(fa) {
    if (!file.exists(fa)) stop("FASTA not found: ", fa, call. = FALSE)
    fai <- paste0(fa, ".fai")
    if (!file.exists(fai)) Rsamtools::indexFa(fa)
    ff <- Rsamtools::FaFile(fa)
    idx <- Rsamtools::scanFaIndex(ff)
    list(fafile = ff, chroms = GenomeInfoDb::seqnames(idx))
  }

  # ---- helper: build TxDb and group CDS ----
  # - Always muffle the known harmless genome-metadata warning.
  # - If verbose=TRUE, print it immediately as an informational note.
  .cds_groups <- function(gff, by = "gene") {
    if (!file.exists(gff)) stop("GFF not found: ", gff, call. = FALSE)

    meta_pat <- "genome version information is not available for this TxDb object"

    txdb <- withCallingHandlers(
      txdbmaker::makeTxDbFromGFF(gff, format = "gff3", circ_seqs = character()),
      warning = function(w) {
        msg <- conditionMessage(w)
        is_meta <- grepl(meta_pat, msg, fixed = TRUE)

        if (is_meta) {
          if (isTRUE(verbose)) vmsg("  [TxDb note] ", msg)
          invokeRestart("muffleWarning")
        }
      }
    )

    if (isTRUE(verbose)) {
      gi <- try(GenomeInfoDb::genome(txdb), silent = TRUE)
      if (inherits(gi, "try-error") || all(is.na(gi)) || !any(nzchar(trimws(gi)))) {
        vmsg("  - TxDb genome metadata not provided by GFF; safe to ignore for CDS extraction")
      }
    }

    if (by == "gene") {
      GenomicFeatures::cdsBy(txdb, by = "gene")
    } else {
      GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
    }
  }

  # ---- helper: filter out groups whose seqlevels aren't in FASTA ----
  .keep_valid_groups <- function(cds_groups, genome_chroms) {
    keep <- vapply(cds_groups, function(gr) {
      if (length(gr) == 0) return(FALSE)
      all(as.character(GenomeInfoDb::seqnames(gr)) %in% as.character(genome_chroms))
    }, logical(1))
    cds_groups[keep]
  }

  # ---- helper: translate CDS -> AA with options ----
  .translate_cds <- function(dna, genetic_code = 1L, keep_internal_stops = FALSE) {
    aa <- Biostrings::translate(dna, if.fuzzy.codon = "X", genetic.code = genetic_code)

    aa_chr <- as.character(aa)
    has_internal_stop <- vapply(aa_chr, function(s) {
      n <- nchar(s)
      if (n <= 1) return(FALSE)
      grepl("*", substr(s, 1, n - 1), fixed = TRUE)
    }, logical(1))

    dropped <- sum(has_internal_stop)

    if (!keep_internal_stops && dropped > 0) {
      aa <- aa[!has_internal_stop]
      aa_chr <- aa_chr[!has_internal_stop]
    }

    aa_chr <- sub("\\*$", "", aa_chr)
    aa_out <- Biostrings::AAStringSet(aa_chr)
    names(aa_out) <- names(aa)
    attr(aa_out, "dropped_with_internal_stops") <- dropped
    aa_out
  }

  # ---- read comparison_file (whitespace-delimited, header or not) ----
  .read_comparisons <- function(x) {
    req <- c("comparison_name","query_fasta","query_gff","subject_fasta","subject_gff")

    if (is.data.frame(x)) {
      df <- x
    } else {
      read_ws <- function(hdr) {
        utils::read.table(
          x,
          header = hdr,
          sep = "",
          quote = "\"",
          stringsAsFactors = FALSE,
          comment.char = "",
          strip.white = TRUE,
          blank.lines.skip = TRUE,
          check.names = FALSE
        )
      }
      df1 <- try(read_ws(TRUE), silent = TRUE)
      if (!inherits(df1, "try-error") && all(req %in% names(df1))) {
        df <- df1
      } else {
        df2 <- read_ws(FALSE)
        if (ncol(df2) < 5) {
          stop("comparison_file must have 5 columns or a header with: ", paste(req, collapse = ", "),
               call. = FALSE)
        }
        names(df2)[1:5] <- req
        df <- df2
      }
    }

    miss <- setdiff(req, names(df))
    if (length(miss)) {
      stop("comparison_file missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    df[, req, drop = FALSE]
  }

  # ---- core extraction for one genome ----
  .extract_one <- function(fasta_path, gff_path, out_base) {
    res <- list(
      cds_path = paste0(out_base, cds_suffix),
      prot_path = paste0(out_base, protein_suffix),
      cds_written = FALSE,
      prot_written = FALSE,
      dropped_with_internal_stops = 0L
    )

    if (!overwrite && file.exists(res$cds_path) &&
        (!export_proteins || file.exists(res$prot_path))) {
      vmsg("[OK] Existing outputs found for ", basename(out_base),
           " (use overwrite=TRUE to regenerate).")
      return(res)
    }

    vmsg("  - Indexing/opening FASTA: ", basename(fasta_path))
    h <- .open_indexed_fasta(fasta_path)
    on.exit(try(close(h$fafile), silent = TRUE), add = TRUE)

    vmsg("  - Parsing GFF and grouping CDS by ", group_by)
    cds_groups <- .cds_groups(gff_path, by = group_by)

    cds_groups2 <- .keep_valid_groups(cds_groups, h$chroms)
    if (length(cds_groups2) == 0L) {
      warning("No CDS groups matched chromosomes in FASTA for: ", out_base, call. = FALSE)
      return(res)
    }

    vmsg("  - Extracting CDS sequences (", length(cds_groups2), " groups)")
    cds <- GenomicFeatures::extractTranscriptSeqs(h$fafile, cds_groups2)

    if (is.null(names(cds)) || all(!nzchar(names(cds)))) {
      nm <- names(cds_groups2)
      if (!is.null(nm) && length(nm) == length(cds) && any(nzchar(nm))) {
        names(cds) <- nm
      }
    }

    Biostrings::writeXStringSet(cds, filepath = res$cds_path, compress = FALSE)
    res$cds_written <- TRUE

    if (export_proteins) {
      vmsg("  - Translating CDS -> proteins (genetic code ", genetic_code, ")")
      aa <- .translate_cds(cds, genetic_code = genetic_code, keep_internal_stops = keep_internal_stops)
      res$dropped_with_internal_stops <- attr(aa, "dropped_with_internal_stops")
      Biostrings::writeXStringSet(aa, filepath = res$prot_path, compress = FALSE)
      res$prot_written <- TRUE
    }

    res
  }

  # ---- orchestrate one comparison (single mode only; batch uses unique tasks) ----
  run_extraction <- function(comparison_name, query_fasta, subject_fasta, query_gff, subject_gff) {
    comp_dir <- file.path(output_dir, comparison_name)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

    q_base <- file.path(comp_dir, tools::file_path_sans_ext(basename(query_fasta)))
    s_base <- if (!is.null(subject_fasta) && nzchar(subject_fasta)) {
      file.path(comp_dir, tools::file_path_sans_ext(basename(subject_fasta)))
    } else {
      NULL
    }

    vmsg("== ", comparison_name, " ==")
    q_res <- .extract_one(query_fasta, query_gff, q_base)

    s_res <- NULL
    has_subject <- !is.null(subject_fasta) && nzchar(subject_fasta) &&
      !is.null(subject_gff) && nzchar(subject_gff)

    if (has_subject) {
      s_res <- .extract_one(subject_fasta, subject_gff, s_base)
    } else {
      vmsg("  - Subject not provided: running in single-genome mode")
    }

    list(
      comparison_name = comparison_name,
      query = q_res,
      subject = s_res,
      group_by = group_by,
      export_proteins = export_proteins,
      cds_suffix = cds_suffix,
      protein_suffix = protein_suffix
    )
  }

  # ---- batch vs single ----
  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)

    # Build unique genome extraction tasks
    q <- data.frame(
      comparison_name = df$comparison_name,
      role = "query",
      fasta = df$query_fasta,
      gff = df$query_gff,
      stringsAsFactors = FALSE
    )
    s <- data.frame(
      comparison_name = df$comparison_name,
      role = "subject",
      fasta = df$subject_fasta,
      gff = df$subject_gff,
      stringsAsFactors = FALSE
    )
    tasks <- rbind(q, s)

    tasks <- tasks[!is.na(tasks$fasta) & nzchar(tasks$fasta) &
                     !is.na(tasks$gff) & nzchar(tasks$gff), , drop = FALSE]

    tasks$comp_dir <- file.path(output_dir, tasks$comparison_name)
    tasks$out_base <- file.path(tasks$comp_dir, tools::file_path_sans_ext(basename(tasks$fasta)))
    tasks$cds_path <- paste0(tasks$out_base, cds_suffix)
    tasks$prot_path <- paste0(tasks$out_base, protein_suffix)

    # Deduplicate by the actual output file target
    dedup_key <- if (isTRUE(export_proteins)) tasks$prot_path else tasks$cds_path
    tasks_unique <- tasks[!duplicated(dedup_key), , drop = FALSE]

    message(sprintf("[extract_cds] threads=%d (pid=%d)", threads, Sys.getpid()))
    if (threads > 1L && .Platform$OS.type == "windows") {
      warning("[extract_cds] threads>1 requested on Windows; running sequentially.", call. = FALSE)
      threads <- 1L
    }

    # Ensure comp dirs exist
    for (d in unique(tasks_unique$comp_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

    run_one_task <- function(k) {
      t <- tasks_unique[k, , drop = FALSE]
      comp <- t$comparison_name
      comp_dir <- t$comp_dir

      out_stub <- basename(t$out_base)
      log_file <- file.path(comp_dir, sprintf("%s_%s_extract_cds.log", comp, out_stub))

      # Clean, non-jumbled progress line in the main console:
      message(sprintf("[extract_cds] %s %s (logging to %s/%s_%s_extract_cds.log)",
                      comp, t$role, comp_dir, comp, out_stub))

      .dndsr_with_log(
        log_file = log_file,
        tag = "extract_cds",
        header = c(
          sprintf("[extract_cds] comp=%s role=%s", comp, t$role),
          sprintf("[extract_cds] fasta=%s", basename(t$fasta)),
          sprintf("[extract_cds] gff=%s", basename(t$gff)),
          sprintf("[extract_cds] pid=%d threads=%d", Sys.getpid(), threads)
        ),
        expr = {
          .extract_one(t$fasta, t$gff, t$out_base)
        }
      )
    }

    if (threads > 1L && nrow(tasks_unique) > 1L && .Platform$OS.type != "windows") {
      res_list <- parallel::mclapply(
        seq_len(nrow(tasks_unique)),
        run_one_task,
        mc.cores = min(threads, nrow(tasks_unique)),
        mc.preschedule = FALSE
      )
    } else {
      if (threads > 1L && .Platform$OS.type == "windows") {
        warning("[extract_cds] threads>1 requested on Windows; running sequentially.", call. = FALSE)
      }
      res_list <- lapply(seq_len(nrow(tasks_unique)), run_one_task)
    }

    # Map results by out_base for fast per-comparison assembly
    res_by_outbase <- setNames(res_list, tasks_unique$out_base)

    out <- lapply(seq_len(nrow(df)), function(i) {
      row <- df[i, , drop = FALSE]
      comp <- row$comparison_name
      comp_dir <- file.path(output_dir, comp)

      q_base <- file.path(comp_dir, tools::file_path_sans_ext(basename(row$query_fasta)))
      q_res <- res_by_outbase[[q_base]]
      if (is.null(q_res)) {
        q_res <- list(
          cds_path = paste0(q_base, cds_suffix),
          prot_path = paste0(q_base, protein_suffix),
          cds_written = FALSE,
          prot_written = FALSE,
          dropped_with_internal_stops = 0L
        )
      }

      sf <- row$subject_fasta
      sg <- row$subject_gff
      has_subject <- !is.na(sf) && nzchar(sf) && !is.na(sg) && nzchar(sg)

      s_res <- NULL
      if (has_subject) {
        s_base <- file.path(comp_dir, tools::file_path_sans_ext(basename(sf)))
        s_res <- res_by_outbase[[s_base]]
        if (is.null(s_res)) {
          s_res <- list(
            cds_path = paste0(s_base, cds_suffix),
            prot_path = paste0(s_base, protein_suffix),
            cds_written = FALSE,
            prot_written = FALSE,
            dropped_with_internal_stops = 0L
          )
        }
      }

      list(
        comparison_name = comp,
        query = q_res,
        subject = s_res,
        group_by = group_by,
        export_proteins = export_proteins,
        cds_suffix = cds_suffix,
        protein_suffix = protein_suffix
      )
    })

    message("All extract_cds tasks complete.")
    return(out)
  }

  # ---- single mode ----
  if (is.null(comparison_name) || !nzchar(comparison_name) ||
      is.null(query_fasta) || !nzchar(query_fasta) ||
      is.null(query_gff) || !nzchar(query_gff)) {
    stop("Provide comparison_file or at least comparison_name, query_fasta, and query_gff.",
         call. = FALSE)
  }

  comp_dir <- file.path(output_dir, comparison_name)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  out_stub <- tools::file_path_sans_ext(basename(query_fasta))
  log_file <- file.path(comp_dir, sprintf("%s_%s_extract_cds.log", comparison_name, out_stub))

  message(sprintf("[extract_cds] %s (logging to %s/%s_%s_extract_cds.log)",
                  comparison_name, comp_dir, comparison_name, out_stub))

  .dndsr_with_log(
    log_file = log_file,
    tag = "extract_cds",
    header = c(
      sprintf("[extract_cds] comp=%s", comparison_name),
      sprintf("[extract_cds] pid=%d threads=%d", Sys.getpid(), threads)
    ),
    expr = {
      run_extraction(comparison_name, query_fasta, subject_fasta, query_gff, subject_gff)
    }
  )
}
