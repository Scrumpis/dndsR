#' Extract CDS (and optionally proteins) from genome+GFF using Bioconductor
#'
#' @param comparison_name Unique identifier for a comparison (e.g., "CheAl_v_CheFo").
#' @param query_fasta Path to the query genome FASTA (required unless using comparison_file).
#' @param subject_fasta Path to the subject genome FASTA (optional in single-genome mode).
#' @param query_gff Path to the query GFF3.
#' @param subject_gff Path to the subject GFF3 (optional in single-genome mode).
#' @param output_dir Output directory; one subdir per comparison.
#' @param overwrite Overwrite existing outputs if TRUE.
#' @param verbose Print progress.
#' @param comparison_file Optional TSV: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#' @param group_by "gene" or "tx" for how CDS are grouped/named.
#' @param export_proteins If TRUE, also write translated proteins (*.faa).
#' @param genetic_code Integer NCBI code for translation (1=Standard).
#' @param keep_internal_stops If FALSE, drop sequences with internal stops (post-translation).
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
                        verbose = TRUE,
                        comparison_file = NULL,
                        group_by = c("gene", "tx"),
                        export_proteins = FALSE,
                        genetic_code = 1L,
                        keep_internal_stops = FALSE) {

  group_by <- match.arg(group_by)

  # Require namespaces you actually use
  requireNamespace("GenomicFeatures", quietly = TRUE)
  requireNamespace("Biostrings", quietly = TRUE)
  requireNamespace("Rsamtools", quietly = TRUE)
  requireNamespace("GenomeInfoDb", quietly = TRUE)
  requireNamespace("tools", quietly = TRUE)

  # Helper: ensure FASTA is indexed; return open FaFile and seqnames
  .open_indexed_fasta <- function(fa) {
    fai <- paste0(fa, ".fai")
    if (!file.exists(fai)) Rsamtools::indexFa(fa)
    ff <- Rsamtools::FaFile(fa)
    # BiocGenerics::open(ff)   # <- remove
    idx <- Rsamtools::scanFaIndex(ff)
    list(fafile = ff, chroms = GenomeInfoDb::seqnames(idx))
  }

  # Helper: build TxDb and choose CDS grouping
  .cds_groups <- function(gff, by = "gene") {
    txdb <- GenomicFeatures::makeTxDbFromGFF(gff, format = "gff3", circ_seqs = character())
    if (by == "gene") {
      GenomicFeatures::cdsBy(txdb, by = "gene")            # ← no use.names here
    } else {
      GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
    }
  }


  # Helper: filter out groups whose seqlevels aren't in FASTA
  .keep_valid_groups <- function(cds_groups, genome_chroms) {
    keep <- vapply(cds_groups, function(gr) {
      if (length(gr) == 0) return(FALSE)
      all(as.character(GenomeInfoDb::seqnames(gr)) %in% as.character(genome_chroms))
    }, logical(1))
    cds_groups[keep]
  }

  # Helper: translate CDS DNAStringSet to AAStringSet with options
  .translate_cds <- function(dna, genetic_code = 1L, keep_internal_stops = FALSE) {
    aa <- Biostrings::translate(dna, if.fuzzy.codon = "X", genetic.code = genetic_code)
    # identify internal stops
    has_internal_stop <- vapply(as.character(aa), function(s) {
      any(grepl("\\*", substr(s, 1, max(0, nchar(s) - 1)), fixed = TRUE))
    }, logical(1))
    if (!keep_internal_stops && any(has_internal_stop)) {
      aa <- aa[!has_internal_stop]
      attr(aa, "dropped_with_internal_stops") <- sum(has_internal_stop)
    } else {
      attr(aa, "dropped_with_internal_stops") <- 0L
    }
    # strip trailing stop '*'
    aa <- Biostrings::AAStringSet(sub("\\*$", "", as.character(aa)))
    aa
  }

  # Core extraction for one genome (returns list of created/skipped paths)
  .extract_one <- function(fasta_path, gff_path, out_base, overwrite, verbose) {
    res <- list(cds_path = paste0(out_base, "_CDS.fasta"),
                prot_path = paste0(out_base, "_proteins.faa"),
                cds_written = FALSE, prot_written = FALSE,
                dropped_with_internal_stops = 0L)

    if (!overwrite && file.exists(res$cds_path) && (!export_proteins || file.exists(res$prot_path))) {
      if (verbose) message("[OK] Existing outputs found for ", out_base, " (use overwrite=TRUE to regenerate).")
      return(res)
    }

    if (verbose) message("  - Indexing/opening FASTA: ", basename(fasta_path))
    h <- .open_indexed_fasta(fasta_path)
    #on.exit(BiocGenerics::close(h$fafile), add = TRUE)

    if (verbose) message("  - Parsing GFF and grouping CDS by ", group_by)
    cds_groups <- .cds_groups(gff_path, by = group_by)

    # Filter by available chromosomes
    cds_groups2 <- .keep_valid_groups(cds_groups, h$chroms)
    if (length(cds_groups2) == 0L) {
      warning("No CDS groups matched chromosomes in FASTA for: ", out_base)
      return(res)
    }

    if (verbose) message("  - Extracting CDS sequences (", length(cds_groups2), " groups)")
    cds <- GenomicFeatures::extractTranscriptSeqs(h$fafile, cds_groups2)

    Biostrings::writeXStringSet(cds, filepath = res$cds_path, compress = FALSE)
    res$cds_written <- TRUE

    if (export_proteins) {
      if (verbose) message("  - Translating CDS → proteins (genetic code ", genetic_code, ")")
      aa <- .translate_cds(cds, genetic_code = genetic_code, keep_internal_stops = keep_internal_stops)
      res$dropped_with_internal_stops <- attr(aa, "dropped_with_internal_stops")
      Biostrings::writeXStringSet(aa, filepath = res$prot_path, compress = FALSE)
      res$prot_written <- TRUE
    }

    res
  }

  # Orchestrate one comparison (query-only allowed)
  run_extraction <- function(comparison_name, query_fasta, subject_fasta, query_gff, subject_gff) {
    comp_dir <- file.path(output_dir, comparison_name)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

    q_base <- file.path(comp_dir, tools::file_path_sans_ext(basename(query_fasta)))
    s_base <- if (!is.null(subject_fasta)) file.path(comp_dir, tools::file_path_sans_ext(basename(subject_fasta))) else NULL

    if (verbose) message("== ", comparison_name, " ==")
    q_res <- .extract_one(query_fasta, query_gff, q_base, overwrite, verbose)

    s_res <- NULL
    if (!is.null(subject_fasta) && !is.null(subject_gff)) {
      s_res <- .extract_one(subject_fasta, subject_gff, s_base, overwrite, verbose)
    } else if (verbose) {
      message("  - Subject not provided: running in single-genome mode")
    }

    list(
      comparison_name = comparison_name,
      query = q_res,
      subject = s_res,
      group_by = group_by,
      export_proteins = export_proteins
    )
  }

  # Batch mode
  if (!is.null(comparison_file)) {
    comparisons <- read.table(comparison_file,
                              header = FALSE, stringsAsFactors = FALSE, quote = "\"",
                              col.names = c("comparison_name", "query_fasta", "query_gff", "subject_fasta", "subject_gff"))

    required_cols <- c("comparison_name", "query_fasta", "query_gff", "subject_fasta", "subject_gff")
    if (!all(required_cols %in% colnames(comparisons))) {
      stop("comparison_file must contain columns: ", paste(required_cols, collapse = ", "))
    }

    return(lapply(seq_len(nrow(comparisons)), function(i) {
      row <- comparisons[i, ]
      run_extraction(row$comparison_name,
                     row$query_fasta,
                     row$subject_fasta,
                     row$query_gff,
                     row$subject_gff)
    }))
  }

  # Single / single-genome mode
  if (is.null(comparison_name) || is.null(query_fasta) || is.null(query_gff)) {
    stop("Provide a comparison_file or at least comparison_name, query_fasta, and query_gff.")
  }

  return(run_extraction(comparison_name, query_fasta, subject_fasta, query_gff, subject_gff))
}
