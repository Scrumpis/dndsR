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
#' @param overwrite Overwrite existing outputs if TRUE.
#' @param verbose Print progress messages if TRUE.
#' @param comparison_file Optional comparison table path or data.frame.
#' @param group_by "gene" or "tx" for how CDS are grouped/named.
#' @param export_proteins If TRUE, also write translated proteins as AA FASTA.
#' @param genetic_code Integer NCBI code for translation (1 = Standard).
#' @param keep_internal_stops If FALSE, drop sequences with internal stops (post-translation).
#' @param cds_suffix Filename suffix for CDS FASTA outputs. Default "_CDS.fasta".
#' @param protein_suffix Filename suffix for AA FASTA outputs. Default "_AA.fasta".
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
                        keep_internal_stops = FALSE,
                        cds_suffix = "_CDS.fasta",
                        protein_suffix = "_AA.fasta") {

  group_by <- match.arg(group_by)

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
  .cds_groups <- function(gff, by = "gene") {
    if (!file.exists(gff)) stop("GFF not found: ", gff, call. = FALSE)
    txdb <- txdbmaker::makeTxDbFromGFF(gff, format = "gff3", circ_seqs = character())
    if (by == "gene") {
      GenomicFeatures::cdsBy(txdb, by = "gene")  # names may exist but are not guaranteed
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

    # detect internal stops: any '*' before the last position
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

    # strip trailing stop '*', preserve names
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

    # ensure required cols exist (subject cols may be NA/blank)
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

    # Ensure FaFile is closed even if something errors
    on.exit(try(close(h$fafile), silent = TRUE), add = TRUE)

    vmsg("  - Parsing GFF and grouping CDS by ", group_by)
    cds_groups <- .cds_groups(gff_path, by = group_by)

    cds_groups2 <- .keep_valid_groups(cds_groups, h$chroms)
    if (length(cds_groups2) == 0L) {
      warning("No CDS groups matched chromosomes in FASTA for: ", out_base)
      return(res)
    }

    vmsg("  - Extracting CDS sequences (", length(cds_groups2), " groups)")
    cds <- GenomicFeatures::extractTranscriptSeqs(h$fafile, cds_groups2)

    # If cdsBy(by="gene") didn’t set names, try to propagate list names
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

  # ---- orchestrate one comparison ----
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

  # ---- batch mode ----
  if (!is.null(comparison_file)) {
    comparisons <- .read_comparisons(comparison_file)

    return(lapply(seq_len(nrow(comparisons)), function(i) {
      row <- comparisons[i, , drop = FALSE]

      # allow blank/NA subjects for query-only
      sf <- row$subject_fasta
      sg <- row$subject_gff
      if (is.na(sf) || !nzchar(sf)) sf <- NULL
      if (is.na(sg) || !nzchar(sg)) sg <- NULL

      run_extraction(
        comparison_name = row$comparison_name,
        query_fasta     = row$query_fasta,
        subject_fasta   = sf,
        query_gff       = row$query_gff,
        subject_gff     = sg
      )
    }))
  }

  # ---- single / single-genome mode ----
  if (is.null(comparison_name) || !nzchar(comparison_name) ||
      is.null(query_fasta) || !nzchar(query_fasta) ||
      is.null(query_gff) || !nzchar(query_gff)) {
    stop("Provide comparison_file or at least comparison_name, query_fasta, and query_gff.",
         call. = FALSE)
  }

  run_extraction(comparison_name, query_fasta, subject_fasta, query_gff, subject_gff)
}
