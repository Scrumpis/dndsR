#' Extract CDS regions from genome and GFF files using Bioconductor
#'
#' @param comparison_name A unique string identifier for the comparison (sp1_v_sp2; CheAlvCheFo)
#' @param query_fasta Path to the query genome FASTA
#' @param subject_fasta Path to the subject genome FASTA
#' @param query_gff Path to the query GFF3 gene annotation
#' @param subject_gff Path to the subject GFF3 gene annotation
#' @param output_dir Directory where output CDS files should be written (a new directory is created here for each comparison)
#' @param overwrite Logical; whether to overwrite existing output files
#' @param verbose Logical; whether to print progress messages
#' @param comparison_file Optional TSV file with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff
#'
#' @return A list of named lists with paths to query and subject CDS files
#' @export
extract_cds <- function(comparison_name = NULL,
                        query_fasta = NULL,
                        subject_fasta = NULL,
                        query_gff = NULL,
                        subject_gff = NULL,
                        output_dir = ".",
                        overwrite = FALSE,
                        verbose = TRUE,
                        comparison_file = NULL) {

  # Require Bioconductor packages
  requireNamespace("GenomicFeatures", quietly = TRUE)
  requireNamespace("Biostrings", quietly = TRUE)

  # Helper function to extract CDS natively
  extract_native_cds <- function(fasta_path, gff_path, output_path) {
    # Open FASTA as indexed FaFile object
    genome <- Rsamtools::FaFile(fasta_path)
    open(genome)  # Do NOT prefix with Rsamtools::

    # Get chromosome names from FASTA index
    fasta_index <- Rsamtools::scanFaIndex(genome)
    genome_chroms <- GenomeInfoDb::seqnames(fasta_index)

    # Load GFF and create TxDb
    txdb <- GenomicFeatures::makeTxDbFromGFF(gff_path, format = "gff3", circ_seqs = character())
    cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)

    # Identify chromosomes in CDS not present in FASTA
    cds_seqlevels <- unique(as.character(GenomeInfoDb::seqnames(unlist(cds_by_tx))))
    missing_chroms <- setdiff(cds_seqlevels, genome_chroms)

    if (length(missing_chroms) > 0) {
      warning("The following chromosomes are present in the GFF but missing from the FASTA and will be skipped: ",
              paste(missing_chroms, collapse = ", "))

      # Keep only transcript groups with valid chromosomes
      valid_tx <- names(cds_by_tx)[
        vapply(cds_by_tx, function(gr) {
          gr_seqnames <- as.character(GenomeInfoDb::seqnames(gr))
          genome_chroms_chr <- as.character(genome_chroms)  # ensure this is also character
          if (length(gr_seqnames) == 0 || length(genome_chroms_chr) == 0) return(FALSE)
          all(gr_seqnames %in% genome_chroms_chr)
        }, logical(1))
      ]

      cds_by_tx <- cds_by_tx[valid_tx]

      if (length(cds_by_tx) == 0) {
        warning("No CDS regions matched chromosomes in the FASTA. Skipping: ", output_path)
        close(genome)
        return(invisible(NULL))
      }
    }

    # Extract sequences and write to file
    cds_seqs <- GenomicFeatures::extractTranscriptSeqs(genome, cds_by_tx)
    Biostrings::writeXStringSet(cds_seqs, filepath = output_path)

    close(genome)
  }


  # Extract CDS
  run_extraction <- function(comparison_name, query_fasta, subject_fasta, query_gff, subject_gff) {
    comp_dir <- file.path(output_dir, comparison_name)
    if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)

    #query_cds_path <- file.path(comp_dir, "query_CDS.fasta")
    #subject_cds_path <- file.path(comp_dir, "subject_CDS.fasta")
    query_base <- tools::file_path_sans_ext(basename(query_fasta))
    subject_base <- tools::file_path_sans_ext(basename(subject_fasta))

    query_cds_path <- file.path(comp_dir, paste0(query_base, "_CDS.fasta"))
    subject_cds_path <- file.path(comp_dir, paste0(subject_base, "_CDS.fasta"))

    query_skipped <- FALSE
    subject_skipped <- FALSE

    # Extract query CDS if needed
    if (!file.exists(query_cds_path)) {
      if (verbose) message("Extracting query CDS for ", comparison_name)
      extract_native_cds(query_fasta, query_gff, query_cds_path)
    } else {
      if (verbose) message("Query CDS already exists for ", comparison_name, "; skipping.")
      query_skipped <- TRUE
    }

    # Extract subject CDS if needed
    if (!file.exists(subject_cds_path)) {
      if (verbose) message("Extracting subject CDS for ", comparison_name)
      extract_native_cds(subject_fasta, subject_gff, subject_cds_path)
    } else {
      if (verbose) message("Subject CDS already exists for ", comparison_name, "; skipping.")
      subject_skipped <- TRUE
    }

    print("CDS extraction complete")

    return(list(
      query_cds = query_cds_path,
      subject_cds = subject_cds_path,
      query_skipped = query_skipped,
      subject_skipped = subject_skipped,
      comparison_name = comparison_name
    ))
  }

  # Batch mode
  if (!is.null(comparison_file)) {
    comparisons <- read.table(comparison_file,
                              header = FALSE,
                              stringsAsFactors = FALSE,
                              quote = "\"",
                              col.names = c("comparison_name", "query_fasta", "query_gff", "subject_fasta", "subject_gff"))


    required_cols <- c("comparison_name", "query_fasta", "query_gff", "subject_fasta", "subject_gff")
    if (!all(required_cols %in% colnames(comparisons))) {
      stop("The comparison_file must contain columns: ", paste(required_cols, collapse = ", "))
    }

    results <- lapply(seq_len(nrow(comparisons)), function(i) {
      row <- comparisons[i, ]
      run_extraction(row$comparison_name,
                     row$query_fasta,
                     row$subject_fasta,
                     row$query_gff,
                     row$subject_gff)
    })
    return(results)
  }

  # Single mode
  if (is.null(comparison_name) || is.null(query_fasta) || is.null(subject_fasta) ||
      is.null(query_gff) || is.null(subject_gff)) {
    stop("Must provide either a valid comparison_file or all required arguments for a single comparison.")
  }

  return(run_extraction(comparison_name, query_fasta, subject_fasta, query_gff, subject_gff))
}
