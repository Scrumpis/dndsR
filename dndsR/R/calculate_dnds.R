#' Calculate dN/dS values using orthologr
#'
#' This function calculates dN/dS values using the \code{orthologr} package.
#' It can be run in batch mode (via a file-of-file-names) or in single mode
#' using individual file paths.
#'
#' @param file_of_file_names A data.frame or path to a TSV/space-delimited file
#'        with columns: comparison_name, subject FASTA, query FASTA, subject GFF, query GFF.
#'        If NULL, must provide individual file paths below.
#' @param comparison_name Name for the comparison folder/output file (for single mode).
#' @param query_fasta,query_gff,subject_fasta,subject_gff Individual file paths (for single mode).
#' @param output_dir Base working directory (default: current working directory).
#' @param comp_cores Number of compute cores to use (default: 10).
#' @param aligner Alignment tool (default: "diamond").
#' @param sensitivity_mode Sensitivity mode for alignment (default: "fast").
#' @param dnds_method Method for dN/dS calculation (default: "Comeron").
#' @param ... Additional parameters passed to orthologr::dNdS().
#'
#' @return Writes TSV files with dN/dS results to respective comparison folders.
#' @export
calculate_dnds <- function(file_of_file_names = NULL,
                           comparison_name = NULL,
                           subject_fasta = NULL,
                           query_fasta = NULL,
                           subject_gff = NULL,
                           query_gff = NULL,
                           output_dir = getwd(),
                           comp_cores = 4,
                           aligner = "diamond",
                           sensitivity_mode = "fast",
                           dnds_method = "Comeron",
                           ...) {
  if (!requireNamespace("orthologr", quietly = TRUE)) {
    stop("The orthologr package is required but not installed.")
  }

  run_dnds <- function(comparison_basename, subject_fa, query_fa, subject_gff, query_gff) {
    comp_dir <- file.path(output_dir, comparison_basename)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

    # Derive CDS filenames based on base names of FASTA files
    query_base <- tools::file_path_sans_ext(basename(query_fa))
    subject_base <- tools::file_path_sans_ext(basename(subject_fa))

    query_CDS <- file.path(comp_dir, paste0(query_base, "_CDS.fasta"))
    subject_CDS <- file.path(comp_dir, paste0(subject_base, "_CDS.fasta"))
    output_file <- file.path(comp_dir, paste0(comparison_basename, "_dnds.tsv"))

    if (!file.exists(query_CDS) || !file.exists(subject_CDS)) {
      warning("CDS files not found in ", comp_dir, ". Expected: ",
              basename(query_CDS), " and ", basename(subject_CDS),
              ". Skipping ", comparison_basename, ".")
      return(NULL)
    }

    if (file.exists(output_file)) {
      message("Skipping ", comparison_basename, " as output file already exists.")
      return(NULL)
    }

    message("Running dNdS for: ", comparison_basename)

    result <- orthologr::dNdS(
      query_file = query_CDS,
      subject_file = subject_CDS,
      aligner = aligner,
      sensitivity_mode = sensitivity_mode,
      seq_type = "cds",
      format = "fasta",
      ortho_detection = "RBH",
      delete_corrupt_cds = TRUE,
      eval = "1E-5",
      aa_aln_type = "pairwise",
      aa_aln_tool = "NW",
      codon_aln_tool = "pal2nal",
      dnds_est.method = dnds_method,
      comp_cores = comp_cores,
      quiet = TRUE,
      clean_folders = FALSE,
      print_citation = FALSE,
      ...
    )

    write.table(result, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Finished: ", comparison_basename)
  }

  # Batch mode
  if (!is.null(file_of_file_names)) {
    if (is.character(file_of_file_names)) {
      file_of_file_names <- read.table(file_of_file_names, header = FALSE, stringsAsFactors = FALSE)
    }

    for (i in seq_len(nrow(file_of_file_names))) {
      run_dnds(
        comparison_basename = file_of_file_names[i, 1],
        query_fa = normalizePath(file_of_file_names[i, 2]),
        query_gff   = normalizePath(file_of_file_names[i, 3]),
        subject_fa = normalizePath(file_of_file_names[i, 4]),
        subject_gff   = normalizePath(file_of_file_names[i, 5])
      )
    }
  } else {
    # Single mode
    if (any(sapply(list(comparison_name, subject_fasta, query_fasta, subject_gff, query_gff), is.null))) {
      stop("In single mode, please supply comparison_name, subject_fasta, query_fasta, subject_gff, and query_gff.")
    }

    run_dnds(
      comparison_basename = comparison_name,
      subject_fa = normalizePath(subject_fasta),
      query_fa   = normalizePath(query_fasta),
      subject_gff = normalizePath(subject_gff),
      query_gff   = normalizePath(query_gff)
    )
  }

  message("All dN/dS calculations complete.")
}
