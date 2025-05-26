#' Calculate dN/dS values using orthologr
#'
#' This function reads a file-of-file-names (fofn) and calculates dN/dS values for
#' each comparison using the \code{orthologr} package.
#'
#' @param file_of_file_names A data.frame or path to a TSV/space-delimited file with columns:
#' comparison name, subject FASTA, query FASTA, subject GFF, query GFF.
#' @param base_dir The base working directory (default: current working directory).
#' @param comp_cores Number of compute cores to use (default: 10).
#' @return TSV files with dN/dS results written to respective comparison folders.
#' @export
calculate_dnds <- function(file_of_file_names,
                                  base_dir = getwd(),
                                  comp_cores = 10) {
  if (is.character(file_of_file_names)) {
    file_of_file_names <- read.table(file_of_file_names, header = FALSE, stringsAsFactors = FALSE)
  }

  # Load required package
  if (!requireNamespace("orthologr", quietly = TRUE)) {
    stop("The orthologr package is required but not installed.")
  }

  for (i in 1:nrow(file_of_file_names)) {
    comparison_basename <- file_of_file_names[i, 1]
    subject_fasta_path <- normalizePath(file_of_file_names[i, 2])
    query_fasta_path   <- normalizePath(file_of_file_names[i, 3])
    subject_gff_path   <- normalizePath(file_of_file_names[i, 4])
    query_gff_path     <- normalizePath(file_of_file_names[i, 5])

    comp_dir <- file.path(base_dir, comparison_basename)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

    query_CDS <- file.path(comp_dir, "query_CDS.fasta")
    subject_CDS <- file.path(comp_dir, "subject_CDS.fasta")

    output_file <- file.path(comp_dir, paste0(comparison_basename, ".tsv"))

    if (file.exists(output_file)) {
      message("Skipping ", comparison_basename, " as output file already exists.")
      next
    }

    message("Running dNdS for: ", comparison_basename)

    result <- orthologr::dNdS(
      query_file = query_CDS,
      subject_file = subject_CDS,
      aligner = "diamond",
      sensitivity_mode = "fast",
      seq_type = "cds",
      format = "fasta",
      ortho_detection = "RBH",
      delete_corrupt_cds = TRUE,
      eval = "1E-5",
      aa_aln_type = "pairwise",
      aa_aln_tool = "NW",
      codon_aln_tool = "pal2nal",
      dnds_est.method = "Comeron",
      comp_cores = comp_cores,
      quiet = TRUE,
      clean_folders = FALSE,
      print_citation = TRUE
    )

    write.table(result, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Finished: ", comparison_basename)
  }

  message("All dN/dS calculations complete.")
}
