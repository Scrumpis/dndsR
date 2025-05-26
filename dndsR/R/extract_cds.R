#' Extract CDS regions from genome and GFF files using gffread
#'
#' @param comparison_name A unique string identifier for the comparison
#' @param query_fasta Path to the query genome FASTA
#' @param subject_fasta Path to the subject genome FASTA
#' @param query_gff Path to the query GFF3 annotation
#' @param subject_gff Path to the subject GFF3 annotation
#' @param output_dir Directory where output CDS files should be written
#' @param overwrite Logical; whether to overwrite existing output files
#' @param verbose Logical; whether to print progress messages
#'
#' @return A named list with full paths to the query and subject CDS FASTA files
#' @export
extract_cds <- function(comparison_name,
                                query_fasta,
                                subject_fasta,
                                query_gff,
                                subject_gff,
                                output_dir = ".",
                                overwrite = FALSE,
                                verbose = TRUE) {
  # Create subdirectory
  comp_dir <- file.path(output_dir, comparison_name)
  if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)

  # Output files
  query_cds_path <- file.path(comp_dir, "query_CDS.fasta")
  subject_cds_path <- file.path(comp_dir, "subject_CDS.fasta")
  summary_file <- file.path(comp_dir, paste0(comparison_name, ".tsv"))

  if (file.exists(summary_file) && !overwrite) {
    if (verbose) message("Skipping ", comparison_name, " (already processed).")
    return(list(query_cds = query_cds_path,
                subject_cds = subject_cds_path,
                skipped = TRUE))
  }

  # Run gffread for subject
  cmd_subject <- sprintf("gffread -x '%s' -g '%s' '%s'",
                         subject_cds_path,
                         normalizePath(subject_fasta),
                         normalizePath(subject_gff))
  system(cmd_subject)

  # Run gffread for query
  cmd_query <- sprintf("gffread -x '%s' -g '%s' '%s'",
                       query_cds_path,
                       normalizePath(query_fasta),
                       normalizePath(query_gff))
  system(cmd_query)

  if (verbose) {
    message("Finished CDS extraction for: ", comparison_name)
    message("Query CDS: ", query_cds_path)
    message("Subject CDS: ", subject_cds_path)
  }

  return(list(query_cds = query_cds_path,
              subject_cds = subject_cds_path,
              skipped = FALSE))
}
