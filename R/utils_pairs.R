#' Identity of an unordered genome pair in a comparison file
#'
#' Two rows name the same pair when their query and subject FASTA files are the same
#' two files, in either order. Paths are normalised so a relative and an absolute
#' spelling of one file compare equal. Used by `calculate_dnds(unique_pairs = TRUE)`
#' to run each unordered pair once in a multi-species comparison file.
#'
#' @param row One row of a comparison file, as a list with `query_fasta` and
#'   `subject_fasta`.
#' @return A single string, equal for both orders of a pair, or `NA` when the row has
#'   no FASTA pair (the `seq_list_file` mode).
#' @keywords internal
.dnds_pair_identity <- function(row) {
  q <- row[["query_fasta"]]; s <- row[["subject_fasta"]]
  if (is.null(q) || is.null(s)) return(NA_character_)
  # normalizePath(mustWork = FALSE) leaves a missing relative path relative, so
  # resolve those against the working directory by hand.
  abs_path <- function(x) {
    p <- normalizePath(x, winslash = "/", mustWork = FALSE)
    if (!grepl("^(/|[A-Za-z]:/)", p)) p <- file.path(getwd(), p)
    p
  }
  q <- as.character(q)[1]; s <- as.character(s)[1]
  if (!nzchar(trimws(q)) || !nzchar(trimws(s))) return(NA_character_)
  q <- abs_path(q)
  s <- abs_path(s)
  if (!nzchar(q) || !nzchar(s)) return(NA_character_)
  paste(sort(c(q, s)), collapse = "\r")
}
