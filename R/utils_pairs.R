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
  # .abs_path resolves a missing relative path against the working directory, which
  # normalizePath(mustWork = FALSE) does not.
  abs_path <- .abs_path
  q <- as.character(q)[1]; s <- as.character(s)[1]
  if (!nzchar(trimws(q)) || !nzchar(trimws(s))) return(NA_character_)
  q <- abs_path(q)
  s <- abs_path(s)
  if (!nzchar(q) || !nzchar(s)) return(NA_character_)
  paste(sort(c(q, s)), collapse = "\r")
}

#' Read a comparison file
#'
#' The comparison-file reader `calculate_dnds()` has always used, shared so that
#' every command can see the optional `snp_file` and `snp_side` columns. A header
#' row is recognised by its column names; a headerless file is read as
#' `comparison_name query_fasta query_gff subject_fasta subject_gff`, then
#' `snp_file` and `snp_side`.
#'
#' @param x Path to a whitespace-delimited comparison file, or a data.frame.
#' @return A data.frame, with attribute `headerless_mk` set when columns 6-7 of a
#'   headerless file were read as McDonald-Kreitman columns.
#' @keywords internal
.read_comparison_table <- function(x) {

  # ---- helper: read whitespace-delimited ----
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

  # ---- data.frame passthrough ----
  if (is.data.frame(x)) {
    df <- x
  } else {

    # ---- peek first non-empty line to detect header ----
    first_line <- ""
    con <- file(x, open = "r")
    on.exit(close(con), add = TRUE)

    # first non-blank line; the old loop never ran, so a header was read as data
    repeat {
      ln <- readLines(con, n = 1, warn = FALSE)
      if (length(ln) == 0) break
      first_line <- trimws(ln)
      if (nzchar(first_line)) break
    }

    has_header <- FALSE
    if (nzchar(first_line)) {
      toks <- strsplit(first_line, "\\s+")[[1]]
      has_header <- any(toks %in% c(
        "comparison_name",
        "query_fasta", "subject_fasta",
        "query_gff", "subject_gff"
      ))
    }

    df <- read_ws(has_header)
  }

  # ---- ensure comparison_name exists ----
  if (!("comparison_name" %in% names(df))) {
    if (ncol(df) >= 1) {
      names(df)[1] <- "comparison_name"
    }
  }

  if (!("comparison_name" %in% names(df))) {
    stop(
      "comparison_file must include a 'comparison_name' column ",
      "(or be headerless with comparison_name in column 1).",
      call. = FALSE
    )
  }

  # ensure required columnts exist
  needed <- c("comparison_name", "query_fasta", "subject_fasta", "query_gff", "subject_gff")

  # headerless schema inference: comparison_name, query_fasta, query_gff, subject_fasta,
  # subject_gff, and optionally snp_file and snp_side
  headerless_like <- !any(grepl("^(query_|subject_)", names(df)))
  inferred_headerless <- FALSE
  if (!all(needed %in% names(df)) && headerless_like && ncol(df) >= 5) {
    names(df)[1:5] <- c("comparison_name", "query_fasta", "query_gff", "subject_fasta", "subject_gff")
    if (ncol(df) >= 6L) names(df)[6] <- "snp_file"
    if (ncol(df) >= 7L) names(df)[7] <- "snp_side"
    inferred_headerless <- ncol(df) >= 6L
  }

  if (!all(needed %in% names(df))) {
    stop(
      paste0(
        "comparison_file must include columns:\n",
        "* comparison_name, query_fasta, query_gff, subject_fasta, subject_gff\n",
        "For FASTA-only batch execution, use seq_list_file with columns:\n",
        "* comparison_name, query_seq, subject_seq"
      ),
      call. = FALSE
    )
  }

  attr(df, "headerless_mk") <- inferred_headerless
  df
}

#' Path made absolute, for comparing two spellings of one file
#' @keywords internal
.abs_path <- function(x) {
  p <- normalizePath(x, winslash = "/", mustWork = FALSE)
  if (!grepl("^(/|[A-Za-z]:/)", p)) p <- file.path(getwd(), p)
  p
}

#' McDonald-Kreitman inputs of one comparison row
#'
#' @param row One row of a comparison file, as a list.
#' @return A string naming the SNP file and the genome it is on, or `NA` when the
#'   row has no MK inputs.
#' @keywords internal
.mk_inputs_identity <- function(row) {
  f <- row[["snp_file"]]; side <- row[["snp_side"]]
  if (is.null(f) || is.null(side)) return(NA_character_)
  f <- as.character(f)[1]; side <- tolower(trimws(as.character(side)[1]))
  if (is.na(f) || !nzchar(trimws(f))) return(NA_character_)
  genome <- if (identical(side, "query")) row[["query_fasta"]] else if (identical(side, "subject")) row[["subject_fasta"]] else side
  paste(.abs_path(f), .abs_path(as.character(genome)[1]), sep = "\r")
}

#' Comparisons that `unique_pairs = TRUE` skips
#'
#' A row is skipped when an earlier, kept row names the same two genome FASTA
#' files, in either order, and the row adds no McDonald-Kreitman inputs of its own:
#' it has no `snp_file`, or the same SNP file on the same genome as a kept row. A
#' row with its own MK inputs is kept, so its MK table is written, at the cost of a
#' second dN/dS run for that pair. Every batch command applies the same rule to
#' the same file, so later steps skip exactly the rows `calculate_dnds` skipped.
#'
#' @param comparison_file Path or data.frame, as for the batch commands.
#' @param unique_pairs Logical; when not `TRUE`, nothing is skipped.
#' @param command Name used in the message.
#' @return Character vector of skipped comparison names (possibly empty).
#' @keywords internal
.unique_pairs_skip <- function(comparison_file, unique_pairs, command = "dndsR") {
  if (!isTRUE(unique_pairs) || is.null(comparison_file)) return(character(0))
  df <- .read_comparison_table(comparison_file)
  seen <- list()            # pair identity -> MK identities already covered ("" = none)
  skip <- character(0); kept_mk <- character(0)
  for (i in seq_len(nrow(df))) {
    row <- as.list(df[i, , drop = FALSE])
    pair <- .dnds_pair_identity(row)
    if (is.na(pair)) next
    mk <- .mk_inputs_identity(row)
    mk_key <- if (is.na(mk)) "" else mk
    if (is.null(seen[[pair]])) {
      seen[[pair]] <- mk_key
    } else if (mk_key == "" || mk_key %in% seen[[pair]]) {
      skip <- c(skip, as.character(row[["comparison_name"]]))
    } else {
      seen[[pair]] <- c(seen[[pair]], mk_key)
      kept_mk <- c(kept_mk, as.character(row[["comparison_name"]]))
    }
  }
  if (length(skip)) {
    message(command, ": unique_pairs = TRUE: skipping ", length(skip), " comparison(s) whose genome pair ",
            "was already listed, in either order, with no new McDonald-Kreitman inputs: ",
            paste(skip, collapse = ", "))
  }
  if (length(kept_mk)) {
    message(command, ": unique_pairs = TRUE: keeping ", paste(kept_mk, collapse = ", "),
            ", which repeat an earlier genome pair but carry their own McDonald-Kreitman inputs")
  }
  skip
}

#' Drop the rows `.unique_pairs_skip()` named
#' @keywords internal
.drop_comparisons <- function(df, skip) {
  if (!length(skip)) return(df)
  df[!as.character(df$comparison_name) %in% skip, , drop = FALSE]
}
