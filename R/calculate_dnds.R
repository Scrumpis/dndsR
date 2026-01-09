#' Calculate dN/dS in batch or single mode using orthologr
#'
#' Wraps \code{orthologr::dNdS()} to provide batch execution, sensible defaults,
#' reproducible container-friendly behavior, and standardized output locations.
#'
#' This function supports two input styles:
#' \itemize{
#'   \item \strong{Single mode}: provide \code{query_seq} and \code{subject_seq} paths
#'     (CDS or protein FASTAs depending on \code{seq_type}).
#'   \item \strong{Batch mode}: provide comparison_file.
#' }
#'
#' @param comparison_file Path to a whitespace-delimited (tabs/spaces) file OR a data.frame
#'   describing comparisons (batch mode).
#'   Must include columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#'   Additional columns are allowed and ignored.
#'
#'   Notes:
#'   - query_gff/subject_gff are accepted for consistency with other dndsR steps but are not used here.
#'   - For a simpler FASTA-only batch wrapper (no GFF columns), use seq_list_file (see below).
#'
#' @param seq_list_file Path to a whitespace-delimited (tabs/spaces) file OR a data.frame
#'   describing FASTA pairs to run with orthologr in batch.
#'   Must include columns: comparison_name, query_seq, subject_seq.
#'   These files are used directly as inputs to orthologr::dNdS() (CDS by default; set seq_type="protein" for proteins).
#'
#' @param comparison_name Name for the comparison (single mode). Used as the per-comparison
#'   folder name under \code{output_dir} and as the output TSV basename.
#' @param query_seq Path to the query sequence FASTA for \code{orthologr::dNdS()}
#'   (single mode). Required for single mode.
#' @param subject_seq Path to the subject sequence FASTA for \code{orthologr::dNdS()}
#'   (single mode). Required for single mode.
#'
#' @param query_fasta Path to the original query FASTA used only to derive the expected
#'   per-comparison sequence FASTA in batch mode.
#' @param subject_fasta Path to the original subject FASTA used only to derive the expected
#'   per-comparison sequence FASTA in batch mode.
#'
#' @param output_dir Root directory containing per-comparison folders. Each run uses
#'   \code{file.path(output_dir, comparison_name)}.
#' @param overwrite Logical; if \code{FALSE}, skip a comparison when the output TSV already exists.
#'   Default \code{FALSE}.
#'
#' @param seq_type Character; sequence type passed to \code{orthologr::dNdS()}.
#'   Must be \code{"cds"} or \code{"protein"}. When \code{"cds"}, input files must be CDS nucleotide
#'   FASTAs. When \code{"protein"}, input files must be protein FASTAs.
#' @param cds_suffix Filename suffix used in batch mode when \code{seq_type="cds"}.
#'   Default \code{"_CDS.fasta"}.
#' @param protein_suffix Filename suffix used in batch mode when \code{seq_type="protein"}.
#'   Default \code{"_AA.fasta"}.
#'
#' @param comp_cores Integer number of CPU cores to pass to \code{orthologr::dNdS()} as
#'   \code{comp_cores}. Default 4.
#' @param aligner Character. Aligner passed to \code{orthologr::dNdS()} (e.g. \code{"diamond"} or \code{"blast"}).
#'   If \code{"diamond"}, the \code{diamond} binary must be available on \code{PATH}.
#' @param sensitivity_mode Character. Passed to \code{orthologr::dNdS()} as \code{sensitivity_mode}.
#' @param dnds_method Character. Passed to \code{orthologr::dNdS()} as \code{dnds_est.method}.
#'
#' @param ... Additional arguments forwarded to \code{orthologr::dNdS()}. Only named arguments
#'   matching \code{formals(orthologr::dNdS)} are forwarded. If an argument overlaps a default
#'   set by this wrapper, it will override the wrapper default, except for \code{query_file},
#'   \code{subject_file}, and \code{seq_type} which are controlled by this wrapper.
#'
#' @details
#' \strong{Output path:} for each comparison, results are written to
#' \code{file.path(output_dir, comparison_name, paste0(comparison_name, "_dnds.tsv"))}.
#'
#' \strong{Batch mode file resolution:} when \code{query_fasta}/\code{subject_fasta} are used,
#' the input FASTAs are resolved as:
#' \itemize{
#'   \item \code{file.path(output_dir, comparison_name, paste0(basename(query_fasta), <suffix>))}
#'   \item \code{file.path(output_dir, comparison_name, paste0(basename(subject_fasta), <suffix>))}
#' }
#' where \code{<suffix>} is \code{cds_suffix} or \code{protein_suffix} depending on \code{seq_type}.
#'
#' This function requires the \pkg{orthologr} package. For Biostrings >= 2.77.1, alignment helpers
#' were moved into \pkg{pwalign}; this function applies a small shim to keep older orthologr code working.
#'
#' @return (Invisibly) a character vector of output TSV paths (batch mode), or a single output TSV path
#'   (single mode). Comparisons skipped due to missing inputs are omitted from batch outputs.
#' @export
calculate_dnds <- function(comparison_file = NULL,
                           seq_list_file = NULL,
                           comparison_name = NULL,
                           query_seq = NULL,
                           subject_seq = NULL,
                           query_fasta = NULL,
                           subject_fasta = NULL,
                           output_dir = getwd(),
                           overwrite = FALSE,
                           seq_type = c("cds", "protein"),
                           cds_suffix = "_CDS.fasta",
                           protein_suffix = "_AA.fasta",
                           comp_cores = 4,
                           aligner = "diamond",
                           sensitivity_mode = "fast",
                           dnds_method = "Comeron",
                           ...) {

  seq_type <- match.arg(seq_type)
  if (!is.null(comparison_file) && !is.null(seq_list_file)) {
    stop("Supply only one of comparison_file (batch mode) or seq_list_file (FASTA-only list mode).", call. = FALSE)
  }

  # ---- dependency guards ----
  .need_pkg <- function(pkg, msg = NULL) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (is.null(msg)) {
        msg <- sprintf(
          "Package '%s' is required for this step. Install it or run inside the dndsR container.",
          pkg
        )
      }
      stop(msg, call. = FALSE)
    }
  }

  .need_pkg(
    "orthologr",
    paste0(
      "This step needs 'orthologr'. Install via remotes::install_github('drostlab/orthologr') ",
      "or run inside the dndsR container (recommended)."
    )
  )
  .need_pkg(
    "pwalign",
    paste0(
      "Biostrings >= 2.77.1 moved pairwiseAlignment()/pattern()/subject()/writePairwiseAlignments() ",
      "into the 'pwalign' package."
    )
  )
  # ---- Biostrings >= 2.77.1 compatibility shim (CRAN-safe; no namespace edits) ----
  .orthologr_dNdS <- orthologr::dNdS

  if (requireNamespace("Biostrings", quietly = TRUE) &&
      utils::packageVersion("Biostrings") >= "2.77.1") {

    .need_pkg(
      "pwalign",
      "Biostrings >= 2.77.1 moved pairwiseAlignment helpers into the 'pwalign' package."
    )

    ortho_ns <- asNamespace("orthologr")
    shim_env <- new.env(parent = ortho_ns)

    # Helpers that moved out of Biostrings
    moved <- c("pairwiseAlignment", "pattern", "subject", "writePairwiseAlignments")
    for (nm in moved) {
      assign(nm, get(nm, envir = asNamespace("pwalign"), inherits = FALSE), envir = shim_env)
    }

    # Rewrite any explicit Biostrings:: calls to pwalign:: (and also catch bare calls)
    .patch_fun_text <- function(fun, from_pkg = "Biostrings", to_pkg = "pwalign") {
      if (!is.function(fun)) return(fun)

      txt <- paste(deparse(fun), collapse = "\n")

      # quick exit if nothing relevant
      if (!grepl(paste0(from_pkg, "::"), txt, fixed = TRUE) &&
          !grepl("pairwiseAlignment\\s*\\(", txt) &&
          !grepl("writePairwiseAlignments\\s*\\(", txt) &&
          !grepl("\\bpattern\\s*\\(", txt) &&
          !grepl("\\bsubject\\s*\\(", txt)) {
        return(fun)
      }

      # Replace explicit namespace-qualified calls first
      for (nm in moved) {
        txt <- gsub(paste0(from_pkg, "::", nm), paste0(to_pkg, "::", nm), txt, fixed = TRUE)
      }

      # Replace unqualified calls (if they exist) to force pwalign
      txt <- gsub("(?<![[:alnum:]_:.])pairwiseAlignment\\s*\\(", "pwalign::pairwiseAlignment(", txt, perl = TRUE)
      txt <- gsub("(?<![[:alnum:]_:.])writePairwiseAlignments\\s*\\(", "pwalign::writePairwiseAlignments(", txt, perl = TRUE)
      txt <- gsub("(?<![[:alnum:]_:.])pattern\\s*\\(", "pwalign::pattern(", txt, perl = TRUE)
      txt <- gsub("(?<![[:alnum:]_:.])subject\\s*\\(", "pwalign::subject(", txt, perl = TRUE)

      # Re-eval as a function and keep formals identical
      patched <- eval(parse(text = txt))
      environment(patched) <- shim_env
      patched
    }

    # Patch EVERY orthologr function that references Biostrings::pairwiseAlignment (or bare calls)
    ortho_objs <- ls(envir = ortho_ns, all.names = TRUE)
    for (nm in ortho_objs) {
      obj <- get(nm, envir = ortho_ns, inherits = FALSE)
      if (is.function(obj)) {
        patched <- .patch_fun_text(obj)
        # only install if actually changed / relevant (cheap heuristic)
        # install into shim_env so calls from patched dNdS resolve here
        assign(nm, patched, envir = shim_env)
      }
    }

    # Patch dNdS itself and run that
    f <- .patch_fun_text(orthologr::dNdS)
    environment(f) <- shim_env
    .orthologr_dNdS <- f

    if (isTRUE(getOption("dndsR.verbose_shims", FALSE))) {
      message("Using pwalign shim for orthologr under Biostrings >= 2.77.1 (auto-patched orthologr functions).")
    }
  }

  # ---- Optional: if using DIAMOND, confirm the binary is on PATH ----
  if (identical(tolower(aligner), "diamond")) {
    if (Sys.which("diamond") == "") {
      stop(
        "aligner='diamond' but the DIAMOND binary was not found on PATH.\n",
        "Install DIAMOND >= 2.1.8 or run inside the dndsR container.",
        call. = FALSE
      )
    }
  }

  # ---- helper: read comparison table (whitespace-delimited) ----
  .read_comparisons <- function(x) {

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
  
      while (length(first_line) == 0) {
        first_line <- readLines(con, n = 1)
        if (length(first_line) == 0) break
        first_line <- trimws(first_line)
        if (first_line == "") first_line <- ""
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

    # headerless 5-column schema inference: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff
    headerless_like <- !any(grepl("^(query_|subject_)", names(df)))
    if (!all(needed %in% names(df)) && headerless_like && ncol(df) >= 5) {
      names(df)[1:5] <- c("comparison_name", "query_fasta", "query_gff", "subject_fasta", "subject_gff")
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
    
    df
  }

  # ---- helper: read FASTA-only list (comparison_name, query_seq, subject_seq) ----
  .read_seq_list <- function(x) {

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

    if (is.data.frame(x)) {
      df <- x
    } else {

      first_line <- ""
      con <- file(x, open = "r")
      on.exit(close(con), add = TRUE)

      while (length(first_line) == 0) {
        first_line <- readLines(con, n = 1)
        if (length(first_line) == 0) break
        first_line <- trimws(first_line)
        if (first_line == "") first_line <- ""
      }

      has_header <- FALSE
      if (nzchar(first_line)) {
        toks <- strsplit(first_line, "\\s+")[[1]]
        has_header <- any(toks %in% c("comparison_name", "query_seq", "subject_seq"))
      }

      df <- read_ws(has_header)
    }

    if (!("comparison_name" %in% names(df)) && ncol(df) >= 1) names(df)[1] <- "comparison_name"
    if (!("query_seq" %in% names(df)) && ncol(df) >= 2) names(df)[2] <- "query_seq"
    if (!("subject_seq" %in% names(df)) && ncol(df) >= 3) names(df)[3] <- "subject_seq"

    needed <- c("comparison_name", "query_seq", "subject_seq")
    if (!all(needed %in% names(df))) {
      stop(
        "seq_list_file must include columns: comparison_name, query_seq, subject_seq (headerless 3-column files are OK).",
        call. = FALSE
      )
    }

    df
  }

  # ---- helper: resolve input FASTA paths for a comparison ----
  .resolve_inputs <- function(comp, comp_dir, row) {
    # single mode wins if provided (and non-empty)
    if ("query_seq" %in% names(row) && "subject_seq" %in% names(row)) {
      qs <- row[["query_seq"]]
      ss <- row[["subject_seq"]]
      if (!is.na(qs) && nzchar(qs) && !is.na(ss) && nzchar(ss)) {
        return(list(
          query_file = normalizePath(qs),
          subject_file = normalizePath(ss),
          mode = "single"
        ))
      }
    }

    # batch mode: derive from query_fasta/subject_fasta basename + suffix
    if (!("query_fasta" %in% names(row)) || !("subject_fasta" %in% names(row))) {
      return(list(query_file = NA_character_, subject_file = NA_character_, mode = "none"))
    }

    qfa <- row[["query_fasta"]]
    sfa <- row[["subject_fasta"]]
    if (is.na(qfa) || !nzchar(qfa) || is.na(sfa) || !nzchar(sfa)) {
      return(list(query_file = NA_character_, subject_file = NA_character_, mode = "none"))
    }

    q_base <- tools::file_path_sans_ext(basename(qfa))
    s_base <- tools::file_path_sans_ext(basename(sfa))

    suffix <- if (seq_type == "cds") cds_suffix else protein_suffix

    q_in <- file.path(comp_dir, paste0(q_base, suffix))
    s_in <- file.path(comp_dir, paste0(s_base, suffix))

    list(
      query_file = normalizePath(q_in, mustWork = FALSE),
      subject_file = normalizePath(s_in, mustWork = FALSE),
      mode = "batch"
    )
  }

  # ---- helper: one comparison ----
  run_dnds <- function(comparison_basename, row) {
    comp_dir <- file.path(output_dir, comparison_basename)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

    out_tsv <- file.path(comp_dir, paste0(comparison_basename, "_dnds.tsv"))

    if (!isTRUE(overwrite) && file.exists(out_tsv) && file.info(out_tsv)$size > 0) {
      message("Skipping ", comparison_basename, " (exists).")
      return(out_tsv)
    }

    inp <- .resolve_inputs(comparison_basename, comp_dir, row)
    q_in <- inp$query_file
    s_in <- inp$subject_file

    if (!nzchar(q_in) || !nzchar(s_in) || !file.exists(q_in) || !file.exists(s_in)) {
      warning(
        "Sequence files not found for ", comparison_basename, " (mode=", inp$mode, ").\n",
        "* query_file:   ", q_in, "\n",
        "* subject_file: ", s_in, "\n",
        "Skipping."
      )
      return(NULL)
    }

    # defaults you provide
    args <- list(
      query_file         = normalizePath(q_in),
      subject_file       = normalizePath(s_in),
      aligner            = aligner,
      sensitivity_mode   = sensitivity_mode,
      seq_type           = seq_type,
      format             = "fasta",
      ortho_detection    = "RBH",
      delete_corrupt_cds = TRUE,
      eval               = "1E-5",
      aa_aln_type        = "pairwise",
      aa_aln_tool        = "NW",
      codon_aln_tool     = "pal2nal",
      dnds_est.method    = dnds_method,
      comp_cores         = as.integer(comp_cores),
      quiet              = TRUE,
      clean_folders      = FALSE,
      print_citation     = FALSE
    )

    # user-supplied overrides / additions
    user <- list(...)
    if (length(user)) {
      nms <- names(user)

      # Drop unnamed/invalid extras
      if (is.null(nms)) {
        user <- list()
      } else {
        keep <- nzchar(nms) & !is.na(nms)
        if (!all(keep)) {
          message("Ignoring ", sum(!keep), " unnamed/invalid extra argument(s) passed via '...'.")
          user <- user[keep]
          nms  <- names(user)
        }
      }

      if (length(user)) {
        # Prevent desync: these are controlled by wrapper
        if ("seq_type" %in% names(user)) {
          message("Ignoring ...$seq_type (use calculate_dnds(seq_type=...) instead).")
          user$seq_type <- NULL
        }
        if ("query_file" %in% names(user)) {
          message("Ignoring ...$query_file (use query_seq/subject_seq in single mode instead).")
          user$query_file <- NULL
        }
        if ("subject_file" %in% names(user)) {
          message("Ignoring ...$subject_file (use query_seq/subject_seq in single mode instead).")
          user$subject_file <- NULL
        }

        # Only forward arguments orthologr::dNdS knows
        dnds_formals <- names(formals(orthologr::dNdS))
        nms <- names(user)
        user <- user[nms %in% dnds_formals]

        overlap <- intersect(names(user), names(args))
        if (length(overlap)) {
          message("Overriding defaults: ", paste(overlap, collapse = ", "))
          args[overlap] <- user[overlap]
        }

        extra <- setdiff(names(user), names(args))
        if (length(extra)) {
          args <- c(args, user[extra])
        }
      }
    }

    message(
      "[calculate_dnds] Running dN/dS for: ", comparison_basename,
      " (seq_type=", seq_type, ", mode=", inp$mode, ")"
    )

    res <- do.call(.orthologr_dNdS, args)

    utils::write.table(
      res,
      file = out_tsv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    message("[calculate_dnds] Finished: ", comparison_basename)
    out_tsv
  }
  
  # ---- batch vs single ----
  out_paths <- character(0)

  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)

    for (i in seq_len(nrow(df))) {
      comp <- df$comparison_name[i]
      row  <- as.list(df[i, , drop = FALSE])

      out <- run_dnds(comparison_basename = comp, row = row)
      if (!is.null(out)) out_paths <- c(out_paths, out)
    }

    message("All dN/dS calculations complete.")
    return(invisible(out_paths))
  }

  if (!is.null(seq_list_file)) {
    df <- .read_seq_list(seq_list_file)

    for (i in seq_len(nrow(df))) {
      comp <- df$comparison_name[i]
      row  <- as.list(df[i, , drop = FALSE])

      out <- run_dnds(comparison_basename = comp, row = row)
      if (!is.null(out)) out_paths <- c(out_paths, out)
    }

    message("All dN/dS calculations complete.")
    return(invisible(out_paths))
  }

  # ---- single mode ----
  if (is.null(comparison_name) || !nzchar(comparison_name)) {
    stop("In single mode, supply comparison_name.", call. = FALSE)
  }

  # single mode if query_seq+subject_seq supplied; otherwise batch mode needs query_fasta+subject_fasta
  if (!is.null(query_seq) && !is.null(subject_seq) && nzchar(query_seq) && nzchar(subject_seq)) {
    row <- list(query_seq = query_seq, subject_seq = subject_seq)
  } else {
    if (is.null(query_fasta) || is.null(subject_fasta) || !nzchar(query_fasta) || !nzchar(subject_fasta)) {
      stop(
        paste0(
          "In single mode, supply either:\n",
          "* query_seq + subject_seq (single mode), OR\n",
          "* query_fasta + subject_fasta (batch mode)."
        ),
        call. = FALSE
      )
    }
    row <- list(query_fasta = query_fasta, subject_fasta = subject_fasta)
  }

  out <- run_dnds(comparison_basename = comparison_name, row = row)
  message("All dN/dS calculations complete.")
  invisible(out)
}
