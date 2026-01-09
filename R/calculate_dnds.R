#' Calculate dN/dS in batch or single mode using orthologr
#'
#' Wraps \code{orthologr::dNdS()} to provide batch execution, sensible defaults,
#' reproducible container-friendly behavior, and standardized output locations.
#'
#' This function supports two input styles:
#' \itemize{
#'   \item \strong{Direct mode}: provide \code{query_seq} and \code{subject_seq} paths
#'     (CDS or protein FASTAs depending on \code{seq_type}).
#'   \item \strong{Pipeline mode}: provide \code{query_fasta} and \code{subject_fasta} only
#'     to derive expected per-comparison sequence FASTAs within \code{output_dir} using
#'     standardized suffixes (e.g. \code{_CDS.fasta} or \code{_AA.fasta}).
#' }
#'
#' @param comparison_file Path to a whitespace-delimited (tabs/spaces) file OR a data.frame
#'   describing comparisons (batch mode).
#'   Must include \code{comparison_name} and either:
#'   \itemize{
#'     \item \code{query_seq} and \code{subject_seq} (direct mode), OR
#'     \item \code{query_fasta} and \code{subject_fasta} (pipeline mode).
#'   }
#'   Additional columns are allowed and ignored. For backward compatibility, files with
#'   \code{query_gff}/\code{subject_gff} columns are accepted but those columns are not used.
#'
#' @param comparison_name Name for the comparison (single mode). Used as the per-comparison
#'   folder name under \code{output_dir} and as the output TSV basename.
#' @param query_seq Path to the query sequence FASTA used directly by \code{orthologr::dNdS()}
#'   (single mode). Required for direct mode.
#' @param subject_seq Path to the subject sequence FASTA used directly by \code{orthologr::dNdS()}
#'   (single mode). Required for direct mode.
#'
#' @param query_fasta Path to the original query FASTA used only to derive the expected
#'   per-comparison sequence FASTA in pipeline mode (single mode).
#' @param subject_fasta Path to the original subject FASTA used only to derive the expected
#'   per-comparison sequence FASTA in pipeline mode (single mode).
#'
#' @param output_dir Root directory containing per-comparison folders. Each run uses
#'   \code{file.path(output_dir, comparison_name)}.
#' @param overwrite Logical; if \code{FALSE}, skip a comparison when the output TSV already exists.
#'   Default \code{FALSE}.
#'
#' @param seq_type Character; sequence type passed to \code{orthologr::dNdS()}.
#'   Must be \code{"cds"} or \code{"protein"}. When \code{"cds"}, input files must be CDS nucleotide
#'   FASTAs. When \code{"protein"}, input files must be protein FASTAs.
#' @param cds_suffix Filename suffix used in pipeline mode when \code{seq_type="cds"}.
#'   Default \code{"_CDS.fasta"}.
#' @param protein_suffix Filename suffix used in pipeline mode when \code{seq_type="protein"}.
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
#' \strong{Pipeline mode file resolution:} when \code{query_fasta}/\code{subject_fasta} are used,
#' the input FASTAs are resolved as:
#' \itemize{
#'   \item \code{file.path(output_dir, comparison_name, paste0(basename(query_fasta_noext), <suffix>))}
#'   \item \code{file.path(output_dir, comparison_name, paste0(basename(subject_fasta_noext), <suffix>))}
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
  # ---- Biostrings >= 2.77.1 compatibility shim (pairwiseAlignment moved to pwalign) ----
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    bs_ver <- utils::packageVersion("Biostrings")

    if (bs_ver >= "2.77.1") {
      # Biostrings >= 2.77.1 moved these helpers into pwalign.
      # Some orthologr versions still call them unqualified (via imports),
      # so we rebind them inside orthologr's namespace to pwalign.
      shim_funs <- c(
        "pairwiseAlignment",
        "pattern",
        "subject",
        "writePairwiseAlignments"
      )

      ortho_ns   <- asNamespace("orthologr")
      pwalign_ns <- asNamespace("pwalign")

      for (fn in shim_funs) {
        # Only patch if:
        # 1) pwalign provides it, AND
        # 2) orthologr already has a binding for it (we can't add new symbols to a locked namespace)
        if (exists(fn, envir = pwalign_ns, inherits = FALSE) &&
            exists(fn, envir = ortho_ns,   inherits = FALSE)) {

          was_locked <- bindingIsLocked(fn, ortho_ns)
          if (was_locked) unlockBinding(fn, ortho_ns)

          assign(fn, get(fn, envir = pwalign_ns, inherits = FALSE), envir = ortho_ns)

          if (was_locked) lockBinding(fn, ortho_ns)
        }
      }

      if (isTRUE(getOption("dndsR.verbose_shims", FALSE))) message(
        "Biostrings ", bs_ver,
        " detected; using pwalign for pairwiseAlignment helpers."
      )
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
          "query_seq", "subject_seq",
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
  
    # ---- determine mode ----
    has_direct   <- all(c("query_seq", "subject_seq") %in% names(df))
    has_pipeline <- all(c("query_fasta", "subject_fasta") %in% names(df))
  
    # ---- headerless schema inference ----
    if (!has_direct && !has_pipeline) {
  
      headerless_like <- !any(grepl("^(query_|subject_)", names(df)))
  
      if (!is.data.frame(x) && headerless_like) {
  
        if (ncol(df) >= 5) {
          names(df)[1:5] <- c(
            "comparison_name",
            "query_fasta", "query_gff",
            "subject_fasta", "subject_gff"
          )
          has_pipeline <- TRUE
  
        } else if (ncol(df) >= 3) {
          names(df)[1:3] <- c(
            "comparison_name",
            "query_seq", "subject_seq"
          )
          has_direct <- TRUE
  
        } else {
          stop(
            paste0(
              "comparison_file must include either:\n",
              "* 3 columns: comparison_name, query_seq, subject_seq (direct mode), OR\n",
              "* 5 columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff (pipeline mode).\n",
              "Additional columns are allowed."
            ),
            call. = FALSE
          )
        }
  
      } else {
        stop(
          paste0(
            "comparison_file must include either:\n",
            "* query_seq and subject_seq (direct mode), OR\n",
            "* query_fasta and subject_fasta (pipeline mode).\n",
            "Additional columns are allowed."
          ),
          call. = FALSE
        )
      }
    }
  
    df
  }

  # ---- helper: resolve input FASTA paths for a comparison ----
  .resolve_inputs <- function(comp, comp_dir, row) {
    # Direct mode wins if provided (and non-empty)
    if ("query_seq" %in% names(row) && "subject_seq" %in% names(row)) {
      qs <- row[["query_seq"]]
      ss <- row[["subject_seq"]]
      if (!is.na(qs) && nzchar(qs) && !is.na(ss) && nzchar(ss)) {
        return(list(
          query_file = normalizePath(qs),
          subject_file = normalizePath(ss),
          mode = "direct"
        ))
      }
    }

    # Pipeline mode: derive from query_fasta/subject_fasta basename + suffix
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
      mode = "pipeline"
    )
  }

  # ---- helper: one comparison ----
  run_dnds <- function(comparison_basename, row) {
    comp_dir <- file.path(output_dir, comparison_basename)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

    # Per-comparison log file (matches go_enrichment style)
    log_file <- file.path(comp_dir, sprintf("%s_dnds.log", comparison_basename))

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
          message("Ignoring ...$query_file (use query_seq/subject_seq in direct mode instead).")
          user$query_file <- NULL
        }
        if ("subject_file" %in% names(user)) {
          message("Ignoring ...$subject_file (use query_seq/subject_seq in direct mode instead).")
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

    message(sprintf(
      "[calculate_dnds] %s (logging to %s)",
      comparison_basename, log_file
    ))

    .dndsr_with_log(
      log_file = log_file,
      tag = "calculate_dnds",
      header = c(
        sprintf("[calculate_dnds] comp=%s", comparison_basename),
        sprintf("[calculate_dnds] pid=%d", Sys.getpid()),
        sprintf("[calculate_dnds] seq_type=%s mode=%s", seq_type, inp$mode),
        sprintf(
          "[calculate_dnds] comp_cores=%d aligner=%s sensitivity_mode=%s dnds_method=%s",
          as.integer(comp_cores), aligner, sensitivity_mode, dnds_method
        ),
        sprintf("[calculate_dnds] query_file=%s", normalizePath(q_in)),
        sprintf("[calculate_dnds] subject_file=%s", normalizePath(s_in)),
        sprintf("[calculate_dnds] out_tsv=%s", out_tsv)
      ),
      expr = {
        message(
          "Running dN/dS for: ", comparison_basename,
          " (seq_type=", seq_type, ", mode=", inp$mode, ")"
        )

        res <- do.call(orthologr::dNdS, args)

        utils::write.table(
          res,
          file = out_tsv,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )

        message("Finished: ", comparison_basename)
        out_tsv
      }
    )
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

  # ---- single mode ----
  if (is.null(comparison_name) || !nzchar(comparison_name)) {
    stop("In single mode, supply comparison_name.", call. = FALSE)
  }

  # Direct mode if query_seq+subject_seq supplied; otherwise pipeline mode needs query_fasta+subject_fasta
  if (!is.null(query_seq) && !is.null(subject_seq) && nzchar(query_seq) && nzchar(subject_seq)) {
    row <- list(query_seq = query_seq, subject_seq = subject_seq)
  } else {
    if (is.null(query_fasta) || is.null(subject_fasta) || !nzchar(query_fasta) || !nzchar(subject_fasta)) {
      stop(
        paste0(
          "In single mode, supply either:\n",
          "* query_seq + subject_seq (direct mode), OR\n",
          "* query_fasta + subject_fasta (pipeline mode)."
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
