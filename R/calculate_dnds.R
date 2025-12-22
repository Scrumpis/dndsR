#' Calculate dN/dS from pre-extracted CDS FASTAs using orthologr
#'
#' Runs \code{orthologr::dNdS()} on per-comparison CDS FASTA files that are expected to
#' already exist on disk. Supports batch mode via a comparison table or single-run mode
#' via explicit arguments.
#'
#' @param comparison_file Path to a whitespace-delimited (tabs/spaces) file OR a data.frame
#'   with columns: \code{comparison_name}, \code{query_fasta}, \code{query_gff},
#'   \code{subject_fasta}, \code{subject_gff}. If provided, batch mode is used.
#' @param comparison_name Name for the comparison (single mode). Used as the per-comparison
#'   folder name under \code{output_dir} and as the output TSV basename.
#' @param query_fasta Path to the original query FASTA (single mode). Used only to derive
#'   the expected CDS filename \code{<basename>_CDS.fasta} under the comparison folder.
#' @param subject_fasta Path to the original subject FASTA (single mode). Used only to derive
#'   the expected CDS filename \code{<basename>_CDS.fasta} under the comparison folder.
#' @param query_gff Path to query GFF3 (single mode). Currently not used by this function,
#'   but accepted for consistency with comparison tables.
#' @param subject_gff Path to subject GFF3 (single mode). Currently not used by this function,
#'   but accepted for consistency with comparison tables.
#' @param output_dir Root directory containing per-comparison folders. Each run uses
#'   \code{file.path(output_dir, comparison_name)}.
#' @param comp_cores Integer number of CPU cores to pass to \code{orthologr::dNdS()}
#'   as \code{comp_cores}. Default 4.
#' @param aligner Character. Aligner passed to \code{orthologr::dNdS()} (e.g. \code{"diamond"}).
#'   If \code{"diamond"}, the \code{diamond} binary must be available on \code{PATH}.
#' @param sensitivity_mode Character. Passed to \code{orthologr::dNdS()} as \code{sensitivity_mode}
#'   (e.g. \code{"fast"}).
#' @param dnds_method Character. Passed to \code{orthologr::dNdS()} as \code{dnds_est.method}
#'   (e.g. \code{"Comeron"}).
#' @param ... Additional arguments forwarded to \code{orthologr::dNdS()}. Only named arguments
#'   that match \code{formals(orthologr::dNdS)} are forwarded. If an argument overlaps a default
#'   set by this wrapper, it will override the wrapper default.
#'
#' @details
#' Expected input CDS FASTAs (per comparison) are:
#' \itemize{
#'   \item \code{file.path(output_dir, comparison_name, paste0(basename(query_fasta_noext), "_CDS.fasta"))}
#'   \item \code{file.path(output_dir, comparison_name, paste0(basename(subject_fasta_noext), "_CDS.fasta"))}
#' }
#' If either CDS file is missing, the comparison is skipped (returns \code{NULL} for that run).
#' If the output TSV already exists, the run is skipped and the existing path is returned.
#'
#' This function requires the \pkg{orthologr} package. For Biostrings >= 2.77.1,
#' alignment helpers were moved into \pkg{pwalign}; this function applies a small shim
#' to keep older \pkg{orthologr} code working.
#'
#' @return (Invisibly) a character vector of output TSV paths (batch mode), or a single
#'   output TSV path (single mode). Comparisons that are skipped due to missing CDS files
#'   are omitted.
#'
#' @export
calculate_dnds <- function(comparison_file = NULL,
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
    "This step needs 'orthologr'. Install via remotes::install_github('drostlab/orthologr') "
    |> paste0("or run inside the dndsR container (recommended).")
  )
  .need_pkg(
    "pwalign",
    "Biostrings >= 2.77.1 moved pairwiseAlignment()/pattern()/subject()/writePairwiseAlignments() into the 'pwalign' package."
  )

  # Shim: if Biostrings made core alignment helpers defunct, re-point them to pwalign's impl
  .fix_pwalign_shim <- function() {
    if (!requireNamespace("Biostrings", quietly = TRUE) ||
        !requireNamespace("pwalign",    quietly = TRUE)) {
      return(invisible())
    }

    bs_ver <- utils::packageVersion("Biostrings")
    if (bs_ver < "2.77.1") return(invisible())

    ns <- asNamespace("Biostrings")

    rebind_if_exists <- function(name, fun) {
      if (!is.function(fun)) return(invisible())
      # Only touch it if Biostrings actually has this symbol
      if (!exists(name, envir = ns, inherits = FALSE)) return(invisible())

      locked <- tryCatch(bindingIsLocked(name, ns),
                         error = function(e) FALSE)
      if (isTRUE(locked)) {
        unlockBinding(name, ns)
        # (Optional) you could relock here with on.exit, but not required.
      }

      assign(name, fun, envir = ns)
      invisible()
    }

    # Core functions that Biostrings says moved to pwalign
    rebind_if_exists("pairwiseAlignment",       pwalign::pairwiseAlignment)
    rebind_if_exists("pattern",                 pwalign::pattern)
    rebind_if_exists("subject",                 pwalign::subject)
    rebind_if_exists("writePairwiseAlignments", pwalign::writePairwiseAlignments)

    invisible()
  }

  .fix_pwalign_shim()


  # Optional but helpful: if using DIAMOND, confirm the binary is on PATH
  if (identical(tolower(aligner), "diamond")) {
    ok <- !inherits(suppressWarnings(try(system2("diamond", "--version"), silent = TRUE)), "try-error")
    if (inherits(ok, "try-error")) {
      stop(
        paste0(
          "aligner='diamond' but the DIAMOND binary was not found on PATH.\n",
          "* If you are using the container, this should already be installed and on PATH.\n",
          "* Otherwise, install DIAMOND >= 2.1.8 and ensure 'diamond' is callable from your shell."
        ),
        call. = FALSE
      )
    }
  }

  # ---- helper: read comparison_file with ANY whitespace as delimiter ----
  .read_comparisons <- function(x) {
    req <- c("comparison_name","query_fasta","query_gff","subject_fasta","subject_gff")
    if (is.data.frame(x)) {
      stopifnot(all(req %in% names(x)))
      return(x[, req, drop = FALSE])
    }
    read_ws <- function(hdr) {
      utils::read.table(x,
                        header = hdr,
                        sep = "",               # any whitespace (tabs OR spaces)
                        quote = "\"",           # allow quoting paths if they contain spaces
                        stringsAsFactors = FALSE,
                        comment.char = "",
                        strip.white = TRUE,
                        blank.lines.skip = TRUE,
                        check.names = FALSE)
    }
    # Try headered first
    df1 <- try(read_ws(TRUE), silent = TRUE)
    if (!inherits(df1, "try-error") && all(req %in% names(df1))) {
      return(df1[, req, drop = FALSE])
    }
    # Fallback: headerless
    df2 <- read_ws(FALSE)
    if (ncol(df2) < 5) {
      stop("comparison_file must have 5 columns or a header with: ",
           paste(req, collapse = ", "), call. = FALSE)
    }
    names(df2)[1:5] <- req
    df2[, req, drop = FALSE]
  }

  # ---- helper: one comparison ----
  run_dnds <- function(comparison_basename, query_fa, subject_fa, query_gff, subject_gff, ...) {
    comp_dir <- file.path(output_dir, comparison_basename)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

    q_base <- tools::file_path_sans_ext(basename(query_fa))
    s_base <- tools::file_path_sans_ext(basename(subject_fa))
    query_CDS   <- file.path(comp_dir, paste0(q_base, "_CDS.fasta"))
    subject_CDS <- file.path(comp_dir, paste0(s_base, "_CDS.fasta"))
    out_tsv     <- file.path(comp_dir, paste0(comparison_basename, "_dnds.tsv"))

    if (!file.exists(query_CDS) || !file.exists(subject_CDS)) {
      warning("CDS files not found in ", comp_dir, ". Expected: ",
              basename(query_CDS), " and ", basename(subject_CDS),
              ". Skipping ", comparison_basename, ".")
      return(NULL)
    }
    if (file.exists(out_tsv)) {
      message("Skipping ", comparison_basename, " (exists).")
      return(out_tsv)
    }

    # defaults you provide
    args <- list(
      query_file         = normalizePath(query_CDS),
      subject_file       = normalizePath(subject_CDS),
      aligner            = aligner,
      sensitivity_mode   = sensitivity_mode,
      seq_type           = "cds",
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
        # user-supplied overrides / additions
    user <- list(...)
    if (length(user)) {
      # Drop unnamed / invalid extras (e.g. stray alist() from launcher/optparse glue)
      nms <- names(user)
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
        # Only forward arguments that orthologr::dNdS actually knows about
        dnds_formals <- names(formals(orthologr::dNdS))
        user <- user[nms %in% dnds_formals]

        overlap <- intersect(names(user), names(args))
        if (length(overlap)) {
          message("Overriding defaults: ", paste(overlap, collapse = ", "))
          args[overlap] <- user[overlap]  # override defaults
        }

        extra <- setdiff(names(user), names(args))
        if (length(extra)) {
          args <- c(args, user[extra])    # add additional valid dNdS args
        }
      }
    }
    message("Running dN/dS for: ", comparison_basename)
    res <- do.call(orthologr::dNdS, args)
    utils::write.table(res, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Finished: ", comparison_basename)
    out_tsv
  }

  # ---- batch vs single ----
  out_paths <- character(0)

  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)
    for (i in seq_len(nrow(df))) {
      out <- run_dnds(
        comparison_basename = df$comparison_name[i],
        query_fa   = normalizePath(df$query_fasta[i]),
        subject_fa = normalizePath(df$subject_fasta[i]),
        query_gff  = normalizePath(df$query_gff[i]),
        subject_gff= normalizePath(df$subject_gff[i])
      )
      if (!is.null(out)) out_paths <- c(out_paths, out)
    }
    message("All dN/dS calculations complete.")
    return(invisible(out_paths))
  }

  # Single mode
  if (any(vapply(list(comparison_name, subject_fasta, query_fasta, subject_gff, query_gff),
                 is.null, logical(1)))) {
    stop("In single mode, supply comparison_name, subject_fasta, query_fasta, subject_gff, and query_gff.", call. = FALSE)
  }
  out <- run_dnds(
    comparison_basename = comparison_name,
    query_fa   = normalizePath(query_fasta),
    subject_fa = normalizePath(subject_fasta),
    query_gff  = normalizePath(query_gff),
    subject_gff= normalizePath(subject_gff)
  )
  message("All dN/dS calculations complete.")
  invisible(out)
}
