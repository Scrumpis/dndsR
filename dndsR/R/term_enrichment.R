#' Term enrichment (Fisher test) for query/subject terms in dNdS results
#'
#' Reads <comp>/<comp>_dnds_annot.tsv and tests enrichment of terms (q_* / s_* columns)
#' among positively selected pairs (dNdS > pos_threshold) vs the filtered background.
#'
#' @param dnds_annot_file Path to a single <comp>_dnds_annot.tsv (single mode).
#' @param comparison_file Path to whitespace-delimited file (tabs/spaces; header or not)
#'   with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#'   If provided, batch mode reads:
#'   file.path(output_dir, comparison_name, paste0(comparison_name, "_dnds_annot.tsv")).
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#' @param terms Character vector of term types to test (e.g., c("ipr","go","kegg")).
#'   Default NULL = auto-detect among {ipr,go,kegg,panther,pfam,tigrfam,cog}.
#'   Pass custom names explicitly via `terms=` (columns must be present as q_* / s_*).
#' @param sides Character vector among c("query","subject"). Default both.
#' @param pos_threshold Numeric. dNdS > pos_threshold defines "positive" (default 1).
#' @param max_dnds Numeric. Drop rows with dNdS >= max_dnds (default 10) or NA dNdS.
#' @param filter_expr Optional character with a logical expression evaluated in the data
#'   (e.g., "q_seqname == s_seqname").
#' @param make_plots Logical; if TRUE, write a top-N bubble plot per result (default TRUE).
#' @param top_n Integer; number of rows for plot (default 20).
#' @param drop_rows_without_term Logical; if TRUE (default), rows with no term of the tested
#'   type are removed from BOTH positives and background (annotation-aware universe).
#' @param min_total Minimum total occurrences (pos+nonpos) required for a term (default 0).
#' @param min_pos   Minimum positive occurrences required for a term (default 0).
#' @param fdr_method One of "BH","IHW","qvalue","none". Defaults to "BH".
#'                   If "IHW"/"qvalue" is chosen but the package is missing, falls back to "BH".
#' @param alpha FDR level for IHW weighting (default 0.05).
#' @return In single mode: (invisibly) list of data.frames by side_term.
#'         In batch mode: (invisibly) vector of output TSV paths.
#' @export
term_enrichment <- function(dnds_annot_file = NULL,
                            comparison_file = NULL,
                            output_dir = getwd(),
                            terms = NULL,
                            sides = c("query","subject"),
                            pos_threshold = 1,
                            max_dnds = 10,
                            filter_expr = NULL,
                            make_plots = TRUE,
                            top_n = 20,
                            drop_rows_without_term = TRUE,
                            min_total = 0,
                            min_pos = 0,
                            fdr_method = c("BH","IHW","qvalue","none"),
                            alpha = 0.05,
                            # deprecated alias for backward compatibility:
                            dnds_file = NULL) {

  # ---- deprecation shim ----
  if (!is.null(dnds_file) && is.null(dnds_annot_file)) {
    warning("`dnds_file` is deprecated; use `dnds_annot_file` instead.", call. = FALSE)
    dnds_annot_file <- dnds_file
  }

  sides <- match.arg(sides, choices = c("query","subject"), several.ok = TRUE)
  fdr_method <- match.arg(fdr_method)

  .read_ws <- function(path, header_try = TRUE) {
    utils::read.table(path, header = header_try, sep = "", quote = "\"",
                      stringsAsFactors = FALSE, comment.char = "",
                      strip.white = TRUE, blank.lines.skip = TRUE, check.names = FALSE)
  }

  .read_comparisons <- function(x) {
    req <- c("comparison_name","query_fasta","query_gff","subject_fasta","subject_gff")
    if (is.data.frame(x)) { stopifnot(all(req %in% names(x))); return(x[, req, drop = FALSE]) }
    df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
    if (!inherits(df1, "try-error") && all(req %in% names(df1))) return(df1[, req, drop = FALSE])
    df2 <- .read_ws(x, header_try = FALSE)
    if (ncol(df2) < 5) stop("comparison_file must have 5 cols or a header with: ", paste(req, collapse = ", "))
    names(df2)[1:5] <- req
    df2[, req, drop = FALSE]
  }

  .split_terms_unique <- function(x) {
    if (is.null(x)) return(character(0))
    s <- as.character(x)[1]
    if (is.na(s) || !nzchar(s)) return(character(0))
    unique(strsplit(s, ";", fixed = TRUE)[[1]])
  }

  # Auto-detect only among known term-types; skip non-character columns
  .available_terms <- function(df) {
    allowed <- c("ipr","go","kegg","panther","pfam","tigrfam","cog")
    have <- function(suf) {
      qn <- paste0("q_", suf); sn <- paste0("s_", suf)
      (qn %in% names(df) && is.character(df[[qn]])) ||
        (sn %in% names(df) && is.character(df[[sn]]))
    }
    tolower(allowed[vapply(allowed, have, logical(1))])
  }

  .apply_filter <- function(d, filter_expr) {
    keep <- !is.na(d$dNdS) & d$dNdS < max_dnds
    if (!is.null(filter_expr) && nzchar(filter_expr)) {
      ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)
      if (!inherits(ok, "try-error")) keep <- keep & isTRUE(as.vector(ok))
    }
    d[keep, , drop = FALSE]
  }

  # Build one enrichment table for a term-column (vector of semicolon-strings)
  .fisher_one_term <- function(vec_all, vec_pos, drop_empty = TRUE) {
    vec_all <- as.character(vec_all)
    vec_pos <- as.character(vec_pos)

    if (drop_empty) {
      vec_all <- vec_all[!is.na(vec_all) & nzchar(vec_all)]
      vec_pos <- vec_pos[!is.na(vec_pos) & nzchar(vec_pos)]
    }

    n_pos <- length(vec_pos)
    n_all <- length(vec_all)
    n_bg  <- n_all - n_pos
    if (n_pos == 0L || n_all == 0L || n_bg < 0L) return(NULL)

    count_terms <- function(v) {
      if (!length(v)) return(integer(0))
      tabs <- table(unlist(lapply(v, .split_terms_unique), use.names = FALSE))
      tabs[names(tabs) != ""]
    }

    all_tab <- count_terms(vec_all)
    pos_tab <- count_terms(vec_pos)
    if (length(all_tab) == 0) return(NULL)

    res <- do.call(rbind, lapply(names(all_tab), function(tt) {
      a <- as.integer(ifelse(tt %in% names(pos_tab), pos_tab[[tt]], 0L))     # pos with term
      c <- as.integer(all_tab[[tt]] - a)                                     # non-pos with term
      b <- n_pos - a                                                         # pos without term
      d <- n_bg  - c                                                         # non-pos without term
      ft <- try(stats::fisher.test(matrix(c(a,b,c,d), nrow = 2),
                                   alternative = "greater"), silent = TRUE)
      data.frame(
        term         = tt,
        pos_count    = a,
        nonpos_count = c,
        pos_total    = n_pos,
        nonpos_total = n_bg,
        odds_ratio   = if (inherits(ft, "try-error")) NA_real_ else as.numeric(ft$estimate),
        p_value      = if (inherits(ft, "try-error")) NA_real_ else ft$p.value,
        stringsAsFactors = FALSE
      )
    }))
    if (is.null(res) || !nrow(res)) return(NULL)
    res$total_count <- res$pos_count + res$nonpos_count
    res
  }

  .write_plot <- function(df, ylab, out_svg) {
    if (!make_plots || !requireNamespace("ggplot2", quietly = TRUE) || !nrow(df)) return(invisible(NULL))
    top <- df[order(df$p_adj, -df$pos_count, df$term), , drop = FALSE]
    top <- utils::head(top, top_n)
    gg <- ggplot2::ggplot(top, ggplot2::aes(x = pos_count, y = stats::reorder(term, -p_adj),
                                            size = pos_count / (pos_total + nonpos_total),
                                            color = p_adj)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradient(low = "red", high = "blue") +
      ggplot2::labs(x = "Positive count", y = ylab, size = "Pos / All", color = "adj p") +
      ggplot2::theme_minimal(base_size = 13)
    ggplot2::ggsave(out_svg, gg, width = 11, height = 9)
    invisible(NULL)
  }

  .adjust_pvals <- function(res, method, alpha) {
    if (!nrow(res)) { res$p_adj <- numeric(0); return(res) }
    if (method == "BH") {
      res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
    } else if (method == "qvalue") {
      if (requireNamespace("qvalue", quietly = TRUE)) {
        qobj <- qvalue::qvalue(res$p_value)
        res$p_adj <- as.numeric(qobj$qvalues)
      } else {
        warning("qvalue not installed; falling back to BH.")
        res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
      }
    } else if (method == "IHW") {
      if (requireNamespace("IHW", quietly = TRUE)) {
        ihw_obj <- IHW::ihw(res$p_value ~ res$total_count, alpha = alpha)
        res$p_adj <- as.numeric(IHW::adj_pvalues(ihw_obj))
      } else {
        warning("IHW not installed; falling back to BH.")
        res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
      }
    } else { # "none"
      res$p_adj <- res$p_value
    }
    res
  }

  .enrich_side_term <- function(d, side = c("query","subject"), term = "ipr",
                                comp = "comparison", comp_dir = ".", save = TRUE) {
    side <- match.arg(side)
    prefix <- if (side == "query") "q_" else "s_"
    col <- paste0(prefix, tolower(term))
    if (!col %in% names(d) || !is.character(d[[col]])) return(NULL)

    df  <- .apply_filter(d, filter_expr)
    if (!nrow(df)) return(NULL)
    pos <- df[df$dNdS > pos_threshold, , drop = FALSE]

    res <- .fisher_one_term(df[[col]], pos[[col]], drop_empty = drop_rows_without_term)
    if (is.null(res) || !nrow(res)) return(NULL)

    # Rare-term prefilters
    keep <- (res$total_count >= min_total) & (res$pos_count >= min_pos)
    res <- res[keep, , drop = FALSE]
    if (!nrow(res)) return(NULL)

    # FDR adjustment
    res <- .adjust_pvals(res, fdr_method, alpha)

    res$side      <- side
    res$term_type <- tolower(term)

    if (save) {
      out_tsv <- file.path(comp_dir, sprintf("%s_%s_%s_enrichment.tsv",
                                             comp, if (side=="query") "q" else "s", toupper(term)))
      utils::write.table(res[order(res$p_adj, -res$pos_count, res$term), ],
                         file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
      out_svg <- sub("_enrichment.tsv$", "_enrichment_top20.svg", out_tsv)
      .write_plot(res, sprintf("%s (%s)", toupper(term), side), out_svg)
      return(out_tsv)
    }
    res
  }

  .run_one_comparison <- function(comp_name, comp_dir) {
    in_file <- file.path(comp_dir, paste0(comp_name, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) stop("Annotated dN/dS file not found: ", in_file)
    d <- utils::read.table(in_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    term_set <- if (is.null(terms) || !length(terms)) .available_terms(d) else tolower(terms)
    if (!length(term_set)) { message("No term columns found in ", in_file, "; skipping."); return(character(0)) }
    out_paths <- character(0)
    for (tm in term_set) for (sd in sides) {
      p <- .enrich_side_term(d, side = sd, term = tm, comp = comp_name, comp_dir = comp_dir, save = TRUE)
      if (!is.null(p)) out_paths <- c(out_paths, p)
    }
    if (!length(out_paths)) message("No enrichments produced for ", comp_name)
    out_paths
  }

  # ---- batch vs single ----
  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)
    outs <- character(0)
    for (i in seq_len(nrow(df))) {
      comp <- df$comparison_name[i]
      comp_dir <- file.path(output_dir, comp)
      dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
      outs <- c(outs, .run_one_comparison(comp, comp_dir))
    }
    message("All term enrichments complete.")
    return(invisible(outs))
  }

  # single mode
  if (is.null(dnds_annot_file)) stop("Provide either comparison_file (batch) OR dnds_annot_file (single).")
  if (!file.exists(dnds_annot_file)) stop("dnds_annot_file not found: ", dnds_annot_file)
  d <- utils::read.table(dnds_annot_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  comp_dir  <- dirname(dnds_annot_file)
  comp_base <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file))
  comp_name <- comp_base

  term_set <- if (is.null(terms) || !length(terms)) .available_terms(d) else tolower(terms)
  if (!length(term_set)) stop("No term columns present in: ", dnds_annot_file)

  outs <- character(0)
  for (tm in term_set) for (sd in sides) {
    p <- .enrich_side_term(d, side = sd, term = tm, comp = comp_name, comp_dir = comp_dir, save = TRUE)
    if (!is.null(p)) outs <- c(outs, p)
  }
  message("Term enrichment complete for: ", dnds_annot_file)
  invisible(outs)
}
