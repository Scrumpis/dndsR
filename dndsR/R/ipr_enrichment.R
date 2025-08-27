#' InterPro (IPR) enrichment (Fisher) for query/subject terms in dNdS results
#'
#' Reads <comp>/<comp>_dnds_annot.tsv and tests enrichment of IPR terms (q_ipr / s_ipr)
#' among positively selected pairs (dNdS > pos_threshold) vs the filtered background.
#' Adds InterPro metadata (ENTRY_TYPE/ENTRY_NAME) and supports pooled vs per-type analyses.
#'
#' @param dnds_annot_file Path to a single <comp>_dnds_annot.tsv (single mode).
#' @param comparison_file Path to whitespace-delimited file (tabs/spaces; header or not)
#'   with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#'   If provided, batch mode reads:
#'   file.path(output_dir, comparison_name, paste0(comparison_name, "_dnds_annot.tsv")).
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#' @param sides Character vector among c("query","subject"). Default both.
#' @param pos_threshold Numeric. dNdS > pos_threshold defines "positive" (default 1).
#' @param max_dnds Numeric. Drop rows with dNdS >= max_dnds (default 10) or NA dNdS.
#' @param filter_expr Optional character with a logical expression evaluated in the data
#'   (e.g., "q_seqname == s_seqname").
#' @param make_plots Logical; if TRUE, write a top-N bubble plot per result (default TRUE).
#' @param top_n Integer; number of rows for plot (default 20).
#' @param drop_rows_without_term Logical; if TRUE (default), rows with no IPR term are
#'   removed from BOTH positives and background (annotation-aware universe).
#' @param min_total Minimum total occurrences (pos+nonpos) required for a term (default 0).
#' @param min_pos   Minimum positive occurrences required for a term (default 0).
#' @param fdr_method One of "BH","IHW","qvalue","none". Defaults to "BH".
#'                   If "IHW"/"qvalue" is chosen but the package is missing, falls back to "BH".
#' @param alpha FDR level for IHW weighting (default 0.05).
#'
#' @param include_types Optional character vector of InterPro ENTRY_TYPE values to keep
#'   in pooled mode (e.g., c("Domain","Homologous_superfamily")). Ignored if stratified.
#' @param stratify_by_type Logical. If TRUE, run separate analyses per ENTRY_TYPE with
#'   type-specific backgrounds (recommended for apples-to-apples). Default FALSE.
#' @param types Character vector of ENTRY_TYPEs to analyze when stratified; default NULL
#'   = infer from data present in q_ipr/s_ipr.
#' @param adjust_scope When stratified, "global" (one BH across all tests) or "per_type".
#'
#' @param entries_source Where to load InterPro entries from:
#'   "auto" (bundled static; default), "local" (use entries_path),
#'   "remote" (download; falls back to bundled), or "none" (skip metadata).
#' @param entries_path Optional path to a TSV copy of InterPro entry.list
#'   (columns: ENTRY_AC, ENTRY_TYPE, ENTRY_NAME). Used for "auto"/"local".
#' @param entries_url  Remote URL for current InterPro entry.list (default EBI).
#' @param entries_timeout_s Numeric timeout (seconds) for remote fetch (default 20).
#' @param keep_unmatched Keep IPRs not found in entry.list when filtering (default TRUE).
#'
#' @return In single mode: (invisibly) list of output TSV paths.
#'         In batch mode:  (invisibly) vector of output TSV paths across comparisons.
#'         Each TSV includes: IPR, ENTRY_TYPE, ENTRY_NAME, pos_count, nonpos_count,
#'         pos_total, nonpos_total, odds_ratio, p_value, p_adj, total_count, label,
#'         side, comparison.
#' @export
ipr_enrichment <- function(dnds_annot_file = NULL,
                           comparison_file = NULL,
                           output_dir = getwd(),
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
                           include_types = NULL,
                           stratify_by_type = FALSE,
                           types = NULL,
                           adjust_scope = c("global","per_type"),
                           entries_source = c("auto","local","remote","none"),
                           entries_path   = NULL,
                           entries_url    = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list",
                           entries_timeout_s = 20,
                           keep_unmatched = TRUE) {

  # ---- arg setup ----
  sides <- match.arg(sides, choices = c("query","subject"), several.ok = TRUE)
  fdr_method <- match.arg(fdr_method)
  adjust_scope <- match.arg(adjust_scope)
  entries_source <- match.arg(entries_source)

  # ---- helpers (copied/adapted from term_enrichment) ----
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

  .apply_filter <- function(d, filter_expr) {
    keep <- !is.na(d$dNdS) & d$dNdS < max_dnds
    if (!is.null(filter_expr) && nzchar(filter_expr)) {
      ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)
      if (!inherits(ok, "try-error")) keep <- keep & isTRUE(as.vector(ok))
    }
    d[keep, , drop = FALSE]
  }

  .split_terms_unique <- function(x) {
    if (is.null(x)) return(character(0))
    s <- as.character(x)[1]
    if (is.na(s) || !nzchar(s)) return(character(0))
    unique(strsplit(s, ";", fixed = TRUE)[[1]])
  }

  .count_terms_from_vec <- function(v) {
    if (!length(v)) return(integer(0))
    tabs <- table(unlist(lapply(v, .split_terms_unique), use.names = FALSE))
    tabs[names(tabs) != ""]
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

  .write_plot <- function(df, ylab, out_svg) {
    if (!make_plots || !requireNamespace("ggplot2", quietly = TRUE) || !nrow(df)) return(invisible(NULL))
    top <- df[order(df$p_adj, -df$pos_count, df$IPR), , drop = FALSE]
    top <- utils::head(top, top_n)
    yvar <- if ("label" %in% names(top)) "label" else "IPR"
    gg <- ggplot2::ggplot(top, ggplot2::aes(x = pos_count,
                                            y = stats::reorder(top[[yvar]], -top$p_adj),
                                            size = pos_count / (pos_total + nonpos_total),
                                            color = p_adj)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradient(low = "red", high = "blue") +
      ggplot2::labs(x = "Positive count", y = ylab, size = "Pos / All", color = "adj p") +
      ggplot2::theme_minimal(base_size = 13)
    ggplot2::ggsave(out_svg, gg, width = 11, height = 9)
    invisible(NULL)
  }

  # ---- InterPro entry.list I/O ----
  .read_entries <- function(path) {
    utils::read.table(path, header = TRUE, sep = "\t", quote = "",
                      stringsAsFactors = FALSE, comment.char = "",
                      check.names = FALSE)
  }

  .bundled_path <- function() {
    p <- system.file("extdata", "interpro_entry.list.tsv", package = "dndsR")
    if (!nzchar(p) || !file.exists(p)) return(NULL)
    p
  }

  .load_entries <- function() {
    if (entries_source == "none") return(NULL)

    if (!is.null(entries_path) && file.exists(entries_path) &&
        entries_source %in% c("auto","local")) {
      return(.read_entries(entries_path))
    }

    if (entries_source == "local") {
      if (is.null(entries_path) || !file.exists(entries_path)) {
        stop("entries_source='local' but entries_path is missing or does not exist.")
      }
      return(.read_entries(entries_path))
    }

    if (entries_source == "auto") {
      bp <- .bundled_path()
      if (!is.null(bp)) return(.read_entries(bp))
      # fall through to remote if bundled missing
    }

    # remote
    tf <- tempfile(fileext = ".tsv")
    old_timeout <- getOption("timeout"); on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = max(getOption("timeout"), entries_timeout_s))
    ok <- try(utils::download.file(entries_url, tf, quiet = TRUE), silent = TRUE)
    if (!inherits(ok, "try-error") && isTRUE(ok == 0)) {
      return(.read_entries(tf))
    }

    if (entries_source == "remote") {
      bp <- .bundled_path()
      if (!is.null(bp)) return(.read_entries(bp))
      stop("Failed to load InterPro entry.list from remote and no bundled copy found.")
    }
    NULL
  }

  .build_ipr_maps <- function(entries) {
    if (is.null(entries)) return(list(type_by_ipr = NULL, name_by_ipr = NULL))
    ac <- as.character(entries$ENTRY_AC)
    ty <- as.character(entries$ENTRY_TYPE)
    nm <- as.character(entries$ENTRY_NAME)
    list(type_by_ipr = stats::setNames(ty, ac),
         name_by_ipr = stats::setNames(nm, ac))
  }

  .filter_terms_string_by_types <- function(s, allowed_types, type_by_ipr, keep_unmatched = TRUE) {
    if (is.na(s) || !nzchar(s)) return(s)
    parts <- unique(strsplit(s, ";", fixed = TRUE)[[1]])
    if (!length(parts)) return("")
    keep <- logical(length(parts))
    for (i in seq_along(parts)) {
      ipr <- parts[i]
      ty  <- if (!is.null(type_by_ipr)) type_by_ipr[[ipr]] else NA_character_
      if (is.na(ty)) {
        keep[i] <- isTRUE(keep_unmatched) # optionally keep if unknown
      } else {
        keep[i] <- ty %in% allowed_types
      }
    }
    parts2 <- parts[keep]
    if (!length(parts2)) return("")
    paste(parts2, collapse = ";")
  }

  .sanitize_type <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    gsub("_+", "_", gsub("^_|_$", "", x))
  }

  # ---- core fisher from vectors (after any filtering) ----
  .fisher_from_vectors <- function(vec_all, vec_pos, drop_empty) {
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

    all_tab <- .count_terms_from_vec(vec_all)
    pos_tab <- .count_terms_from_vec(vec_pos)
    if (length(all_tab) == 0) return(NULL)

    res <- do.call(rbind, lapply(names(all_tab), function(tt) {
      a <- as.integer(ifelse(tt %in% names(pos_tab), pos_tab[[tt]], 0L))     # pos with term
      c <- as.integer(all_tab[[tt]] - a)                                     # non-pos with term
      b <- n_pos - a                                                         # pos without term
      d <- n_bg  - c                                                         # non-pos without term
      ft <- try(stats::fisher.test(matrix(c(a,b,c,d), nrow = 2),
                                   alternative = "greater"), silent = TRUE)
      data.frame(
        IPR          = tt,
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

  # ---- InterPro metadata attach ----
  .attach_metadata_to_res <- function(res, type_by_ipr, name_by_ipr) {
    if (is.null(res) || !nrow(res)) return(res)
    if (is.null(type_by_ipr) && is.null(name_by_ipr)) {
      res$ENTRY_TYPE <- NA_character_
      res$ENTRY_NAME <- NA_character_
    } else {
      ty <- type_by_ipr[res$IPR]; nm <- name_by_ipr[res$IPR]
      res$ENTRY_TYPE <- unname(ifelse(is.na(ty), NA_character_, ty))
      res$ENTRY_NAME <- unname(ifelse(is.na(nm), NA_character_, nm))
    }
    res$label <- ifelse(is.na(res$ENTRY_NAME) | !nzchar(res$ENTRY_NAME),
                        res$IPR,
                        paste0(res$IPR, " — ", res$ENTRY_NAME))
    res
  }

  # ---- one side, one (possibly filtered) IPR vector ----
  .enrich_side_ipr <- function(d, side, comp, comp_dir,
                               type_by_ipr, name_by_ipr,
                               pooled_allowed_types = NULL,
                               stratified_types = NULL) {

    prefix <- if (side == "query") "q_" else "s_"
    col <- paste0(prefix, "ipr")
    if (!col %in% names(d) || !is.character(d[[col]])) return(character(0))

    df <- .apply_filter(d, filter_expr)
    if (!nrow(df)) return(character(0))
    pos <- df[df$dNdS > pos_threshold, , drop = FALSE]

    # --- STRATIFIED mode: per ENTRY_TYPE with type-specific background
    if (!is.null(stratified_types)) {
      out_paths <- character(0)
      # Prepare a filtered copy per type
      for (tp in stratified_types) {
        # Filter each row's term-string to only this type
        vec_all <- vapply(df[[col]], .filter_terms_string_by_types,
                          FUN.VALUE = character(1),
                          allowed_types = tp, type_by_ipr = type_by_ipr,
                          keep_unmatched = keep_unmatched)
        vec_pos <- vapply(pos[[col]], .filter_terms_string_by_types,
                          FUN.VALUE = character(1),
                          allowed_types = tp, type_by_ipr = type_by_ipr,
                          keep_unmatched = keep_unmatched)

        res <- .fisher_from_vectors(vec_all, vec_pos, drop_rows_without_term)
        if (is.null(res) || !nrow(res)) next

        # rare-term prefilters
        keep <- (res$total_count >= min_total) & (res$pos_count >= min_pos)
        res <- res[keep, , drop = FALSE]
        if (!nrow(res)) next

        res <- .attach_metadata_to_res(res, type_by_ipr, name_by_ipr)
        res$ENTRY_TYPE <- tp

        # We'll set p_adj later (global/per_type) in the caller
        # stash the slice by type
        out_paths <- c(out_paths, list(res))
      }
      return(out_paths)  # list of data.frames to be merged at caller
    }

    # --- POOLED mode (optionally filter to a set of types)
    if (!is.null(pooled_allowed_types)) {
      vec_all <- vapply(df[[col]], .filter_terms_string_by_types,
                        FUN.VALUE = character(1),
                        allowed_types = pooled_allowed_types,
                        type_by_ipr = type_by_ipr,
                        keep_unmatched = keep_unmatched)
      vec_pos <- vapply(pos[[col]], .filter_terms_string_by_types,
                        FUN.VALUE = character(1),
                        allowed_types = pooled_allowed_types,
                        type_by_ipr = type_by_ipr,
                        keep_unmatched = keep_unmatched)
    } else {
      vec_all <- as.character(df[[col]])
      vec_pos <- as.character(pos[[col]])
    }

    res <- .fisher_from_vectors(vec_all, vec_pos, drop_rows_without_term)
    if (is.null(res) || !nrow(res)) return(character(0))

    keep <- (res$total_count >= min_total) & (res$pos_count >= min_pos)
    res <- res[keep, , drop = FALSE]
    if (!nrow(res)) return(character(0))

    res <- .attach_metadata_to_res(res, type_by_ipr, name_by_ipr)
    res
  }

  # ---- entries + maps ----
  entries <- .load_entries()
  maps <- .build_ipr_maps(entries)
  type_by_ipr <- maps$type_by_ipr
  name_by_ipr <- maps$name_by_ipr

  # ---- runner per comparison ----
  .run_one_comparison <- function(comp_name, comp_dir) {
    in_file <- file.path(comp_dir, paste0(comp_name, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) stop("Annotated dN/dS file not found: ", in_file)
    d <- utils::read.table(in_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

    outs <- character(0)

    # STRATIFIED: figure out which types to run
    strat_types <- NULL
    if (isTRUE(stratify_by_type)) {
      if (is.null(type_by_ipr)) stop("stratify_by_type=TRUE requires InterPro metadata (entries_source != 'none').")
      if (is.null(types)) {
        # infer from IPRs present in q_ipr and s_ipr
        all_terms <- character(0)
        if ("q_ipr" %in% names(d) && is.character(d$q_ipr)) {
          all_terms <- c(all_terms, unlist(strsplit(paste(na.omit(d$q_ipr)), ";", fixed = TRUE)))
        }
        if ("s_ipr" %in% names(d) && is.character(d$s_ipr)) {
          all_terms <- c(all_terms, unlist(strsplit(paste(na.omit(d$s_ipr)), ";", fixed = TRUE)))
        }
        all_terms <- unique(all_terms[nzchar(all_terms)])
        strat_types <- sort(unique(type_by_ipr[all_terms]))
        strat_types <- strat_types[!is.na(strat_types)]
      } else {
        strat_types <- intersect(types, unique(type_by_ipr))
      }
      if (!length(strat_types)) {
        message("No InterPro ENTRY_TYPEs detected in ", in_file, "; skipping stratified analysis.")
        return(character(0))
      }
    }

    for (sd in sides) {
      if (isTRUE(stratify_by_type)) {
        # Collect slices, combine, then FDR (global/per_type) and write per-type files
        slices <- .enrich_side_ipr(d, side = sd, comp = comp_name, comp_dir = comp_dir,
                                   type_by_ipr = type_by_ipr, name_by_ipr = name_by_ipr,
                                   stratified_types = strat_types)
        slices <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, slices)
        if (!length(slices)) next
        allres <- do.call(rbind, slices)
        if (!nrow(allres)) next

        # FDR scope
        if (adjust_scope == "global") {
          allres$p_adj <- stats::p.adjust(allres$p_value, method = "BH")
        } else {
          allres$p_adj <- NA_real_
          for (tp in unique(allres$ENTRY_TYPE)) {
            idx <- which(allres$ENTRY_TYPE == tp)
            allres$p_adj[idx] <- stats::p.adjust(allres$p_value[idx], method = "BH")
          }
        }

        # annotate and write per type
        allres$side <- sd
        allres$comparison <- comp_name

        for (tp in unique(allres$ENTRY_TYPE)) {
          sub <- allres[allres$ENTRY_TYPE == tp, , drop = FALSE]
          if (!nrow(sub)) next
          sub <- sub[order(sub$p_adj, -sub$pos_count, sub$IPR), , drop = FALSE]
          out_tsv <- file.path(comp_dir, sprintf("%s_%s_IPR_%s_enrichment.tsv",
                                                 comp_name, if (sd=="query") "q" else "s", .sanitize_type(tp)))
          utils::write.table(sub, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
          out_svg <- sub("_enrichment.tsv$", "_enrichment_top20.svg", out_tsv)
          .write_plot(sub, sprintf("IPR (%s, %s)", tp, sd), out_svg)
          outs <- c(outs, out_tsv)
        }

      } else {
        # pooled run (optionally filter to include_types)
        pooled_types <- include_types
        if (!is.null(pooled_types) && is.null(type_by_ipr)) {
          stop("include_types requires InterPro metadata (entries_source != 'none').")
        }
        res <- .enrich_side_ipr(d, side = sd, comp = comp_name, comp_dir = comp_dir,
                                type_by_ipr = type_by_ipr, name_by_ipr = name_by_ipr,
                                pooled_allowed_types = pooled_types)
        if (is.character(res)) next  # nothing to write
        if (is.null(res) || !nrow(res)) next

        # FDR + annotate + write
        res <- .adjust_pvals(res, fdr_method, alpha)
        res$side <- sd
        res$comparison <- comp_name
        res <- res[order(res$p_adj, -res$pos_count, res$IPR), , drop = FALSE]

        out_tsv <- file.path(comp_dir, sprintf("%s_%s_IPR_enrichment.tsv",
                                               comp_name, if (sd=="query") "q" else "s"))
        utils::write.table(res, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
        out_svg <- sub("_enrichment.tsv$", "_enrichment_top20.svg", out_tsv)
        .write_plot(res, sprintf("IPR (%s)", sd), out_svg)
        outs <- c(outs, out_tsv)
      }
    }

    if (!length(outs)) message("No IPR enrichments produced for ", comp_name)
    outs
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
    message("All IPR enrichments complete.")
    return(invisible(outs))
  }

  # single mode
  if (is.null(dnds_annot_file)) stop("Provide either comparison_file (batch) OR dnds_annot_file (single).")
  if (!file.exists(dnds_annot_file)) stop("dnds_annot_file not found: ", dnds_annot_file)
  d <- utils::read.table(dnds_annot_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  comp_dir  <- dirname(dnds_annot_file)
  comp_base <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file))
  comp_name <- comp_base

  outs <- .run_one_comparison(comp_name, comp_dir)
  message("IPR enrichment complete for: ", dnds_annot_file)
  invisible(outs)
}
