#' Term enrichment (Fisher test) for query/subject term columns in dN/dS results
#'
#' Reads a `<comparison>/<comparison>_dnds_annot.tsv` file and tests enrichment of
#' annotation terms found in `q_*` / `s_*` columns among positively selected rows
#' (`dNdS > pos_threshold`) versus the filtered background.
#'
#' Term "types" are inferred from suffixes of `q_*`/`s_*` columns (e.g., `q_pfam`,
#' `s_kegg`). By default, `ipr` and `go` are excluded from this generic routine.
#' (Use dedicated functions for IPR/GO if you have them.)
#'
#' For each term type and side, runs a one-sided Fisher exact test ("greater") on
#' term presence/absence counts and reports odds ratio, p-value, and adjusted p-value.
#' Optionally joins metadata (NAME/DEFINITION) and supports term-ID exclusion and
#' descendant exclusion using a parent-child tree.
#'
#' @param dnds_annot_file Path to a single `<comp>_dnds_annot.tsv` file (single mode).
#'   If provided, `comparison_file` must be NULL.
#' @param comparison_file Path to whitespace-delimited file (tabs/spaces; header or not)
#'   with columns: `comparison_name, query_fasta, query_gff, subject_fasta, subject_gff`.
#'   In batch mode, each comparison is read from:
#'   `file.path(output_dir, comparison_name, paste0(comparison_name, "_dnds_annot.tsv"))`.
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#' @param terms Optional character vector of term types to test (suffixes like
#'   `c("kegg","pfam")`). Default NULL = auto-detect from `q_*`/`s_*` columns.
#' @param exclude_terms Character vector of term types to exclude entirely from testing
#'   (applies to auto-detect and explicit `terms`). Default `c("ipr","go")`.
#'   Set to NULL to allow them.
#' @param sides Character vector among `c("query","subject")`. Default both.
#' @param pos_threshold Numeric. `dNdS > pos_threshold` defines "positive" (default 1).
#' @param max_dnds Numeric. Drops rows with `dNdS >= max_dnds` or NA `dNdS` (default 10).
#' @param filter_expr Optional character; a logical expression evaluated in the data
#'   to further filter rows (e.g., `"q_seqname == s_seqname"`). Default NULL.
#' @param make_plots Logical; if TRUE, write a top-N bubble plot per result (default TRUE).
#' @param top_n Integer; number of rows to include in the plot (default 20).
#' @param drop_rows_without_term Logical; if TRUE (default), rows with no term in the
#'   tested column are removed from BOTH positives and background before counting.
#' @param min_total Minimum total occurrences (pos+nonpos) required for a term (default 0).
#' @param min_pos Minimum positive occurrences required for a term (default 0).
#' @param fdr_method One of `"BH"`, `"IHW"`, `"qvalue"`, `"none"`. Default `"BH"`.
#' @param alpha FDR level for IHW weighting (default 0.05).
#' @param term_seps Candidate separators to auto-detect per column
#'   (default `c(";", "|", ",")`).
#' @param term_blocklist Character vector of suffixes to ignore as term families.
#'   (e.g. fields like `id`, `name`, `start`, `end`, etc.) Defaults are conservative.
#' @param exclude_ids Either a character vector of term IDs to exclude globally (all types),
#'   or a named list mapping `term_type -> character vector` of IDs to exclude for that type.
#'   Default NULL.
#' @param term_trees Optional tree(s) enabling descendant exclusion. Either:
#'   (1) a path/data.frame edgelist (two columns: parent, child) applied to all types, or
#'   (2) a named list mapping `term_type -> (path or data.frame edgelist)`.
#' @param exclude_descendants Logical; if TRUE, expands `exclude_ids` per type using
#'   descendants found in `term_trees`. Default FALSE.
#' @param exclude_descendants_depth Integer depth limit when expanding descendants.
#'   `1` = direct children; `Inf` = full subtree. Default `Inf`.
#' @param exclude_descendants_limit Hard cap on number of excluded nodes per type to
#'   avoid accidental mass exclusion. Default 5000.
#' @param term_metadata Optional metadata join(s). Either:
#'   (1) a path/data.frame with columns `ID`, `NAME` (optional `DEFINITION`) applied to all types, or
#'   (2) a named list mapping `term_type -> (path or data.frame)`.
#' @param keep_unmatched_ids Logical; if TRUE (default), keep terms lacking metadata.
#'   If FALSE, rows without a metadata NAME are dropped.
#' @param verbose Logical; if TRUE, prints detection and audit messages. Default TRUE.
#'
#' @return (Invisibly) a character vector of output TSV paths.
#'   In batch mode, this includes all comparisons. In single mode, paths for that file.
#'   TSV columns include: term, term_type, side, pos_count, nonpos_count, pos_total,
#'   nonpos_total, odds_ratio, p_value, p_adj, total_count, enrichment, and optionally
#'   term_name, term_def, label.
#'
#' @export
term_enrichment <- function(dnds_annot_file = NULL,
                            comparison_file = NULL,
                            output_dir = getwd(),
                            terms = NULL,
                            exclude_terms = c("ipr","go"),
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
                            term_seps = c(";", "|", ","),
                            term_blocklist = c(
                              "attributes","attribute","attr","notes","note","description",
                              "product","name","id","gene_id","transcript_id","parent",
                              "dbxref","source","target","type","seqname","seqid",
                              "start","end","strand","phase","biotype","class","len"
                            ),
                            exclude_ids = NULL,
                            term_trees = NULL,
                            exclude_descendants = FALSE,
                            exclude_descendants_depth = Inf,
                            exclude_descendants_limit = 5000,
                            term_metadata = NULL,
                            keep_unmatched_ids = TRUE,
                            verbose = TRUE) {

  sides <- match.arg(sides, choices = c("query","subject"), several.ok = TRUE)
  fdr_method <- match.arg(fdr_method)

  .msg <- function(...) if (isTRUE(verbose)) message(...)

  .read_dnds_annot <- function(path) {
    utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                      quote = "", comment.char = "", check.names = FALSE)
  }

  .coerce_qs_cols_to_char <- function(d) {
    qs <- names(d)[startsWith(names(d), "q_") | startsWith(names(d), "s_")]
    for (nm in qs) d[[nm]] <- as.character(d[[nm]])
    d
  }

  .read_ws <- function(path, header_try = TRUE) {
    utils::read.table(path, header = header_try, sep = "", quote = "\"",
                      stringsAsFactors = FALSE, comment.char = "",
                      strip.white = TRUE, blank.lines.skip = TRUE, check.names = FALSE)
  }

  .read_comparisons <- function(x) {
    req <- c("comparison_name","query_fasta","query_gff","subject_fasta","subject_gff")
    if (is.data.frame(x)) {
      stopifnot(all(req %in% names(x)))
      return(x[, req, drop = FALSE])
    }
    df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
    if (!inherits(df1, "try-error") && all(req %in% names(df1))) return(df1[, req, drop = FALSE])
    df2 <- .read_ws(x, header_try = FALSE)
    if (ncol(df2) < 5) stop("comparison_file must have 5 cols or a header with: ", paste(req, collapse = ", "))
    names(df2)[1:5] <- req
    df2[, req, drop = FALSE]
  }

  .apply_filter <- function(d) {
    keep <- !is.na(d$dNdS) & (d$dNdS < max_dnds)
    if (!is.null(filter_expr) && nzchar(filter_expr)) {
      ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)
      if (!inherits(ok, "try-error")) keep <- keep & isTRUE(as.vector(ok))
    }
    d[keep, , drop = FALSE]
  }

  .best_sep <- function(x, seps) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(seps[[1]])
    x <- utils::head(x, 5000)
    scores <- vapply(seps, function(s) sum(lengths(strsplit(x, s, fixed = TRUE))), numeric(1))
    seps[[which.max(scores)]]
  }

  .split_unique <- function(s, sep) {
    s <- as.character(s)[1]
    if (is.na(s) || !nzchar(s)) return(character(0))
    unique(strsplit(s, sep, fixed = TRUE)[[1]])
  }

  .count_terms <- function(v, sep) {
    v <- as.character(v)
    if (!length(v)) return(integer(0))
    tabs <- table(unlist(lapply(v, .split_unique, sep = sep), use.names = FALSE))
    tabs[names(tabs) != ""]
  }

  .build_set <- function(ids) {
    if (is.null(ids) || !length(ids)) return(NULL)
    ids <- unique(as.character(ids))
    stats::setNames(rep(TRUE, length(ids)), ids)
  }

  .load_df <- function(obj) {
    if (is.null(obj)) return(NULL)
    if (is.character(obj) && length(obj) == 1L && file.exists(obj)) {
      return(utils::read.table(obj, header = TRUE, sep = "\t", quote = "",
                               stringsAsFactors = FALSE, check.names = FALSE))
    }
    if (is.data.frame(obj)) return(obj)
    NULL
  }

  .expand_descendants <- function(seeds, tree_df, depth = Inf, limit = Inf) {
    if (is.null(tree_df) || !length(seeds) || ncol(tree_df) < 2) return(seeds)
    parents  <- as.character(tree_df[[1]])
    children <- as.character(tree_df[[2]])
    adj <- split(children, parents)

    seen <- unique(as.character(seeds))
    frontier <- seen
    curd <- 0L

    while (length(frontier) && curd < depth) {
      kids <- unlist(adj[frontier], use.names = FALSE)
      kids <- unique(kids[!is.na(kids) & nzchar(kids)])
      kids <- setdiff(kids, seen)
      if (!length(kids)) break

      seen <- c(seen, kids)

      if (is.finite(limit) && length(seen) > limit) {
        warning(sprintf("Descendant exclusion capped at %d nodes.", limit), call. = FALSE)
        seen <- unique(seen)[seq_len(limit)]
        break
      }

      frontier <- kids
      curd <- curd + 1L
    }

    unique(seen)
  }

  .prepare_meta_maps <- function(meta_spec) {
    df <- .load_df(meta_spec)
    if (is.null(df) || !"ID" %in% names(df)) return(list(name = NULL, def = NULL))
    nm  <- if ("NAME" %in% names(df)) df$NAME else rep(NA_character_, nrow(df))
    def <- if ("DEFINITION" %in% names(df)) df$DEFINITION else rep(NA_character_, nrow(df))
    list(
      name = stats::setNames(as.character(nm),  as.character(df$ID)),
      def  = stats::setNames(as.character(def), as.character(df$ID))
    )
  }

  .per_type_assets <- function(term_type) {
    term_type <- tolower(term_type)

    ids_global <- if (is.character(exclude_ids)) exclude_ids else NULL
    ids_type <- if (is.list(exclude_ids) && !is.null(exclude_ids[[term_type]])) exclude_ids[[term_type]] else NULL
    seeds <- unique(c(ids_global, ids_type))
    seeds <- seeds[!is.na(seeds) & nzchar(seeds)]

    if (isTRUE(exclude_descendants) && length(seeds)) {
      tree_df <- NULL
      if (is.list(term_trees) && !is.null(term_trees[[term_type]])) tree_df <- .load_df(term_trees[[term_type]])
      else tree_df <- .load_df(term_trees)
      if (!is.null(tree_df)) {
        seeds <- .expand_descendants(seeds, tree_df,
                                     depth = exclude_descendants_depth,
                                     limit = exclude_descendants_limit)
      }
    }

    exclude_set <- .build_set(seeds)

    meta_maps <- NULL
    if (is.list(term_metadata) && !is.null(term_metadata[[term_type]])) meta_maps <- .prepare_meta_maps(term_metadata[[term_type]])
    else meta_maps <- .prepare_meta_maps(term_metadata)

    list(exclude_set = exclude_set, meta_maps = meta_maps)
  }

  .filter_terms_string <- function(s, exclude_set, sep) {
    s <- as.character(s)[1]
    if (is.na(s) || !nzchar(s)) return(s)
    parts <- unique(strsplit(s, sep, fixed = TRUE)[[1]])
    if (!length(parts)) return("")
    if (!is.null(exclude_set)) {
      bad <- as.logical(exclude_set[parts]); bad[is.na(bad)] <- FALSE
      parts <- parts[!bad]
    }
    if (!length(parts)) return("")
    paste(parts, collapse = sep)
  }

  .attach_metadata <- function(res, meta_maps) {
    if (is.null(res) || !nrow(res)) return(res)

    nm <- meta_maps$name
    df <- meta_maps$def

    res$term_name <- if (!is.null(nm)) unname(nm[res$term]) else NA_character_
    if (!is.null(df)) res$term_def <- unname(df[res$term])

    if (isFALSE(keep_unmatched_ids)) {
      keep <- !is.na(res$term_name) & nzchar(res$term_name)
      res <- res[keep, , drop = FALSE]
      if (!nrow(res)) return(res)
    }

    res$label <- ifelse(is.na(res$term_name) | !nzchar(res$term_name),
                        res$term, paste0(res$term, " -- ", res$term_name))
    res
  }

  .adjust_pvals <- function(res) {
    if (!nrow(res)) { res$p_adj <- numeric(0); return(res) }

    if (fdr_method == "BH") {
      res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
    } else if (fdr_method == "qvalue") {
      if (requireNamespace("qvalue", quietly = TRUE)) {
        res$p_adj <- as.numeric(qvalue::qvalue(res$p_value)$qvalues)
      } else {
        warning("qvalue not installed; falling back to BH.", call. = FALSE)
        res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
      }
    } else if (fdr_method == "IHW") {
      if (requireNamespace("IHW", quietly = TRUE)) {
        ihw_obj <- IHW::ihw(res$p_value ~ res$total_count, alpha = alpha)
        res$p_adj <- as.numeric(IHW::adj_pvalues(ihw_obj))
      } else {
        warning("IHW not installed; falling back to BH.", call. = FALSE)
        res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
      }
    } else {
      res$p_adj <- res$p_value
    }

    res
  }

  .fisher_enrichment <- function(vec_all, vec_pos, sep, drop_empty) {
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

    all_tab <- .count_terms(vec_all, sep)
    pos_tab <- .count_terms(vec_pos, sep)
    if (!length(all_tab)) return(NULL)

    res <- do.call(rbind, lapply(names(all_tab), function(tt) {
      a <- as.integer(ifelse(tt %in% names(pos_tab), pos_tab[[tt]], 0L))
      c <- as.integer(all_tab[[tt]] - a)
      b <- n_pos - a
      d <- n_bg  - c

      ft <- try(stats::fisher.test(matrix(c(a,b,c,d), nrow = 2), alternative = "greater"), silent = TRUE)

      data.frame(
        term = tt,
        pos_count = a,
        nonpos_count = c,
        pos_total = n_pos,
        nonpos_total = n_bg,
        odds_ratio = if (inherits(ft, "try-error")) NA_real_ else as.numeric(ft$estimate),
        p_value = if (inherits(ft, "try-error")) NA_real_ else ft$p.value,
        stringsAsFactors = FALSE
      )
    }))

    if (is.null(res) || !nrow(res)) return(NULL)

    res$total_count <- res$pos_count + res$nonpos_count

    pc <- pmax(res$pos_count, 0);  pt <- pmax(res$pos_total, 1)
    nc <- pmax(res$nonpos_count, 0); nt <- pmax(res$nonpos_total, 1)
    bg_rate  <- nc / nt
    pos_rate <- pc / pt
    res$enrichment <- ifelse(bg_rate == 0 & pos_rate > 0, Inf,
                             ifelse(bg_rate == 0 & pos_rate == 0, 1, pos_rate / bg_rate))
    res
  }

  .sanitize_type <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    gsub("_+", "_", gsub("^_|_$", "", x))
  }

  .write_plot <- function(df, ylab, out_svg) {
    if (!isTRUE(make_plots) || !requireNamespace("ggplot2", quietly = TRUE) || !nrow(df)) {
      return(invisible(NULL))
    }

    top <- df[order(df$p_adj, -df$enrichment, df$term), , drop = FALSE]
    top <- utils::head(top, top_n)

    yvar <- if ("label" %in% names(top)) "label" else "term"
    top$y_lab <- top[[yvar]]

    gg <- ggplot2::ggplot(
      top,
      ggplot2::aes(
        x = enrichment,
        y = stats::reorder(y_lab, -p_adj),
        size = pos_count,
        color = p_adj
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradient(low = "red", high = "blue") +
      ggplot2::labs(x = "Enrichment (pos/bg)", y = ylab, size = "# pos", color = "adj p") +
      ggplot2::theme_minimal(base_size = 13)

    ggplot2::ggsave(out_svg, gg, width = 11, height = 9)
    invisible(NULL)
  }

  .detect_term_suffixes <- function(d) {
    qcols <- names(d)[startsWith(names(d), "q_")]
    scols <- names(d)[startsWith(names(d), "s_")]
    cand <- unique(tolower(sub("^q_|^s_", "", c(qcols, scols))))

    cand <- setdiff(cand, tolower(term_blocklist))

    keep <- vapply(cand, function(suf) {
      qn <- paste0("q_", suf); sn <- paste0("s_", suf)
      (qn %in% names(d) && is.character(d[[qn]])) ||
        (sn %in% names(d) && is.character(d[[sn]]))
    }, logical(1))

    cand[keep]
  }

  .resolve_term_set <- function(d) {
    base_detect <- .detect_term_suffixes(d)
    term_set <- base_detect

    if (!is.null(exclude_terms) && length(exclude_terms)) {
      term_set <- setdiff(term_set, tolower(exclude_terms))
    }

    if (!is.null(terms) && length(terms)) {
      term_set <- tolower(terms)
      if (!is.null(exclude_terms) && length(exclude_terms)) {
        drop_now <- intersect(term_set, tolower(exclude_terms))
        if (length(drop_now)) {
          .msg("Dropping excluded term types: ", paste(drop_now, collapse = ", "))
          term_set <- setdiff(term_set, drop_now)
        }
      }
    }

    list(base_detect = base_detect, term_set = term_set)
  }

  .enrich_side_term <- function(d, side, term_type, comp_name, comp_dir) {
    prefix <- if (side == "query") "q_" else "s_"
    term_type <- tolower(term_type)
    col <- paste0(prefix, term_type)

    if (!col %in% names(d) || !is.character(d[[col]])) return(NULL)

    df <- .apply_filter(d)
    if (!nrow(df)) return(NULL)
    pos <- df[df$dNdS > pos_threshold, , drop = FALSE]

    sep_here <- .best_sep(df[[col]], seps = term_seps)

    assets <- .per_type_assets(term_type)
    exset  <- assets$exclude_set
    meta   <- assets$meta_maps

    vec_all <- vapply(df[[col]],  .filter_terms_string, FUN.VALUE = character(1),
                      exclude_set = exset, sep = sep_here)
    vec_pos <- vapply(pos[[col]], .filter_terms_string, FUN.VALUE = character(1),
                      exclude_set = exset, sep = sep_here)

    if (isTRUE(verbose)) {
      all_nonempty <- vec_all[!is.na(vec_all) & nzchar(vec_all)]
      pos_nonempty <- vec_pos[!is.na(vec_pos) & nzchar(vec_pos)]
      uniq_all <- length(.count_terms(all_nonempty, sep_here))
      uniq_pos <- length(.count_terms(pos_nonempty, sep_here))
      .msg(sprintf("[audit] %s | %s | %s :: sep='%s' rows(all_with_term)=%d rows(pos_with_term)=%d uniq(all)=%d uniq(pos)=%d",
                   comp_name, side, toupper(term_type), sep_here,
                   length(all_nonempty), length(pos_nonempty), uniq_all, uniq_pos))
    }

    res <- .fisher_enrichment(vec_all, vec_pos, sep = sep_here, drop_empty = drop_rows_without_term)
    if (is.null(res) || !nrow(res)) return(NULL)

    keep <- (res$total_count >= min_total) & (res$pos_count >= min_pos)
    res <- res[keep, , drop = FALSE]
    if (!nrow(res)) return(NULL)

    res <- .adjust_pvals(res)

    res$side <- side
    res$term_type <- term_type
    res <- .attach_metadata(res, meta)
    res <- res[order(res$p_adj, -res$enrichment, res$term), , drop = FALSE]

    out_tsv <- file.path(comp_dir, sprintf("%s_%s_%s_enrichment.tsv",
                                           comp_name,
                                           if (side == "query") "q" else "s",
                                           toupper(.sanitize_type(term_type))))
    utils::write.table(res, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    if (isTRUE(make_plots)) {
      out_svg <- sub("_enrichment.tsv$", "_enrichment_top20.svg", out_tsv)
      .write_plot(res, sprintf("%s (%s)", toupper(term_type), side), out_svg)
    }

    out_tsv
  }

  .run_one_comparison <- function(comp_name, comp_dir) {
    in_file <- file.path(comp_dir, paste0(comp_name, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) stop("Annotated dN/dS file not found: ", in_file)

    d <- .read_dnds_annot(in_file)
    d <- .coerce_qs_cols_to_char(d)

    det <- .resolve_term_set(d)
    .msg(sprintf("[term_enrichment] %s: detected term types = {%s}; testing = {%s}",
                 comp_name,
                 paste(sort(det$base_detect), collapse = ", "),
                 paste(sort(det$term_set), collapse = ", ")))

    if (!length(det$term_set)) {
      .msg("No eligible term columns found in ", in_file, "; skipping.")
      return(character(0))
    }

    out_paths <- character(0)
    for (tm in det$term_set) {
      for (sd in sides) {
        p <- .enrich_side_term(d, side = sd, term_type = tm, comp_name = comp_name, comp_dir = comp_dir)
        if (!is.null(p)) out_paths <- c(out_paths, p)
      }
    }

    if (!length(out_paths)) .msg("No enrichments produced for ", comp_name)
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

    .msg("All term enrichments complete.")
    return(invisible(outs))
  }

  if (is.null(dnds_annot_file)) stop("Provide either comparison_file (batch) OR dnds_annot_file (single).")
  if (!file.exists(dnds_annot_file)) stop("dnds_annot_file not found: ", dnds_annot_file)

  comp_dir  <- dirname(dnds_annot_file)
  comp_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file))

  # Single-mode: allow running in-place without requiring the directory be named exactly comp_name
  # (but still write outputs beside the input file).
  d <- .read_dnds_annot(dnds_annot_file)
  d <- .coerce_qs_cols_to_char(d)

  det <- .resolve_term_set(d)
  .msg(sprintf("[term_enrichment] %s: detected term types = {%s}; testing = {%s}",
               comp_name,
               paste(sort(det$base_detect), collapse = ", "),
               paste(sort(det$term_set), collapse = ", ")))
  if (!length(det$term_set)) stop("No eligible term columns present in: ", dnds_annot_file)

  outs <- character(0)
  for (tm in det$term_set) {
    for (sd in sides) {
      p <- .enrich_side_term(d, side = sd, term_type = tm, comp_name = comp_name, comp_dir = comp_dir)
      if (!is.null(p)) outs <- c(outs, p)
    }
  }

  .msg("Term enrichment complete for: ", dnds_annot_file)
  invisible(outs)
}
