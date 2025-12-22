#' Term enrichment (Fisher test) for query/subject terms in dNdS results
#'
#' Reads <comp>/<comp>_dnds_annot.tsv and tests enrichment of terms (q_* / s_* columns)
#' among positively selected pairs (dNdS > pos_threshold) vs the filtered background.
#'
#' Supports optional per-term-type ID exclusion, subtree (descendant) exclusion using a
#' provided parent->child edgelist, and metadata joins to append names/definitions.
#'
#' @param dnds_annot_file Path to a single <comp>_dnds_annot.tsv (single mode).
#' @param comparison_file Path to whitespace-delimited file (tabs/spaces; header or not)
#'   with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#'   If provided, batch mode reads:
#'   file.path(output_dir, comparison_name, paste0(comparison_name, "_dnds_annot.tsv")).
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#' @param terms Character vector of term types to test (e.g., c("kegg","pfam")).
#'   Default NULL = auto-detect among {kegg, panther, pfam, tigrfam, cog} and any custom
#'   q_*/s_* term columns present. IPR/GO are excluded by default (see `exclude_terms`).
#' @param exclude_terms Character vector of term **types** to exclude entirely from testing
#'   (auto-detect and explicit `terms=`). Default c("ipr","go"). Set to NULL to allow them.
#' @param sides Character vector among c("query","subject"). Default both.
#' @param pos_threshold Numeric. dNdS > pos_threshold defines "positive" (default 1).
#' @param max_dnds Numeric. Drop rows with dNdS >= max_dnds (default 10) or NA dNdS.
#' @param filter_expr Optional character with a logical expression evaluated in the data.
#' @param make_plots Logical; if TRUE, write a top-N bubble plot per result (default TRUE).
#' @param top_n Integer; number of rows for plot (default 20).
#' @param drop_rows_without_term Logical; if TRUE (default), rows with no term of the tested
#'   type are removed from BOTH positives and background (annotation-aware universe).
#' @param min_total Minimum total occurrences (pos+nonpos) required for a term (default 0).
#' @param min_pos   Minimum positive occurrences required for a term (default 0).
#' @param fdr_method One of "BH","IHW","qvalue","none". Defaults to "BH".
#' @param alpha FDR level for IHW weighting (default 0.05).
#' @param term_sep Default separator for term strings in q_*/s_* columns (default ";").
#' @param term_seps Candidate separators to auto-detect per column (default c(";","|",",")).
#' @param term_blocklist Character vector of suffixes to ignore as term families.
#'
#' @param exclude_ids Either a character vector of term IDs (applies to *all* term types),
#'   or a named list mapping term_type -> character vector of IDs to exclude. Default NULL.
#' @param term_trees Optional named list mapping term_type -> path or data.frame of a
#'   parent/child edgelist (two columns). Enables descendant exclusion for that type.
#' @param exclude_descendants Logical; if TRUE, expand `exclude_ids` per type using
#'   descendants from `term_trees[[type]]`. Default FALSE.
#' @param exclude_descendants_depth Integer depth limit when expanding descendants.
#' @param exclude_descendants_limit Hard cap on number of excluded nodes per (type).
#'
#' @param term_metadata Optional named list mapping term_type -> (path or data.frame)
#'   with columns: ID, NAME (optional: DEFINITION). Joined to results to add names/defs.
#' @param keep_unmatched_ids Logical; if TRUE (default), keep terms lacking metadata.
#'
#' @return In single mode: (invisibly) list of data.frames by side_term (if make_plots=FALSE)
#'         or vector of output TSV paths (if make_plots=TRUE). In batch mode: vector of paths.
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
                            term_sep = ";",
                            term_seps = c(";", "|", ","),
                            term_blocklist = c(
                              "attributes","attribute","attr","notes","note","description",
                              "product","name","id","gene_id","transcript_id","parent",
                              "dbxref","source","target","type","seqname","seqid",
                              "start","end","strand","phase","biotype","class", "len"
                            ),
                            exclude_ids = NULL,
                            term_trees = NULL,
                            exclude_descendants = FALSE,
                            exclude_descendants_depth = Inf,
                            exclude_descendants_limit = 5000,
                            term_metadata = NULL,
                            keep_unmatched_ids = TRUE,
                            dnds_file = NULL) {

  # ---- deprecation shim ----
  if (!is.null(dnds_file) && is.null(dnds_annot_file)) {
    warning("`dnds_file` is deprecated; use `dnds_annot_file` instead.", call. = FALSE)
    dnds_annot_file <- dnds_file
  }

  sides <- match.arg(sides, choices = c("query","subject"), several.ok = TRUE)
  fdr_method <- match.arg(fdr_method)

  # ---- utils ----
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

  # ===== per-column separator detection & split/count helpers =====
  .best_sep <- function(x, seps = c(";", "|", ",")) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(seps[[1]])
    x <- utils::head(x, 5000)
    scores <- vapply(seps, function(s) sum(lengths(strsplit(x, s, fixed = TRUE))), numeric(1))
    seps[[which.max(scores)]]
  }
  .split_terms_unique_sep <- function(s, sep) {
    if (is.null(s)) return(character(0))
    s <- as.character(s)[1]
    if (is.na(s) || !nzchar(s)) return(character(0))
    unique(strsplit(s, sep, fixed = TRUE)[[1]])
  }
  .count_terms_from_vec_sep <- function(v, sep) {
    if (!length(v)) return(integer(0))
    tabs <- table(unlist(lapply(v, .split_terms_unique_sep, sep = sep), use.names = FALSE))
    tabs[names(tabs) != ""]
  }

  # (legacy helpers kept for BC)
  .split_terms_unique <- function(s, sep = term_sep) {
    if (is.null(s)) return(character(0))
    s <- as.character(s)[1]
    if (is.na(s) || !nzchar(s)) return(character(0))
    unique(strsplit(s, sep, fixed = TRUE)[[1]])
  }
  .count_terms_from_vec <- function(v) {
    if (!length(v)) return(integer(0))
    tabs <- table(unlist(lapply(v, .split_terms_unique, sep = term_sep), use.names = FALSE))
    tabs[names(tabs) != ""]
  }

  .adjust_pvals <- function(res, method, alpha) {
    if (!nrow(res)) { res$p_adj <- numeric(0); return(res) }
    if (method == "BH") {
      res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
    } else if (method == "qvalue") {
      if (requireNamespace("qvalue", quietly = TRUE)) {
        res$p_adj <- as.numeric(qvalue::qvalue(res$p_value)$qvalues)
      } else { warning("qvalue not installed; falling back to BH.")
        res$p_adj <- stats::p.adjust(res$p_value, method = "BH") }
    } else if (method == "IHW") {
      if (requireNamespace("IHW", quietly = TRUE)) {
        res$p_adj <- as.numeric(IHW::adj_pvalues(IHW::ihw(res$p_value ~ res$total_count, alpha = alpha)))
      } else { warning("IHW not installed; falling back to BH.")
        res$p_adj <- stats::p.adjust(res$p_value, method = "BH") }
    } else res$p_adj <- res$p_value
    res
  }
  .write_plot <- function(df, ylab, out_svg) {
    if (!make_plots || !requireNamespace("ggplot2", quietly = TRUE) || !nrow(df)) return(invisible(NULL))
    top <- df[order(df$p_adj, -df$enrichment, df$term), , drop = FALSE]
    top <- utils::head(top, top_n)
    yvar <- if ("label" %in% names(top)) "label" else "term"
    top$y_lab <- top[[yvar]]
    gg <- ggplot2::ggplot(top, ggplot2::aes(x = enrichment,
                                            y = stats::reorder(y_lab, -p_adj),
                                            size = pos_count,
                                            color = p_adj)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradient(low = "red", high = "blue") +
      ggplot2::labs(x = "Enrichment (pos/bg)", y = ylab, size = "# pos", color = "adj p") +
      ggplot2::theme_minimal(base_size = 13)
    ggplot2::ggsave(out_svg, gg, width = 11, height = 9)
    invisible(NULL)
  }

  # ---- metadata & tree helpers ----
  .build_set <- function(ids) {
    if (is.null(ids) || !length(ids)) return(NULL)
    ids <- unique(as.character(ids))
    stats::setNames(rep(TRUE, length(ids)), ids)
  }
  .load_df <- function(obj) {
    if (is.null(obj)) return(NULL)
    if (is.character(obj) && length(obj) == 1L && file.exists(obj)) {
      return(utils::read.table(obj, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE))
    }
    if (is.data.frame(obj)) return(obj)
    NULL
  }
  .prepare_meta_maps <- function(meta_spec) {
    df <- .load_df(meta_spec)
    if (is.null(df) || !"ID" %in% names(df)) return(list(name = NULL, def = NULL))
    nm <- if ("NAME" %in% names(df)) df$NAME else rep(NA_character_, nrow(df))
    def <- if ("DEFINITION" %in% names(df)) df$DEFINITION else rep(NA_character_, nrow(df))
    list(name = stats::setNames(as.character(nm), df$ID),
         def  = stats::setNames(as.character(def), df$ID))
  }
  .expand_descendants <- function(seeds, tree_df, depth = Inf, limit = Inf) {
    if (is.null(tree_df) || !length(seeds)) return(seeds)
    parents <- as.character(tree_df[[1]]); children <- as.character(tree_df[[2]])
    adj <- split(children, parents)
    seen <- unique(seeds)
    frontier <- unique(seeds)
    curd <- 0L
    while (length(frontier) && curd < depth) {
      kids <- unlist(adj[frontier], use.names = FALSE)
      kids <- unique(kids[!is.na(kids)])
      kids <- setdiff(kids, seen)
      if (!length(kids)) break
      seen <- c(seen, kids)
      if (is.finite(limit) && length(seen) > limit) {
        warning(sprintf("Descendant exclusion capped at %d nodes.", limit))
        seen <- unique(seen)[seq_len(limit)]
        break
      }
      frontier <- kids; curd <- curd + 1L
    }
    unique(seen)
  }

  # Exclusion-aware filtering of a term string
  .filter_terms_string <- function(s, exclude_set = NULL, sep = term_sep) {
    if (is.na(s) || !nzchar(s)) return(s)
    parts <- unique(strsplit(s, sep, fixed = TRUE)[[1]])
    if (!length(parts)) return("")
    if (!is.null(exclude_set)) {
      ex <- as.logical(exclude_set[parts]); ex[is.na(ex)] <- FALSE
      parts <- parts[!ex]
    }
    if (!length(parts)) return("")
    paste(parts, collapse = sep)
  }

  # Core Fisher with per-column separator
  .fisher_one_term <- function(vec_all, vec_pos, drop_empty = TRUE, sep = ";") {
    vec_all <- as.character(vec_all); vec_pos <- as.character(vec_pos)
    if (drop_empty) {
      vec_all <- vec_all[!is.na(vec_all) & nzchar(vec_all)]
      vec_pos <- vec_pos[!is.na(vec_pos) & nzchar(vec_pos)]
    }
    n_pos <- length(vec_pos); n_all <- length(vec_all); n_bg <- n_all - n_pos
    if (n_pos == 0L || n_all == 0L || n_bg < 0L) return(NULL)

    all_tab <- .count_terms_from_vec_sep(vec_all, sep)
    pos_tab <- .count_terms_from_vec_sep(vec_pos, sep)
    if (length(all_tab) == 0) return(NULL)

    res <- do.call(rbind, lapply(names(all_tab), function(tt) {
      a <- as.integer(ifelse(tt %in% names(pos_tab), pos_tab[[tt]], 0L))
      c <- as.integer(all_tab[[tt]] - a)
      b <- n_pos - a
      d <- n_bg  - c
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

    # Enrichment = (pos_count/pos_total) / (nonpos_count/nonpos_total)
    pc <- pmax(res$pos_count, 0);  pt <- pmax(res$pos_total, 1)
    nc <- pmax(res$nonpos_count, 0); nt <- pmax(res$nonpos_total, 1)
    bg_rate  <- nc / nt
    pos_rate <- pc / pt
    res$enrichment <- ifelse(bg_rate == 0 & pos_rate > 0, Inf,
                             ifelse(bg_rate == 0 & pos_rate == 0, 1, pos_rate / bg_rate))
    res
  }

  # Detect available term types
  .available_terms <- function(df) {
    allowed <- setdiff(c("kegg","panther","pfam","tigrfam","cog","ipr","go"),
                       tolower(term_blocklist))
    have <- function(suf) {
      qn <- paste0("q_", suf); sn <- paste0("s_", suf)
      (qn %in% names(df) && is.character(df[[qn]])) ||
        (sn %in% names(df) && is.character(df[[sn]]))
    }
    tolower(allowed[vapply(allowed, have, logical(1))])
  }
  .detect_term_columns <- function(d) {
    qcols <- names(d)[startsWith(names(d), "q_")]
    scols <- names(d)[startsWith(names(d), "s_")]
    cand  <- unique(tolower(sub("^q_|^s_", "", c(qcols, scols))))
    cand  <- setdiff(cand, tolower(term_blocklist))
    keep <- vapply(cand, function(suf) {
      qn <- paste0("q_", suf); sn <- paste0("s_", suf)
      (qn %in% names(d) && is.character(d[[qn]])) ||
        (sn %in% names(d) && is.character(d[[sn]]))
    }, logical(1))
    cand[keep]
  }

  # Prepare per-type exclusion set and metadata maps
  .per_type_assets <- function(term_type) {
    ids_global <- if (is.character(exclude_ids)) exclude_ids else NULL
    ids_type   <- if (is.list(exclude_ids) && !is.null(exclude_ids[[term_type]])) exclude_ids[[term_type]] else NULL
    seeds <- unique(c(ids_global, ids_type))
    if (isTRUE(exclude_descendants) && length(seeds)) {
      tree_df <- if (is.list(term_trees)) .load_df(term_trees[[term_type]]) else .load_df(term_trees)
      if (!is.null(tree_df) && ncol(tree_df) >= 2) {
        seeds <- .expand_descendants(seeds, tree_df,
                                     depth = exclude_descendants_depth,
                                     limit = exclude_descendants_limit)
      }
    }
    exclude_set <- .build_set(seeds)
    meta_maps <- if (is.list(term_metadata)) .prepare_meta_maps(term_metadata[[term_type]]) else .prepare_meta_maps(term_metadata)
    list(exclude_set = exclude_set, meta_maps = meta_maps)
  }

  .attach_metadata <- function(res, meta_maps) {
    if (is.null(res) || !nrow(res)) return(res)
    nm <- meta_maps$name; df <- meta_maps$def
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

  .sanitize_type <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    gsub("_+", "_", gsub("^_|_$", "", x))
  }

  # ===== audit logger: print to screen only, force flush =====
  .audit_side_term <- function(df, pos, col, sep, comp, side, term) {
    all_nonempty <- as.character(df[[col]]); all_nonempty <- all_nonempty[!is.na(all_nonempty) & nzchar(all_nonempty)]
    pos_nonempty <- as.character(pos[[col]]); pos_nonempty <- pos_nonempty[!is.na(pos_nonempty) & nzchar(pos_nonempty)]
    uniq_all <- length(.count_terms_from_vec_sep(all_nonempty, sep))
    uniq_pos <- length(.count_terms_from_vec_sep(pos_nonempty, sep))
    msg <- sprintf("[audit] %s | %s | %s :: sep='%s' rows(all_with_term)=%d rows(pos_with_term)=%d uniq(all)=%d uniq(pos)=%d",
                   comp, side, toupper(term), sep, length(all_nonempty), length(pos_nonempty), uniq_all, uniq_pos)
    cat(msg, "\n")
    flush.console()
  }

  .enrich_side_term <- function(d, side = c("query","subject"), term = "kegg",
                                comp = "comparison", comp_dir = ".", save = TRUE) {
    side <- match.arg(side)
    prefix <- if (side == "query") "q_" else "s_"
    col <- paste0(prefix, tolower(term))
    if (!col %in% names(d) || !is.character(d[[col]])) return(NULL)

    df  <- .apply_filter(d, filter_expr)
    if (!nrow(df)) return(NULL)
    pos <- df[df$dNdS > pos_threshold, , drop = FALSE]

    # choose best separator for this column
    sep_here <- .best_sep(df[[col]], seps = term_seps)

    assets <- .per_type_assets(tolower(term))
    exset  <- assets$exclude_set
    meta   <- assets$meta_maps

    # Filter out excluded IDs directly in the vectors
    vec_all <- vapply(df[[col]],  .filter_terms_string, FUN.VALUE = character(1),
                      exclude_set = exset, sep = sep_here)
    vec_pos <- vapply(pos[[col]], .filter_terms_string, FUN.VALUE = character(1),
                      exclude_set = exset, sep = sep_here)

    # Audit before Fisher (screen only)
    .audit_side_term(df, pos, col, sep_here, comp, side, term)

    res <- .fisher_one_term(vec_all, vec_pos, drop_empty = drop_rows_without_term, sep = sep_here)
    if (is.null(res) || !nrow(res)) {
      message(sprintf("[term_enrichment] No rows after counting for type=%s side=%s", term, side))
      return(NULL)
    }

    # Drop any terms explicitly excluded
    if (!is.null(exset) && nrow(res)) {
      bad <- as.logical(exset[res$term]); bad[is.na(bad)] <- FALSE
      res <- res[!bad, , drop = FALSE]
      if (!nrow(res)) {
        message(sprintf("[term_enrichment] All rows excluded for type=%s side=%s", term, side))
        return(NULL)
      }
    }

    # Rare-term prefilters
    keep <- (res$total_count >= min_total) & (res$pos_count >= min_pos)
    res <- res[keep, , drop = FALSE]
    if (!nrow(res)) {
      message(sprintf("[term_enrichment] All terms filtered (min_total/min_pos) for type=%s side=%s", term, side))
      return(NULL)
    }

    # FDR + metadata
    res <- .adjust_pvals(res, fdr_method, alpha)
    res$side      <- side
    res$term_type <- tolower(term)
    res <- .attach_metadata(res, meta)

    if (save) {
      safe_term <- toupper(.sanitize_type(term))
      out_tsv <- file.path(comp_dir, sprintf("%s_%s_%s_enrichment.tsv",
                                             comp, if (side=="query") "q" else "s", safe_term))
      utils::write.table(res[order(res$p_adj, -res$enrichment, res$term), ],
                         file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
      out_svg <- sub("_enrichment.tsv$", "_enrichment_top20.svg", out_tsv)
      .write_plot(res, sprintf("%s (%s)", toupper(term), side), out_svg)
      return(out_tsv)
    }
    res
  }

  .coerce_qs_cols_to_char <- function(d) {
    qs <- names(d)[startsWith(names(d), "q_") | startsWith(names(d), "s_")]
    for (nm in qs) d[[nm]] <- as.character(d[[nm]])
    d
  }

  .run_one_comparison <- function(comp_name, comp_dir) {
    in_file <- file.path(comp_dir, paste0(comp_name, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) stop("Annotated dN/dS file not found: ", in_file)
    d <- utils::read.table(in_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    d <- .coerce_qs_cols_to_char(d)

    base_detect <- unique(c(.available_terms(d), .detect_term_columns(d)))
    term_set <- base_detect
    if (!is.null(exclude_terms) && length(exclude_terms)) {
      term_set <- setdiff(term_set, tolower(exclude_terms))
    }
    if (!is.null(terms) && length(terms)) {
      term_set <- tolower(terms)
      if (!is.null(exclude_terms) && length(exclude_terms)) {
        drop_now <- intersect(term_set, tolower(exclude_terms))
        if (length(drop_now)) {
          message("Dropping excluded term types: ", paste(drop_now, collapse = ", "))
          term_set <- setdiff(term_set, drop_now)
        }
      }
    }
    message(sprintf("[term_enrichment] %s: detected term types = {%s}; testing = {%s}",
                    comp_name,
                    paste(sort(base_detect), collapse=", "),
                    paste(sort(term_set), collapse=", ")))
    if (!length(term_set)) { message("No eligible term columns found in ", in_file, "; skipping."); return(character(0)) }

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
  d <- .coerce_qs_cols_to_char(d)
  comp_dir  <- dirname(dnds_annot_file)
  comp_base <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file))
  comp_name <- comp_base

  base_detect <- unique(c(.available_terms(d), .detect_term_columns(d)))
  term_set <- base_detect
  if (!is.null(exclude_terms) && length(exclude_terms)) {
    term_set <- setdiff(term_set, tolower(exclude_terms))
  }
  if (!is.null(terms) && length(terms)) {
    term_set <- tolower(terms)
    if (!is.null(exclude_terms) && length(exclude_terms)) {
      drop_now <- intersect(term_set, tolower(exclude_terms))
      if (length(drop_now)) {
        message("Dropping excluded term types: ", paste(drop_now, collapse = ", "))
        term_set <- setdiff(term_set, drop_now)
      }
    }
  }
  message(sprintf("[term_enrichment] %s: detected term types = {%s}; testing = {%s}",
                  comp_name,
                  paste(sort(base_detect), collapse=", "),
                  paste(sort(term_set), collapse=", ")))
  if (!length(term_set)) stop("No eligible term columns present in: ", dnds_annot_file)

  outs <- character(0)
  for (tm in term_set) for (sd in sides) {
    p <- .enrich_side_term(d, side = sd, term = tm, comp = comp_name, comp_dir = comp_dir, save = TRUE)
    if (!is.null(p)) outs <- c(outs, p)
  }
  message("Term enrichment complete for: ", dnds_annot_file)
  invisible(outs)
}
