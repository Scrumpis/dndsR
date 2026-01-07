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
#' @param min_total Minimum total occurrences (pos+nonpos) required for a term (default 2).
#' @param min_pos   Minimum positive occurrences required for a term (default 2).
#' @param fdr_method One of "BH","BY","IHW","qvalue","none". Defaults to "BH".
#' @param alpha FDR level for IHW weighting (default 0.05).
#' @param term_sep Separator used in q_ipr/s_ipr strings (default ";").
#'
#' @param include_types Optional character vector of InterPro ENTRY_TYPE values to keep
#'   in pooled mode (e.g., c("Domain","Homologous_superfamily")). Ignored if stratified.
#' @param stratify_by_type Logical. If TRUE, run separate analyses per ENTRY_TYPE with
#'   type-specific backgrounds (recommended). Default TRUE.
#' @param types Character vector of ENTRY_TYPEs to analyze when stratified; default NULL
#'   = infer from data present in q_ipr/s_ipr (after exclusions).
#' @param adjust_scope When stratified, "global" (one BH across all tests) or "per_type" (not implemented currently).
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
#' @param exclude_ids Character vector of IPR accessions (e.g. "IPR000123") to exclude
#'   globally from both positives and background before enrichment. Preferred name.
#' @param exclude_iprs Deprecated alias of `exclude_ids` for backward compatibility.
#' @param max_prop Maximum allowed proportion of the background annotated set
#'   that may be assigned to a term (default 0.20). Terms with
#'   (pos_count + nonpos_count) / (pos_total + nonpos_total) > max_prop
#'   are filtered out to avoid overly broad terms dominating enrichment.
#'
#' @param term_trees Optional path or data.frame of a parent/child edgelist for IPRs.
#'   (High-priority local override if provided; two columns parent,child).
#' @param exclude_descendants If TRUE, expand `exclude_ids` using tree. Default FALSE.
#' @param exclude_descendants_depth Integer depth limit for descendant exclusion;
#'   1 = direct children only, Inf = entire subtree (default Inf).
#' @param exclude_descendants_limit Hard cap on the number of excluded IPRs to avoid
#'   accidental mass exclusion (default 5000).
#'
#' @param method Enrichment mode: "fisher" (default), "parent_child" (elim-like), "weight01".
#' @param ancestor_novel_frac Keep parent only if >= this fraction of unclaimed genes (default 0.20).
#'
#' @param x_axis_min,x_axis_max Optional fixed x-axis limits for plots.
#' @param x_axis_pad_right Numeric. Extra units to extend the x-axis beyond the
#'   maximum observed enrichment when x_axis_max is not set (default 3).
#'
#' @param interpro_release Optional string (e.g. "94.0") to pin archived InterPro release.
#' @param auto_detect_release Logical; auto-pick best archived release by IPR coverage (default TRUE).
#' @param release_candidates Character vector of releases to test (default c("92.0","94.0","96.0","98.0","100.0")).
#' @param strict_coverage Numeric in \eqn{[0,1]}; fail if best coverage < this (default 0.90).
#' @param tree_source Where to load ParentChildTree from: "auto","local","remote","none".
#' @param tree_path Optional local ParentChildTreeFile.txt or 2-col TSV edgelist.
#' @param tree_url Optional explicit URL to ParentChildTree (raw txt or TSV).
#' @param tree_timeout_s Numeric timeout (seconds) for tree fetch (default 20).
#'
#' @return In single mode: (invisibly) list of output TSV paths.
#'         In batch mode:  (invisibly) vector of output TSV paths across comparisons.
#'         Each TSV includes: IPR, ENTRY_TYPE, ENTRY_NAME, pos_count, nonpos_count,
#'         pos_total, nonpos_total, odds_ratio, p_value, p_adj, total_count, enrichment,
#'         label, side, comparison.
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
                           min_total = 2,
                           min_pos = 2,
                           max_prop = 0.20,
                           fdr_method = c("BH","BY","IHW","qvalue","none"),
                           alpha = 0.05,
                           term_sep = ";",
                           include_types = NULL,
                           stratify_by_type = TRUE,
                           types = NULL,
                           adjust_scope = c("global","per_type"),
                           entries_source = c("auto","local","remote","none"),
                           entries_path   = NULL,
                           entries_url    = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list",
                           entries_timeout_s = 20,
                           keep_unmatched = TRUE,
                           exclude_ids = NULL,
                           exclude_iprs = NULL,
                           term_trees = NULL,
                           exclude_descendants = FALSE,
                           exclude_descendants_depth = Inf,
                           exclude_descendants_limit = 5000,
                           method = c("fisher","parent_child","weight01"),
                           ancestor_novel_frac = 0.20,
                           x_axis_min = NULL,
                           x_axis_max = NULL,
                           x_axis_pad_right = 3,
                           interpro_release = NULL,
                           auto_detect_release = TRUE,
                           release_candidates = c("92.0","94.0","96.0","98.0","100.0"),
                           strict_coverage = 0.90,
                           tree_source = c("auto","local","remote","none"),
                           tree_path = NULL,
                           tree_url = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt",
                           tree_timeout_s = 20) {

  sides <- match.arg(sides, choices = c("query","subject"), several.ok = TRUE)
  fdr_method <- match.arg(fdr_method)
  adjust_scope <- match.arg(adjust_scope)
  entries_source <- match.arg(entries_source)
  tree_source <- match.arg(tree_source)
  method <- match.arg(method)

  # ---- normalize possibly-NA flags from CLI ----
  .normalize_flag <- function(x, default) if (is.na(x)) default else x
  make_plots             <- .normalize_flag(make_plots,             TRUE)
  drop_rows_without_term <- .normalize_flag(drop_rows_without_term, TRUE)
  stratify_by_type       <- .normalize_flag(stratify_by_type,       TRUE)
  keep_unmatched         <- .normalize_flag(keep_unmatched,         TRUE)
  exclude_descendants    <- .normalize_flag(exclude_descendants,    FALSE)

  # Back-compat alias
  if (!is.null(exclude_iprs) && is.null(exclude_ids)) {
    warning("`exclude_iprs` is deprecated; use `exclude_ids`.", call. = FALSE)
    exclude_ids <- exclude_iprs
  }

  # -------------------------- NEW: InterPro release & tree helpers ----------------
  .release_path <- function(rel) sprintf("https://ftp.ebi.ac.uk/pub/databases/interpro/releases/%s", rel)
  .release_urls <- function(rel) {
    b <- .release_path(rel)
    list(entry_list = sprintf("%s/entry.list", b),
         parent_child = sprintf("%s/ParentChildTreeFile.txt", b),
         notes = sprintf("%s/release_notes.txt", b))
  }
  .current_urls <- list(
    entry_list   = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list",
    parent_child = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt",
    notes        = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/release_notes.txt"
  )
  .download_quiet <- function(url, file, timeout_s = 30) {
    old <- getOption("timeout"); on.exit(options(timeout = old), add = TRUE)
    options(timeout = max(getOption("timeout"), timeout_s))
    ok <- try(utils::download.file(url, file, quiet = TRUE, mode = "wb"), silent = TRUE)
    if (inherits(ok, "try-error") || (!is.null(ok) && ok != 0)) return(FALSE)
    file.exists(file) && file.info(file)$size > 0
  }
  .read_tsv_strict <- function(path) {
    utils::read.table(path, header = TRUE, sep = "\t", quote = "",
                      stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
  }

  # ---- NEW: entry.list reader that supports both headered TSV (your bundled file)
  #          and headerless FTP entry.list format. Returns standardized columns. ----
  .read_interpro_entry_list <- function(path) {
    first <- readLines(path, n = 1, warn = FALSE)
    has_header <- length(first) == 1L &&
      grepl("^\\s*ENTRY_AC\\tENTRY_TYPE\\tENTRY_NAME\\s*$", first)

    if (has_header) {
      df <- utils::read.table(
        path, header = TRUE, sep = "\t", quote = "",
        stringsAsFactors = FALSE, comment.char = "", check.names = FALSE
      )
    } else {
      df <- utils::read.table(
        path, header = FALSE, sep = "\t", quote = "",
        stringsAsFactors = FALSE, comment.char = "", check.names = FALSE
      )
      if (ncol(df) < 3) stop("InterPro entry.list parse error: expected >= 3 tab-separated columns.")
      df <- df[, 1:3, drop = FALSE]
      names(df) <- c("ENTRY_AC", "ENTRY_TYPE", "ENTRY_NAME")
    }

    req <- c("ENTRY_AC","ENTRY_TYPE","ENTRY_NAME")
    if (!all(req %in% names(df))) {
      stop("InterPro entry.list missing required columns: ",
           paste(setdiff(req, names(df)), collapse = ", "))
    }
    df
  }

  .collect_iprs_from_df <- function(d, term_sep = ";") {
    cols <- intersect(c("q_ipr","s_ipr"), names(d))
    if (!length(cols)) return(character(0))
    vals <- unlist(lapply(cols, function(cn) na.omit(as.character(d[[cn]]))), use.names = FALSE)
    parts <- unlist(strsplit(vals, term_sep, fixed = TRUE), use.names = FALSE)
    unique(parts[nzchar(parts)])
  }
  .coverage_for_release <- function(iprs, release, timeout_s = 30) {
    u <- .release_urls(release)$entry_list
    tf <- tempfile(fileext = ".tsv")
    if (!.download_quiet(u, tf, timeout_s)) return(NA_real_)
    el <- try(.read_interpro_entry_list(tf), silent = TRUE)
    if (inherits(el, "try-error") || !"ENTRY_AC" %in% names(el)) return(NA_real_)
    if (!length(iprs)) return(1.0)
    round(sum(iprs %in% el$ENTRY_AC) / length(iprs), 6)
  }
  .pick_best_release <- function(iprs, candidates, timeout_s = 30) {
    if (!length(candidates)) return(list(release = NULL, coverage = NA_real_))
    cov <- vapply(candidates, function(r) .coverage_for_release(iprs, r, timeout_s), numeric(1))
    best <- which.max(replace(cov, is.na(cov), -Inf))
    list(release = candidates[best], coverage = cov[best])
  }
  .load_interpro_entry_list <- function(entries_source = c("auto","local","remote","none"),
                                        entries_path = NULL,
                                        entries_url  = NULL,
                                        interpro_release = NULL,
                                        timeout_s = 30) {
    entries_source <- match.arg(entries_source)
    if (entries_source == "none") return(list(df = NULL, provenance = list()))
    if (!is.null(entries_path) && file.exists(entries_path) &&
        entries_source %in% c("auto","local"))
      return(list(df = .read_interpro_entry_list(entries_path),
                  provenance = list(mode = "local", path = normalizePath(entries_path, winslash = "/"))))
    if (entries_source == "local")
      stop("entries_source='local' but entries_path is missing or not found.")
    if (!is.null(entries_url)) {
      tf <- tempfile(fileext = ".tsv")
      if (.download_quiet(entries_url, tf, timeout_s))
        return(list(df = .read_interpro_entry_list(tf),
                    provenance = list(mode = "remote-url", url = entries_url)))
    }
    url <- if (!is.null(interpro_release)) .release_urls(interpro_release)$entry_list else .current_urls$entry_list
    tf <- tempfile(fileext = ".tsv")
    if (.download_quiet(url, tf, timeout_s))
      return(list(df = .read_interpro_entry_list(tf),
                  provenance = list(mode = if (is.null(interpro_release)) "remote-current" else "remote-release",
                                    url = url, release = interpro_release)))
    if (entries_source == "auto") {
      bp <- system.file("extdata", "interpro_entry.list.tsv", package = "dndsR")
      if (nzchar(bp) && file.exists(bp))
        return(list(df = .read_interpro_entry_list(bp),
                    provenance = list(mode = "bundled", path = normalizePath(bp, winslash = "/"))))
    }
    stop("Failed to load InterPro entry.list from remote/bundled sources.")
  }
  .parse_parent_child_tree_txt <- function(path) {
    lines <- readLines(path, warn = FALSE)
    rx <- "(IPR\\d{6,})"
    keep_idx <- grep(rx, lines)
    if (!length(keep_idx)) return(data.frame(parent = character(0), child = character(0), stringsAsFactors = FALSE))
    lines <- lines[keep_idx]
    ipr   <- regmatches(lines, regexpr(rx, lines))
    indent <- vapply(lines, function(s) {
      s2 <- gsub("\t", "  ", s, fixed = TRUE)
      nchar(s2) - nchar(sub("^\\s+", "", s2))
    }, integer(1))
    parent_stack <- character(0); indent_stack <- integer(0)
    edges_p <- character(0); edges_c <- character(0)
    for (i in seq_along(ipr)) {
      cur_ipr <- ipr[i]; cur_ind <- indent[i]
      while (length(indent_stack) && indent_stack[length(indent_stack)] >= cur_ind) {
        indent_stack <- indent_stack[-length(indent_stack)]
        parent_stack <- parent_stack[-length(parent_stack)]
      }
      if (length(parent_stack)) {
        edges_p <- c(edges_p, parent_stack[length(parent_stack)])
        edges_c <- c(edges_c, cur_ipr)
      }
      parent_stack <- c(parent_stack, cur_ipr)
      indent_stack <- c(indent_stack, cur_ind)
    }
    data.frame(parent = edges_p, child = edges_c, stringsAsFactors = FALSE)
  }
  .load_interpro_tree <- function(tree_source = c("auto","local","remote","none"),
                                  term_trees = NULL,
                                  tree_path = NULL,
                                  tree_url = NULL,
                                  interpro_release = NULL,
                                  timeout_s = 30) {
    tree_source <- match.arg(tree_source)
    if (tree_source == "none") {
      if (is.data.frame(term_trees) && ncol(term_trees) >= 2) {
        colnames(term_trees)[1:2] <- c("parent","child")
        return(list(df = term_trees, provenance = list(mode = "provided-df")))
      }
      return(list(df = NULL, provenance = list()))
    }
    if (is.character(term_trees) && length(term_trees) == 1L && file.exists(term_trees) &&
        tree_source %in% c("auto","local")) {
      ed <- .read_tsv_strict(term_trees); colnames(ed)[1:2] <- c("parent","child")
      return(list(df = ed, provenance = list(mode = "provided-tsv", path = normalizePath(term_trees, winslash = "/"))))
    }
    if (!is.null(tree_path) && file.exists(tree_path) && tree_source %in% c("auto","local")) {
      headl <- readLines(tree_path, n = 5, warn = FALSE)
      looks_tsv <- grepl("\t", paste(headl, collapse = "")) &&
        !inherits(try(read.table(text = paste(headl, collapse = "\n"), sep = "\t"), silent = TRUE), "try-error")
      if (looks_tsv) {
        ed <- .read_tsv_strict(tree_path); colnames(ed)[1:2] <- c("parent","child")
        return(list(df = ed, provenance = list(mode = "local-tsv", path = normalizePath(tree_path, winslash = "/"))))
      } else {
        ed <- .parse_parent_child_tree_txt(tree_path)
        return(list(df = ed, provenance = list(mode = "local-txt", path = normalizePath(tree_path, winslash = "/"))))
      }
    }
    if (tree_source == "local")
      stop("tree_source='local' but neither term_trees TSV nor tree_path was supplied.")
    if (!is.null(tree_url)) {
      tf <- tempfile(fileext = ".txt")
      if (.download_quiet(tree_url, tf, timeout_s)) {
        ed <- try(.read_tsv_strict(tf), silent = TRUE)
        if (!inherits(ed, "try-error") && ncol(ed) >= 2) {
          colnames(ed)[1:2] <- c("parent","child")
        } else {
          ed <- .parse_parent_child_tree_txt(tf)
        }
        return(list(df = ed, provenance = list(mode = "remote-url", url = tree_url)))
      }
    }
    url <- if (is.null(interpro_release)) .current_urls$parent_child else .release_urls(interpro_release)$parent_child
    tf <- tempfile(fileext = ".txt")
    if (.download_quiet(url, tf, timeout_s)) {
      ed <- .parse_parent_child_tree_txt(tf)
      return(list(df = ed, provenance = list(mode = if (is.null(interpro_release)) "remote-current" else "remote-release",
                                             url = url, release = interpro_release)))
    }
    if (tree_source == "auto") {
      bp <- system.file("extdata", "InterPro_ParentChildTreeFile.txt", package = "dndsR")
      if (nzchar(bp) && file.exists(bp)) {
        ed <- .parse_parent_child_tree_txt(bp)
        return(list(df = ed, provenance = list(mode = "bundled", path = normalizePath(bp, winslash = "/"))))
      }
    }
    warning("Failed to load InterPro parent-child tree; hierarchy features will be disabled.")
    list(df = NULL, provenance = list(error = "load-failed"))
  }
  .write_provenance <- function(out_dir, comp_name, prov_list) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) return(invisible(NULL))
    fn <- file.path(out_dir, sprintf("%s_interpro_provenance.json", comp_name))
    jsonlite::write_json(prov_list, path = fn, auto_unbox = TRUE, pretty = TRUE)
    invisible(fn)
  }

  .num <- function(x) suppressWarnings(as.numeric(x))

  .detect_current_release <- function(timeout_s = 30) {
    # Try to read the current release notes and parse a version like "100.0"
    tf <- tempfile(fileext = ".txt")
    if (.download_quiet(.current_urls$notes, tf, timeout_s)) {
      txt <- paste(readLines(tf, warn = FALSE), collapse = "\n")
      m <- regmatches(txt, regexpr("([0-9]{2,3}\\.[0-9])", txt))
      if (length(m) == 1L && nzchar(m)) return(m)
    }
    NA_character_
  }

  .seq_desc_releases <- function(start_release, step = 2.0, max_backtrack = 12L) {
    sr <- .num(start_release)
    if (!is.finite(sr)) return(character(0))
    stops <- pmax(sr - step * seq_len(max_backtrack), 0)
    as.character(sprintf("%.1f", stops))
  }

  .coverage_for_current <- function(iprs, timeout_s = 30) {
    tf <- tempfile(fileext = ".tsv")
    if (!.download_quiet(.current_urls$entry_list, tf, timeout_s)) return(NA_real_)
    el <- try(.read_interpro_entry_list(tf), silent = TRUE)
    if (inherits(el, "try-error") || !"ENTRY_AC" %in% names(el)) return(NA_real_)
    if (!length(iprs)) return(1.0)
    round(sum(iprs %in% el$ENTRY_AC) / length(iprs), 6)
  }

  .backtrack_release_until_full <- function(iprs,
                                            start_release,
                                            step = 2.0,
                                            max_backtrack = 12L,
                                            timeout_s = 30) {
    tried <- character(0); covv <- numeric(0)

    # 1) Evaluate the start (often current) first
    rels <- unique(c(start_release, .seq_desc_releases(start_release, step, max_backtrack)))
    for (rel in rels) {
      if (!nzchar(rel)) next
      cvr <- .coverage_for_release(iprs, rel, timeout_s)
      tried <- c(tried, rel); covv <- c(covv, cvr)
      if (isTRUE(cvr == 1.0)) {
        return(list(release = rel, coverage = cvr, tried = tried, coverages = covv))
      }
    }

    # 2) No 100% hit; choose the best non-NA coverage among tried
    if (!length(covv)) return(list(release = NULL, coverage = NA_real_, tried = tried, coverages = covv))
    idx <- which.max(replace(covv, is.na(covv), -Inf))
    list(release = tried[idx], coverage = covv[idx], tried = tried, coverages = covv)
  }

  # -------------------------------------------------------------------------------

  # ---------- generic helpers (existing) ----------
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
  
      if (inherits(ok, "try-error")) {
        stop("filter_expr evaluation failed: ", as.character(ok))
      }
  
      ok <- as.vector(ok)
      if (length(ok) != nrow(d)) {
        stop("filter_expr must evaluate to length nrow(d).")
      }
  
      ok <- as.logical(ok)
      keep <- keep & !is.na(ok) & ok
    }
  
    d[keep, , drop = FALSE]
  }
  .split_terms_unique <- function(x) {
    if (is.null(x)) return(character(0))
    s <- as.character(x)[1]
    if (is.na(s) || !nzchar(s)) return(character(0))
    unique(strsplit(s, term_sep, fixed = TRUE)[[1]])
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
    } else if (method == "BY") {
      res$p_adj <- stats::p.adjust(res$p_value, method = "BY")
    } else if (method == "qvalue") {
      if (requireNamespace("qvalue", quietly = TRUE)) {
        res$p_adj <- as.numeric(qvalue::qvalue(res$p_value)$qvalues)
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
    } else {
      res$p_adj <- res$p_value
    }
    res
  }

  # ---------- logging helper (MATCH go_enrichment behavior) ----------
  .with_log <- function(log_file, tag, header = NULL, expr) {
    # If dndsR internal logger exists, use it. Otherwise just run.
    if (exists(".dndsr_with_log", mode = "function", inherits = TRUE)) {
      .dndsr_with_log(
        log_file = log_file,
        tag = tag,
        header = header,
        expr = expr
      )
    } else {
      force(expr)
    }
  }                                         

  # ---------- font & plotting helpers (MATCH go_enrichment) ----------
  .bundled_arial_path <- function() {
    ttf <- system.file("fonts", "ArialBold.ttf", package = "dndsR")
    if (!nzchar(ttf) || !file.exists(ttf)) return(NULL)
    ttf
  }

  .setup_showtext <- function() {
    ttf <- .bundled_arial_path()
    if (is.null(ttf)) return(NULL)

    if (!requireNamespace("showtext", quietly = TRUE)) return(NULL)
    if (!requireNamespace("sysfonts", quietly = TRUE)) return(NULL)

    ok <- try({
      sysfonts::font_add(family = "Arial Bold", regular = ttf)
      showtext::showtext_auto(TRUE)
      TRUE
    }, silent = TRUE)

    if (isTRUE(ok)) "Arial Bold" else NULL
  }

  .register_svglite_mapping <- function() {
    if (!requireNamespace("svglite", quietly = TRUE)) return(NULL)
    ttf <- .bundled_arial_path()
    if (is.null(ttf)) return(NULL)

    function(file, ...) svglite::svglite(
      file,
      user_fonts = list(`Arial Bold` = ttf),
      ...
    )
  }

  .pick_sans_family <- function() {
    fam <- .setup_showtext()
    if (!is.null(fam)) return(fam)
    "sans"
  }

  .svg_device <- function() {
    dev_map <- .register_svglite_mapping()
    if (!is.null(dev_map)) return(dev_map)
    if (requireNamespace("svglite", quietly = TRUE)) return(function(file, ...) svglite::svglite(file, ...))
    NULL
  }

  .upper_padj <- function(d, alpha) {
    max_p <- suppressWarnings(max(d$p_adj, na.rm = TRUE))
    if (!is.finite(max_p)) max_p <- alpha
    eps <- max(1e-12, alpha * 1e-6)
    max(max_p, alpha + eps)
  }
  .padj_scale <- function(alpha, upper) {
    ggplot2::scale_color_gradientn(
      colours = c("red", "grey80", "steelblue"),
      values  = scales::rescale(c(0, alpha, upper), to = c(0,1), from = c(0, upper)),
      limits  = c(0, upper),
      oob     = scales::squish, name = "adj p"
    )
  }
                                                            
  .write_plot <- function(df, ylab, out_svg, alpha_val = alpha,
                          x_axis_min = x_axis_min,
                          x_axis_max = x_axis_max,
                          x_axis_pad_right = x_axis_pad_right) {
  
    if (!make_plots || !requireNamespace("ggplot2", quietly = TRUE) || !nrow(df)) {
      return(invisible(NULL))
    }
  
    top <- df[order(df$p_adj, -df$enrichment, df$IPR), , drop = FALSE]
    top <- utils::head(top, top_n)
  
    yvar <- if ("label" %in% names(top)) "label" else "IPR"
    top$y_lab <- top[[yvar]]
  
    # ---- Identify non-finite enrichments (Inf/NA) and create plot-safe x ----
    top$is_inf_enrichment <- !is.finite(top$enrichment) | is.na(top$enrichment)
    top$enrichment_plot <- top$enrichment
  
    finite_x <- top$enrichment_plot[is.finite(top$enrichment_plot)]
    max_finite <- if (length(finite_x)) max(finite_x, na.rm = TRUE) else 1
  
    # cap used only when enrichment is non-finite/NA so it can be drawn
    cap_x <- max_finite * 1.05
    if (!is.finite(cap_x) || cap_x <= 0) cap_x <- 1
  
    if (any(top$is_inf_enrichment)) top$enrichment_plot[top$is_inf_enrichment] <- cap_x
  
    # ---- Make it explicit in the label ----
    #top$y_lab_plot <- top$y_lab
    #if (any(top$is_inf_enrichment)) {
    #  top$y_lab_plot[top$is_inf_enrichment] <- paste0(top$y_lab_plot[top$is_inf_enrichment], " (Inf)")
    #}
  
    # ---- Drop only rows that truly can't be mapped (pos_count/p_adj) ----
    keep <- !is.na(top$pos_count) & is.finite(top$pos_count) &
            !is.na(top$p_adj)     & is.finite(top$p_adj)
    top_plot <- top[keep, , drop = FALSE]
    if (!nrow(top_plot)) return(invisible(NULL))
  
    base_family <- .pick_sans_family()
    upper <- .upper_padj(top_plot, alpha_val)
  
    # ---- Base plot: map Inf-flag to shape + stroke, keep color = p_adj ----
    gg <- ggplot2::ggplot(
      top_plot,
      ggplot2::aes(
        x = enrichment_plot,
        y = stats::reorder(y_lab, -p_adj),
        size = pos_count,
        color = p_adj,
        shape = is_inf_enrichment
      )
    ) +
      ggplot2::geom_point(stroke = 1.2) +
      .padj_scale(alpha_val, upper) +
      ggplot2::scale_shape_manual(
        values = c(`FALSE` = 16, `TRUE` = 1),
        labels = c(`FALSE` = "finite", `TRUE` = "inf"),
        name = "enrichment"
      ) +
      ggplot2::scale_size_continuous(name = "# pos") +
      ggplot2::labs(x = "Enrichment (pos/bg)", y = ylab) +
      ggplot2::theme_minimal(base_size = 13, base_family = base_family) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.05, add = 0)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(plot.margin = ggplot2::margin(5.5, 50, 5.5, 5.5))
  
    # ---- Optional fixed x limits (still no hard clipping) ----
    if (!is.null(x_axis_min) || !is.null(x_axis_max)) {
      xmin <- if (is.null(x_axis_min)) -Inf else x_axis_min
      xmax <- if (is.null(x_axis_max))  Inf else x_axis_max
      gg <- gg + ggplot2::coord_cartesian(xlim = c(xmin, xmax), clip = "off")
    }
  
    dev_fun <- .svg_device()
    if (!is.null(dev_fun)) {
      ggplot2::ggsave(out_svg, gg, device = dev_fun, width = 11, height = 9)
    } else {
      ggplot2::ggsave(out_svg, gg, width = 11, height = 9)
    }
  
    invisible(NULL)
  }
                                                            
  # ---------- InterPro entry.list (existing minimal) ----------
  .read_entries <- function(path) .read_tsv_strict(path)
  .bundled_path <- function() {
    p <- system.file("extdata", "interpro_entry.list.tsv", package = "dndsR")
    if (!nzchar(p) || !file.exists(p)) return(NULL)
    p
  }
  .build_ipr_maps <- function(entries) {
    if (is.null(entries)) return(list(type_by_ipr = NULL, name_by_ipr = NULL))
    ac <- as.character(entries$ENTRY_AC)
    ty <- as.character(entries$ENTRY_TYPE)
    nm <- as.character(entries$ENTRY_NAME)
    list(type_by_ipr = stats::setNames(ty, ac),
         name_by_ipr = stats::setNames(nm, ac))
  }

  # ---------- tree/metadata helpers ----------
  .load_df <- function(obj) {
    if (is.null(obj)) return(NULL)
    if (is.character(obj) && length(obj) == 1L && file.exists(obj)) {
      return(utils::read.table(obj, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE))
    }
    if (is.data.frame(obj)) return(obj)
    NULL
  }
  .expand_descendants <- function(seeds, tree_df, depth = Inf, limit = Inf) {
    if (is.null(tree_df) || !length(seeds)) return(seeds)
    parents <- as.character(tree_df[[1]])
    children <- as.character(tree_df[[2]])
    adj <- split(children, parents)
    seen <- unique(seeds); frontier <- unique(seeds); curd <- 0L
    while (length(frontier) && curd < depth) {
      kids <- unique(unlist(adj[frontier], use.names = FALSE))
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
  .build_exclude_set <- function(excl) {
    if (is.null(excl) || !length(excl)) return(NULL)
    x <- unique(as.character(excl))
    stats::setNames(rep(TRUE, length(x)), x)
  }
  .filter_terms_string <- function(s, allowed_types = NULL, type_by_ipr = NULL,
                                   exclude_set = NULL, keep_unmatched = TRUE) {
    if (is.na(s) || !nzchar(s)) return(s)
    parts <- unique(strsplit(s, term_sep, fixed = TRUE)[[1]])
    if (!length(parts)) return("")
    keep <- logical(length(parts))
    for (i in seq_along(parts)) {
      ipr <- parts[i]
      if (!is.null(exclude_set) && isTRUE(exclude_set[ipr])) { keep[i] <- FALSE; next }
      if (!is.null(allowed_types)) {
        ty <- if (!is.null(type_by_ipr)) type_by_ipr[ipr] else NA_character_
        keep[i] <- if (is.na(ty)) isTRUE(keep_unmatched) else ty %in% allowed_types
      } else keep[i] <- TRUE
    }
    parts2 <- parts[keep]
    if (!length(parts2)) return("")
    paste(parts2, collapse = term_sep)
  }
  .sanitize_type <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    gsub("_+", "_", gsub("^_|_$", "", x))
  }

  # ---------- core fisher ----------
  .fisher_from_vectors <- function(vec_all, vec_pos, drop_empty) {
    vec_all <- as.character(vec_all)
    vec_pos <- as.character(vec_pos)
    if (drop_empty) {
      vec_all <- vec_all[!is.na(vec_all) & nzchar(vec_all)]
      vec_pos <- vec_pos[!is.na(vec_pos) & nzchar(vec_pos)]
    }
    n_pos <- length(vec_pos); n_all <- length(vec_all); n_bg <- n_all - n_pos
    if (n_pos == 0L || n_all == 0L || n_bg < 0L) return(NULL)
    all_tab <- .count_terms_from_vec(vec_all)
    pos_tab <- .count_terms_from_vec(vec_pos)
    if (length(all_tab) == 0) return(NULL)
    res <- do.call(rbind, lapply(names(all_tab), function(tt) {
      a <- as.integer(ifelse(tt %in% names(pos_tab), pos_tab[[tt]], 0L))
      c <- as.integer(all_tab[[tt]] - a)
      b <- n_pos - a
      d <- n_bg  - c
      ft <- try(stats::fisher.test(matrix(c(a,b,c,d), nrow = 2), alternative = "greater"), silent = TRUE)
      ci <- c(NA_real_, NA_real_); or <- NA_real_; p <- NA_real_
      if (!inherits(ft, "try-error")) {
        p  <- ft$p.value
        or <- as.numeric(ft$estimate)
        if (!is.null(ft$conf.int)) ci <- as.numeric(ft$conf.int)[1:2]
      }
      data.frame(IPR = tt, pos_count = a, nonpos_count = c,
                 pos_total = n_pos, nonpos_total = n_bg,
                 odds_ratio = or, ci_lower = ci[1], ci_upper = ci[2],
                 p_value = p, stringsAsFactors = FALSE)
    }))
    if (is.null(res) || !nrow(res)) return(NULL)
    res$total_count <- res$pos_count + res$nonpos_count
    pc <- pmax(res$pos_count, 0); pt <- pmax(res$pos_total, 1)
    nc <- pmax(res$nonpos_count, 0); nt <- pmax(res$nonpos_total, 1)
    bg_rate  <- nc / nt; pos_rate <- pc / pt
    res$enrichment <- ifelse(bg_rate == 0 & pos_rate > 0, Inf,
                             ifelse(bg_rate == 0 & pos_rate == 0, 1, pos_rate / bg_rate))
    res
  }

  # ---------- metadata attach ----------
  .attach_metadata_to_res <- function(res, type_by_ipr, name_by_ipr) {
    if (is.null(res) || !nrow(res)) return(res)
    if (is.null(type_by_ipr) && is.null(name_by_ipr)) {
      res$ENTRY_TYPE <- NA_character_; res$ENTRY_NAME <- NA_character_
    } else {
      ty <- type_by_ipr[res$IPR]; nm <- name_by_ipr[res$IPR]
      res$ENTRY_TYPE <- unname(ifelse(is.na(ty), NA_character_, ty))
      res$ENTRY_NAME <- unname(ifelse(is.na(nm), NA_character_, nm))
    }
    res$label <- ifelse(is.na(res$ENTRY_NAME) | !nzchar(res$ENTRY_NAME),
                        res$IPR, paste0(res$IPR, " ", res$ENTRY_NAME))
    res
  }

  # ---------- hierarchy ----------
  .build_term2rows <- function(df, term_col) {
    ids <- seq_len(nrow(df))
    id_col <- if (grepl("^q_", term_col) && "q_id" %in% names(df)) "q_id" else
      if (grepl("^s_", term_col) && "s_id" %in% names(df)) "s_id" else NULL
    if (!is.null(id_col)) ids <- df[[id_col]]
    term2rows <- list()
    for (i in seq_len(nrow(df))) {
      ts <- .split_terms_unique(df[[term_col]][i])
      if (!length(ts)) next
      for (t in ts) term2rows[[t]] <- c(term2rows[[t]], ids[i])
    }
    lapply(term2rows, unique)
  }
  .filter_terms_by_thresholds <- function(term2rows, pos_ids, all_ids, min_total, min_pos, max_prop) {
    keep <- list(); n_all <- length(all_ids)
    for (t in names(term2rows)) {
      members <- intersect(term2rows[[t]], all_ids)
      tot <- length(members)
      pos <- length(intersect(members, pos_ids))
      prop <- if (n_all > 0) tot / n_all else 0
      if (tot >= min_total && pos >= min_pos && prop <= max_prop) keep[[t]] <- members
    }
    keep
  }
  .dag_from_term_trees <- function(term2rows, type_terms, tree_df) {
    if (is.null(tree_df) || ncol(tree_df) < 2) return(NULL)
    parent <- as.character(tree_df[[1]]); child <- as.character(tree_df[[2]])
    parent <- parent[parent %in% type_terms]; child <- child[child %in% type_terms]
    if (!length(parent)) return(NULL)
    split(child, parent)
  }
  .dag_infer_by_subset <- function(term2rows, type_terms) {
    if (!length(type_terms)) return(NULL)
    children <- lapply(type_terms, function(x) character(0)); names(children) <- type_terms
    for (i in seq_along(type_terms)) {
      ti <- type_terms[i]; si <- term2rows[[ti]]
      if (!length(si)) next
      for (j in seq_along(type_terms)) {
        if (i == j) next
        tj <- type_terms[j]; sj <- term2rows[[tj]]
        if (!length(sj)) next
        if (length(sj) <= length(si) && all(sj %in% si) && !all(si %in% sj))
          children[[ti]] <- c(children[[ti]], tj)
      }
    }
    lapply(children, unique)
  }
  .topo_specific_first <- function(children_map) {
    nodes <- unique(c(names(children_map), unlist(children_map, use.names = FALSE)))
    depth <- setNames(integer(length(nodes)), nodes)
    visited <- setNames(logical(length(nodes)), nodes)
    rec <- function(n) {
      if (visited[[n]]) return(depth[[n]])
      visited[[n]] <<- TRUE
      ch <- children_map[[n]]
      if (is.null(ch) || !length(ch)) { depth[[n]] <<- 0L; return(0L) }
      d <- max(vapply(ch, rec, integer(1), USE.NAMES = FALSE)) + 1L
      depth[[n]] <<- d; d
    }
    for (n in nodes) rec(n)
    nodes[order(depth[nodes])]  # leaves (most specific) first
  }
  .parent_child_enrich <- function(df, term_col, type_label, type_by_ipr, name_by_ipr,
                                   min_total, min_pos, max_prop,
                                   fdr_method, alpha,
                                   term_trees_df = NULL) {
    n_all <- nrow(df); if (!n_all) return(NULL)
    pos_id <- which(df$dNdS > pos_threshold); all_id <- seq_len(n_all)
    term2rows_all <- .build_term2rows(df, term_col); if (!length(term2rows_all)) return(NULL)
    term2rows <- .filter_terms_by_thresholds(term2rows_all, pos_id, all_id, min_total, min_pos, max_prop)
    if (!length(term2rows)) return(NULL)
    type_terms <- names(term2rows)

    children_map <- .dag_from_term_trees(term2rows, type_terms, term_trees_df)
    if (is.null(children_map)) children_map <- .dag_infer_by_subset(term2rows, type_terms)
    if (is.null(children_map)) {
      vec_all <- vapply(df[[term_col]], function(s) {
        parts <- .split_terms_unique(s); parts <- parts[parts %in% type_terms]
        if (!length(parts)) "" else paste(unique(parts), collapse = term_sep)
      }, FUN.VALUE = character(1))
      vec_pos <- vec_all[pos_id]
      res <- .fisher_from_vectors(vec_all, vec_pos, drop_rows_without_term)
      if (is.null(res) || !nrow(res)) return(NULL)
      res <- res[res$IPR %in% type_terms, , drop = FALSE]
      res <- .attach_metadata_to_res(res, type_by_ipr, name_by_ipr); res$ENTRY_TYPE <- type_label
      res <- .adjust_pvals(res, fdr_method, alpha)
      return(res[order(res$p_adj, -res$enrichment, res$IPR), , drop = FALSE])
    }

    order_terms <- .topo_specific_first(children_map)
    claimed <- setNames(vector("list", length(order_terms)), order_terms)
    for (k in seq_along(claimed)) claimed[[k]] <- integer(0)
    rows <- list(); n_pos <- length(pos_id); n_bg <- n_all - n_pos
    sig <- setNames(logical(length(order_terms)), order_terms)

    for (t in order_terms) {
      members <- intersect(term2rows[[t]], all_id)
      if (!length(members)) next
      kids <- children_map[[t]]
      if (!is.null(kids) && length(kids)) {
        kids_sig <- kids[kids %in% names(sig)[sig]]
        if (length(kids_sig)) {
          kids_claimed <- unique(unlist(claimed[kids_sig], use.names = FALSE))
          members <- setdiff(members, kids_claimed)
        }
      }
      tot <- length(members); if (tot < min_total) next
      a <- length(intersect(members, pos_id)); c <- tot - a
      b <- n_pos - a; d <- n_bg - c
      if (a < min_pos) next
      ft <- try(stats::fisher.test(matrix(c(a,b,c,d), nrow = 2), alternative = "greater"), silent = TRUE)
      ci <- c(NA_real_, NA_real_); or <- NA_real_; p <- NA_real_
      if (!inherits(ft, "try-error")) {
        p  <- ft$p.value; or <- as.numeric(ft$estimate)
        if (!is.null(ft$conf.int)) ci <- as.numeric(ft$conf.int)[1:2]
      }
      rows[[t]] <- data.frame(IPR = t, pos_count = a, nonpos_count = c,
                              pos_total = n_pos, nonpos_total = n_bg,
                              odds_ratio = or, ci_lower = ci[1], ci_upper = ci[2],
                              p_value = p, total_count = tot, stringsAsFactors = FALSE)
      if (is.finite(p) && !is.na(p) && p <= alpha) sig[t] <- TRUE
    }

    nodes <- unique(c(names(children_map), unlist(children_map, use.names = FALSE)))
    parents_map <- setNames(vector("list", length(nodes)), nodes)
    for (p in names(children_map)) for (c in children_map[[p]]) parents_map[[c]] <- c(parents_map[[c]], p)
    claimed <- setNames(vector("list", length(nodes)), nodes)
    for (k in seq_along(claimed)) claimed[[k]] <- integer(0)
    for (t in order_terms) {
      r <- rows[[t]]; if (is.null(r)) next
      if (!isTRUE(sig[[t]])) next
      members <- term2rows[[t]]; an <- parents_map[[t]]
      if (length(an)) for (a_name in an) claimed[[a_name]] <- unique(c(claimed[[a_name]], members))
    }

    if (!length(rows)) return(NULL)
    res <- do.call(rbind, unname(rows))
    pc <- pmax(res$pos_count, 0); pt <- pmax(res$pos_total, 1)
    nc <- pmax(res$nonpos_count, 0); nt <- pmax(res$nonpos_total, 1)
    bg_rate  <- nc / nt; pos_rate <- pc / pt
    res$enrichment <- ifelse(bg_rate == 0 & pos_rate > 0, Inf,
                             ifelse(bg_rate == 0 & pos_rate == 0, 1, pos_rate / bg_rate))
    res <- .attach_metadata_to_res(res, type_by_ipr, name_by_ipr); res$ENTRY_TYPE <- type_label
    res <- .adjust_pvals(res, fdr_method, alpha)

    if (!is.null(children_map) && length(children_map)) {
      novelty_ok <- rep(TRUE, nrow(res)); names(novelty_ok) <- res$IPR
      for (t in res$IPR) {
        kids <- children_map[[t]]; if (!length(kids)) next
        kids_sig <- intersect(kids, res$IPR[res$p_adj <= alpha]); if (!length(kids_sig)) next
        members <- term2rows[[t]]; kids_union <- unique(unlist(term2rows[kids_sig], use.names = FALSE))
        frac <- length(setdiff(members, kids_union)) / max(length(members), 1)
        if (frac < ancestor_novel_frac) novelty_ok[[t]] <- FALSE
      }
      res <- res[novelty_ok[res$IPR], , drop = FALSE]
    }
    res[order(res$p_adj, -res$enrichment, res$IPR), , drop = FALSE]
  }

  # ---------- Decide release & load InterPro resources ----------
  .collect_all_iprs_from_comp <- function(path, term_sep = ";") {
    if (!file.exists(path)) return(character(0))
    d <- utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    .collect_iprs_from_df(d, term_sep = term_sep)
  }

  chosen_release <- interpro_release
  best_cov <- NA_real_
  backtrack_meta <- list()

  if (isTRUE(auto_detect_release) && is.null(interpro_release)) {
    # Collect IPRs from input(s)
    iprs_all <- character(0)
    if (!is.null(comparison_file)) {
      df_cf <- .read_comparisons(comparison_file)
      for (i in seq_len(nrow(df_cf))) {
        comp <- df_cf$comparison_name[i]
        in_file <- file.path(output_dir, comp, paste0(comp, "_dnds_annot.tsv"))
        if (file.exists(in_file)) iprs_all <- unique(c(iprs_all, .collect_all_iprs_from_comp(in_file, term_sep)))
      }
    } else if (!is.null(dnds_annot_file)) {
      iprs_all <- unique(c(iprs_all, .collect_all_iprs_from_comp(dnds_annot_file, term_sep)))
    }

    # 1) Get "current" release number if possible; else fall back to highest candidate or a sane default
    cur_rel <- .detect_current_release(timeout_s = max(entries_timeout_s, tree_timeout_s))
    if (!nzchar(cur_rel) || !is.finite(.num(cur_rel))) {
      stop("[ipr_enrichment] Auto-detection of current InterPro release failed. ",
           "Please specify 'interpro_release' manually or provide local copies of ",
           "'entry.list' and 'ParentChildTreeFile.txt'.")
    }

    # 2) Evaluate "current" coverage first (even if backtracking later)
    cur_cov <- .coverage_for_current(iprs_all, timeout_s = max(entries_timeout_s, tree_timeout_s))

    # 3) Walk downwards until 100% coverage, else keep best
    step_sz <- 2.0
    max_back <- 12L
    bt <- .backtrack_release_until_full(iprs_all, cur_rel, step = step_sz, max_backtrack = max_back,
                                        timeout_s = max(entries_timeout_s, tree_timeout_s))

    chosen_release <- bt$release
    best_cov       <- bt$coverage
    backtrack_meta <- list(
      current_release_detected = cur_rel,
      current_coverage_estimate = cur_cov,
      attempted_releases = bt$tried,
      attempted_coverages = bt$coverages
    )

    if (is.na(best_cov)) {
      message("[ipr_enrichment] Could not estimate coverage across releases; using current_release endpoints.")
      chosen_release <- NULL
    } else {
      message(sprintf("[ipr_enrichment] Selected InterPro release %s (coverage %.3f)", chosen_release, best_cov))
      if (is.finite(strict_coverage) && best_cov < strict_coverage) {
        stop(sprintf("Coverage %.3f < strict_coverage=%.2f for selected release %s.",
                     best_cov, strict_coverage, chosen_release))
      }
    }
  }


  el <- .load_interpro_entry_list(entries_source = entries_source,
                                  entries_path   = entries_path,
                                  entries_url    = if (!is.null(entries_url)) entries_url else NULL,
                                  interpro_release = chosen_release,
                                  timeout_s = entries_timeout_s)
  entries <- el$df
  maps <- .build_ipr_maps(entries)
  type_by_ipr <- maps$type_by_ipr
  name_by_ipr <- maps$name_by_ipr

  tr <- .load_interpro_tree(tree_source = tree_source,
                            term_trees  = term_trees,
                            tree_path   = tree_path,
                            tree_url    = if (!is.null(tree_url)) tree_url else NULL,
                            interpro_release = chosen_release,
                            timeout_s = tree_timeout_s)
  tree_df <- tr$df

  .prov_common <- list(
    chosen_release = chosen_release,
    auto_detect_release = auto_detect_release,
    best_coverage = best_cov,
    backtrack = backtrack_meta,
    entries_provenance = el$provenance,
    tree_provenance = tr$provenance
  )

  # ---------- exclusions (now using loaded tree_df if asked) ----------
  exclude_seeds <- exclude_ids
  if (isTRUE(exclude_descendants) && length(exclude_seeds)) {
    if (!is.null(tree_df) && ncol(tree_df) >= 2) {
      exclude_seeds <- .expand_descendants(exclude_seeds, tree_df,
                                           depth = exclude_descendants_depth,
                                           limit = exclude_descendants_limit)
    } else {
      warning("exclude_descendants requested but no parent-child tree available; skipping expansion.")
    }
  }
  exclude_set <- .build_exclude_set(exclude_seeds)

  # ---------- pooled helper ----------
  .pooled_enrichment <- function(df, side, allowed_types = NULL) {
    prefix <- if (side == "query") "q_" else "s_"
    col <- paste0(prefix, "ipr")
    if (!col %in% names(df) || !is.character(df[[col]])) return(NULL)
    vec_all <- vapply(df[[col]], .filter_terms_string, FUN.VALUE = character(1),
                      allowed_types = allowed_types, type_by_ipr = type_by_ipr,
                      exclude_set = exclude_set, keep_unmatched = keep_unmatched)
    pos <- df$dNdS > pos_threshold
    vec_pos <- vec_all[pos]
    res <- .fisher_from_vectors(vec_all, vec_pos, drop_rows_without_term)
    if (is.null(res) || !nrow(res)) return(NULL)
    n_all <- length(vec_all)
    prop  <- res$total_count / max(n_all, 1)
    keep  <- (res$total_count >= min_total) & (res$pos_count >= min_pos) & (prop <= max_prop)
    res <- res[keep, , drop = FALSE]
    if (!nrow(res)) return(NULL)
    res <- .attach_metadata_to_res(res, type_by_ipr, name_by_ipr)
    res <- .adjust_pvals(res, fdr_method, alpha)
    res
  }

  # ---------- per-side/type runner ----------
  .enrich_side_ipr <- function(d, side, comp, comp_dir, stratified_types = NULL) {
    prefix <- if (side == "query") "q_" else "s_"
    term_col <- paste0(prefix, "ipr")
    if (!term_col %in% names(d) || !is.character(d[[term_col]])) return(list())
    df <- .apply_filter(d, filter_expr); if (!nrow(df)) return(list())
    filter_fun <- function(v, keep_types = NULL) vapply(v, .filter_terms_string, FUN.VALUE = character(1),
                                                        allowed_types = keep_types, type_by_ipr = type_by_ipr,
                                                        exclude_set = exclude_set, keep_unmatched = keep_unmatched)
    results <- list()
    if (!is.null(stratified_types)) {
      for (tp in stratified_types) {
        df2 <- df; df2[[term_col]] <- filter_fun(df[[term_col]], keep_types = tp)
        if (method == "fisher") {
          res <- .pooled_enrichment(df2, side, allowed_types = tp)
          if (!is.null(res) && nrow(res)) { res$ENTRY_TYPE <- tp; results[[tp]] <- res }
          else message(sprintf("[ipr_enrichment] No rows for type=%s side=%s in %s", tp, side, comp))
        } else {
          res <- .parent_child_enrich(df2, term_col, tp, type_by_ipr, name_by_ipr,
                                      min_total, min_pos, max_prop,
                                      fdr_method, alpha,
                                      term_trees_df = tree_df)
          if (!is.null(res) && nrow(res)) results[[tp]] <- res
          else message(sprintf("[ipr_enrichment] No rows (hierarchy mode) for type=%s side=%s in %s", tp, side, comp))
        }
      }
      return(results)
    }
    df[[term_col]] <- filter_fun(df[[term_col]], keep_types = include_types)
    res <- .pooled_enrichment(df, side, allowed_types = include_types)
    if (is.null(res) || !nrow(res)) { message(sprintf("[ipr_enrichment] No rows for pooled side=%s", side)); return(list()) }
    list(pooled = res)
  }

  # ---------- comparison runner ----------
  .run_one_comparison <- function(comp_name, comp_dir) {
    in_file <- file.path(comp_dir, paste0(comp_name, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) stop("Annotated dN/dS file not found: ", in_file)
    d <- utils::read.table(in_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

    outs <- character(0)
    strat_types <- NULL
    if (isTRUE(stratify_by_type)) {
      if (is.null(type_by_ipr)) stop("stratify_by_type=TRUE requires InterPro metadata (entries_source != 'none').")
      if (is.null(types)) {
        all_terms <- character(0)
        if ("q_ipr" %in% names(d) && is.character(d$q_ipr))
          all_terms <- c(all_terms, unlist(strsplit(paste(na.omit(d$q_ipr)), term_sep, fixed = TRUE)))
        if ("s_ipr" %in% names(d) && is.character(d$s_ipr))
          all_terms <- c(all_terms, unlist(strsplit(paste(na.omit(d$s_ipr)), term_sep, fixed = TRUE)))
        all_terms <- unique(all_terms[nzchar(all_terms)])
        if (!is.null(exclude_ids) && length(exclude_ids)) all_terms <- setdiff(all_terms, exclude_ids)
        strat_types <- sort(unique(type_by_ipr[all_terms])); strat_types <- strat_types[!is.na(strat_types)]
      } else {
        strat_types <- intersect(types, unique(type_by_ipr))
      }
      if (!length(strat_types)) {
        message("No InterPro ENTRY_TYPEs detected in ", in_file, "; skipping stratified analysis.")
        return(character(0))
      }
    }

    for (sd in sides) {
      slices <- .enrich_side_ipr(d, side = sd, comp = comp_name, comp_dir = comp_dir,
                                 stratified_types = if (isTRUE(stratify_by_type)) strat_types else NULL)
      if (!length(slices)) next
      if (isTRUE(stratify_by_type)) {
        allres <- do.call(rbind, Filter(function(x) is.data.frame(x) && nrow(x) > 0, slices))
        if (!nrow(allres)) next
        allres$side <- sd; allres$comparison <- comp_name
        for (tp in unique(allres$ENTRY_TYPE)) {
          sub <- allres[allres$ENTRY_TYPE == tp, , drop = FALSE]
          if (!nrow(sub)) next
          sub <- sub[order(sub$p_adj, -sub$enrichment, sub$IPR), , drop = FALSE]
          out_tsv <- file.path(comp_dir, sprintf("%s_%s_IPR_%s_enrichment.tsv",
                                                 comp_name, if (sd=="query") "q" else "s", .sanitize_type(tp)))
          utils::write.table(sub, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
          out_svg <- sub("_enrichment.tsv$", "_enrichment_top20.svg", out_tsv)
          .write_plot(sub, sprintf("IPR (%s, %s) [%s]", tp, sd, method),
                      out_svg, x_axis_min = x_axis_min, x_axis_max = x_axis_max, x_axis_pad_right = x_axis_pad_right)
          outs <- c(outs, out_tsv)
        }
      } else {
        res <- slices[[1]]; if (is.null(res) || !nrow(res)) next
        res$side <- sd; res$comparison <- comp_name
        res <- res[order(res$p_adj, -res$enrichment, res$IPR), , drop = FALSE]
        out_tsv <- file.path(comp_dir, sprintf("%s_%s_IPR_enrichment.tsv", comp_name, if (sd=="query") "q" else "s"))
        utils::write.table(res, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
        out_svg <- sub("_enrichment.tsv$", "_enrichment_top20.svg", out_tsv)
        .write_plot(res, sprintf("IPR (%s) [%s]", sd, method),
                    out_svg, x_axis_min = x_axis_min, x_axis_max = x_axis_max, x_axis_pad_right = x_axis_pad_right)
        outs <- c(outs, out_tsv)
      }
    }

    if (!length(outs)) message("No IPR enrichments produced for ", comp_name)

    # NEW: provenance sidecar
    prov <- c(.prov_common, list(comparison = comp_name,
                                 timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
    .write_provenance(comp_dir, comp_name, prov)

    outs
  }

  # ---------- batch vs single ----------
  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)
    outs <- character(0)

    for (i in seq_len(nrow(df))) {
      comp <- df$comparison_name[i]
      comp_dir <- file.path(output_dir, comp)
      dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

      message(sprintf("[ipr_enrichment] %s (logging to %s/%s_ipr_enrichment.log)",
                      comp, comp_dir, comp))

      log_file <- file.path(comp_dir, sprintf("%s_ipr_enrichment.log", comp))

      out_i <- .with_log(
        log_file = log_file,
        tag = "ipr_enrichment",
        header = c(
          sprintf("[ipr_enrichment] comp=%s", comp),
          sprintf("[ipr_enrichment] pid=%d", Sys.getpid())
        ),
        expr = {
          .run_one_comparison(comp, comp_dir)
        }
      )

      outs <- c(outs, out_i)
    }

    message("All IPR enrichments complete.")
    return(invisible(outs))
  }

  if (is.null(dnds_annot_file)) stop("Provide either comparison_file (batch) OR dnds_annot_file (single).")
  if (!file.exists(dnds_annot_file)) stop("dnds_annot_file not found: ", dnds_annot_file)
  d <- utils::read.table(dnds_annot_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  comp_dir  <- dirname(dnds_annot_file)
  comp_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file))

  message(sprintf("[ipr_enrichment] %s (logging to %s/%s_ipr_enrichment.log)",
                  comp_name, comp_dir, comp_name))

  log_file <- file.path(comp_dir, sprintf("%s_ipr_enrichment.log", comp_name))

  outs <- .with_log(
    log_file = log_file,
    tag = "ipr_enrichment",
    header = c(sprintf("[ipr_enrichment] comp=%s", comp_name)),
    expr = {
      .run_one_comparison(comp_name, comp_dir)
    }
  )

  message("IPR enrichment complete for: ", dnds_annot_file)
  invisible(outs)
}
