# dnds_state_contrast.R
# dN/dS state contrasts via logistic regression (optional mixed model) + forest plot
#
# Merged "pairwise" + "gene" modes:
#   - level = "pairwise": model states on per-ortholog-pair rows
#   - level = "gene":     stack *_dnds_annot.tsv, convert to gene rows (query/subject),
#                         optionally filter by regions, aggregate across partners per gene,
#                         then model states on gene-level observations
#
# NOTE: Under development; no backwards compatibility guaranteed.

# -----------------------------
# Internal helpers
# -----------------------------

#' Read a whitespace-delimited file with optional header
#'
#' Internal helper: tries header=TRUE first, then header=FALSE.
#'
#' @keywords internal
.read_ws <- function(path, header_try = TRUE) {
  utils::read.table(path,
                    header           = header_try,
                    sep              = "",
                    quote            = "\"",
                    stringsAsFactors = FALSE,
                    comment.char     = "",
                    strip.white      = TRUE,
                    blank.lines.skip = TRUE,
                    check.names      = FALSE)
}

#' Read comparison file used throughout dndsR
#'
#' Expects columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff
#'
#' @keywords internal
.read_comparisons <- function(x) {
  req <- c("comparison_name", "query_fasta", "query_gff", "subject_fasta", "subject_gff")
  if (is.data.frame(x)) {
    stopifnot(all(req %in% names(x)))
    return(x[, req, drop = FALSE])
  }
  df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
  if (!inherits(df1, "try-error") && all(req %in% names(df1))) {
    return(df1[, req, drop = FALSE])
  }
  df2 <- .read_ws(x, header_try = FALSE)
  if (ncol(df2) < 5) {
    stop("comparison_file must have 5 columns or a header with: ",
         paste(req, collapse = ", "))
  }
  names(df2)[1:5] <- req
  df2[, req, drop = FALSE]
}

#' Clean column names
#'
#' Strips CRLF artifacts and trims leading/trailing whitespace from column names.
#'
#' @keywords internal
.clean_colnames <- function(d) {
  nn <- names(d)
  nn <- sub("\r$", "", nn)
  nn <- trimws(nn)
  names(d) <- nn
  d
}

#' Filter dNdS table by NA / max_dnds / optional expression
#'
#' @keywords internal
.filter_dnds <- function(d, dnds_col = "dNdS", filter_expr = NULL, max_dnds = 10) {
  if (!dnds_col %in% names(d)) stop("Missing dnds_col: ", dnds_col)
  d[[dnds_col]] <- suppressWarnings(as.numeric(d[[dnds_col]]))
  keep <- !is.na(d[[dnds_col]]) & is.finite(d[[dnds_col]]) & d[[dnds_col]] < max_dnds

  if (!is.null(filter_expr) && nzchar(filter_expr)) {
    env <- list2env(d, parent = parent.frame())

    expr <- tryCatch(parse(text = filter_expr), error = function(e) e)
    if (inherits(expr, "error")) {
      stop("filter_expr failed to parse: ", conditionMessage(expr))
    }

    ok <- tryCatch(eval(expr, envir = env), error = function(e) e)
    if (inherits(ok, "error")) {
      stop("filter_expr failed to evaluate: ", conditionMessage(ok))
    }

    ok <- as.logical(ok)
    if (length(ok) == 1L) ok <- rep(ok, nrow(d))
    if (length(ok) != nrow(d)) {
      stop("filter_expr must return a logical vector of length 1 or nrow(d).")
    }
    keep <- keep & ok
  }

  d[keep, , drop = FALSE]
}

#' Read and normalize a regions BED-like file (EXPECTED HEADERLESS)
#'
#' Expected (headerless):
#'   col1 = seqname, col2 = start, col3 = end, optional col4 = region_name/label
#'
#' Column selectors may be:
#'   - integer indices (1-based), OR
#'   - column names like "V1", "V2" (if you read elsewhere and name them)
#'
#' Returns data.frame: seqname, start, end, region_name
#'
#' Notes:
#' - Coordinates are treated as inclusive numeric intervals for overlap:
#'     gene_start <= region_end AND gene_end >= region_start
#' - If you provide canonical BED (often 0-based, half-open), convert upstream to match
#'   the coordinate convention in your dN/dS tables.
#'
#' @keywords internal
.read_regions_file <- function(regions_file,
                               region_seq_col   = NULL,
                               region_start_col = NULL,
                               region_end_col   = NULL,
                               region_name_col  = NULL) {
  if (is.null(regions_file)) return(NULL)
  if (!file.exists(regions_file)) stop("regions_file not found: ", regions_file)

  reg_raw <- .read_ws(regions_file, header_try = FALSE)
  if (ncol(reg_raw) < 3) stop("regions_file must have at least 3 columns (headerless): seqname, start, end.")
  names(reg_raw) <- paste0("V", seq_len(ncol(reg_raw)))

  .resolve_idx <- function(sel, default_idx) {
    if (is.null(sel)) return(as.integer(default_idx))
    if (is.numeric(sel) && length(sel) == 1L) {
      idx <- as.integer(sel)
      if (idx < 1L || idx > ncol(reg_raw)) stop("regions_file column index out of range: ", idx)
      return(idx)
    }
    if (is.character(sel) && length(sel) == 1L) {
      nm <- as.character(sel)
      if (!nm %in% names(reg_raw)) stop("regions_file column not found: ", nm, " (expected e.g. 'V1').")
      return(match(nm, names(reg_raw)))
    }
    stop("Invalid regions_file column selector: ", sel)
  }

  seq_idx   <- .resolve_idx(region_seq_col,   1L)
  start_idx <- .resolve_idx(region_start_col, 2L)
  end_idx   <- .resolve_idx(region_end_col,   3L)

  name_idx <- NULL
  if (!is.null(region_name_col)) {
    name_idx <- .resolve_idx(region_name_col, 4L)
  } else if (ncol(reg_raw) >= 4) {
    name_idx <- 4L
  }

  seqv <- as.character(reg_raw[[seq_idx]])
  st0  <- suppressWarnings(as.numeric(reg_raw[[start_idx]]))
  en0  <- suppressWarnings(as.numeric(reg_raw[[end_idx]]))

  if (anyNA(seqv) || any(!nzchar(seqv))) stop("regions_file contains missing/blank seqname values.")
  if (anyNA(st0) || anyNA(en0)) stop("regions_file contains non-numeric start/end values (NA after coercion).")

  lo <- pmin(st0, en0)
  hi <- pmax(st0, en0)
  if (any(lo > hi)) stop("regions_file contains rows with start > end (after normalization).")

  namev <- if (!is.null(name_idx)) as.character(reg_raw[[name_idx]]) else rep("region", nrow(reg_raw))
  namev[is.na(namev) | !nzchar(namev)] <- "region"

  out <- data.frame(
    seqname     = seqv,
    start       = lo,
    end         = hi,
    region_name = namev,
    stringsAsFactors = FALSE
  )
  out
}

#' Label rows as inside/outside regions given explicit coord columns
#'
#' Adds region_status ("region"/"background") based on overlap.
#' Uses data.table::foverlaps if available; otherwise falls back to a loop.
#'
#' Behavior:
#' - Rows with missing/non-finite coords are labeled "background".
#' - start/end are normalized so start <= end.
#'
#' @keywords internal
.label_region_status_coords <- function(d, regions, seq_col, start_col, end_col) {
  if (is.null(regions)) {
    d$region_status <- "global"
    return(d)
  }
  if (!all(c(seq_col, start_col, end_col) %in% names(d))) {
    stop("Missing required columns for region labeling: ",
         paste(c(seq_col, start_col, end_col), collapse = ", "))
  }

  seqv <- as.character(d[[seq_col]])
  st0  <- suppressWarnings(as.numeric(d[[start_col]]))
  en0  <- suppressWarnings(as.numeric(d[[end_col]]))

  st <- pmin(st0, en0)
  en <- pmax(st0, en0)
  ok <- !is.na(seqv) & nzchar(seqv) & is.finite(st) & is.finite(en)

  lab <- rep("background", nrow(d))
  if (!any(ok) || !nrow(regions)) {
    d$region_status <- lab
    return(d)
  }

  has_dt <- requireNamespace("data.table", quietly = TRUE)

  if (has_dt) {
    dt_g <- data.table::data.table(
      idx     = which(ok),
      seqname = seqv[ok],
      start   = st[ok],
      end     = en[ok]
    )
    dt_r <- data.table::data.table(
      seqname = as.character(regions$seqname),
      start   = suppressWarnings(as.numeric(regions$start)),
      end     = suppressWarnings(as.numeric(regions$end))
    )
    dt_r <- dt_r[is.finite(start) & is.finite(end) & !is.na(seqname)]
    if (nrow(dt_r)) {
      data.table::setkey(dt_g, seqname, start, end)
      data.table::setkey(dt_r, seqname, start, end)
      ov <- data.table::foverlaps(dt_g, dt_r, nomatch = 0L)
      in_idx <- unique(ov$idx)
      lab[in_idx] <- "region"
    }
  } else {
    ok_idx <- which(ok)
    for (i in seq_len(nrow(regions))) {
      r <- regions[i, ]
      if (!is.finite(r$start) || !is.finite(r$end) || is.na(r$seqname)) next
      hits <- ok_idx[
        seqv[ok_idx] == r$seqname &
          st[ok_idx] <= r$end &
          en[ok_idx] >= r$start
      ]
      if (length(hits)) lab[hits] <- "region"
    }
  }

  d$region_status <- lab
  d
}

#' Label region_status for standard dndsR pairwise tables using side (query/subject)
#'
#' Requires columns like q_gff_seqname/q_gff_start/q_gff_end or s_* depending on side.
#'
#' @keywords internal
.label_region_status_pairwise <- function(d, regions, side = c("query","subject")) {
  side <- match.arg(side)
  if (is.null(regions)) {
    d$region_status <- "global"
    return(d)
  }
  if (side == "query") {
    seq_col   <- "q_gff_seqname"
    start_col <- "q_gff_start"
    end_col   <- "q_gff_end"
  } else {
    seq_col   <- "s_gff_seqname"
    start_col <- "s_gff_start"
    end_col   <- "s_gff_end"
  }
  .label_region_status_coords(d, regions, seq_col, start_col, end_col)
}

#' Make state labels from numeric dN/dS
#'
#' Uses both neutral bounds and an (optional) explicit positive threshold:
#' - Positive:  x > pos_threshold (if finite), otherwise x > neutral_upper
#' - Purifying: x < neutral_lower
#' - Neutral:   otherwise (neutral band in-between)
#'
#' @keywords internal
.make_states <- function(d, dnds_col, pos_threshold, neutral_lower, neutral_upper) {
  x <- suppressWarnings(as.numeric(d[[dnds_col]]))
  pos_cut <- if (is.finite(pos_threshold)) pos_threshold else neutral_upper

  d$state <- ifelse(x > pos_cut, "Positive",
                    ifelse(x < neutral_lower, "Purifying", "Neutral"))
  d$state <- factor(d$state, levels = c("Positive","Neutral","Purifying"))
  d
}

#' Adjust p-values
#'
#' @keywords internal
.adjust_p <- function(p, method = c("BH","BY","none")) {
  method <- match.arg(method)
  if (method == "none") return(p)
  stats::p.adjust(p, method = method)
}

#' Fit one logistic model and extract group effect for (A vs B)
#'
#' Uses reference=B, coefficient for A is log(odds_A) - log(odds_B)
#'
#' Expects d$y already present (0/1)
#'
#' Adds a guard against (quasi-)complete separation:
#' - returns NULL if y has no variation
#' - returns NULL if either group has only 0s or only 1s for y
#'
#' @keywords internal
.fit_pair <- function(d, group_col, groupA, groupB, random_effect_col = NULL) {
  d2 <- d[d[[group_col]] %in% c(groupA, groupB), , drop = FALSE]
  if (!nrow(d2)) return(NULL)

  if (!"y" %in% names(d2)) stop("Internal error: d$y is required for model fitting.")
  d2$y <- as.integer(d2$y)

  # Must have both outcomes present
  if (length(unique(d2$y[is.finite(d2$y)])) < 2L) return(NULL)

  # Separation guard: each group must have both outcomes (force levels 0/1)
  tt <- table(d2[[group_col]], factor(d2$y, levels = c(0, 1)))
  if (any(tt[, "0"] == 0L | tt[, "1"] == 0L)) return(NULL)

  d2[[group_col]] <- stats::relevel(factor(d2[[group_col]]), ref = groupB)

  use_re <- !is.null(random_effect_col) &&
    nzchar(random_effect_col) &&
    random_effect_col %in% names(d2) &&
    any(!is.na(d2[[random_effect_col]]))

  if (use_re) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("random_effect_col set but package 'lme4' is not installed.")
    }
    d3 <- d2[!is.na(d2[[random_effect_col]]), , drop = FALSE]
    if (!nrow(d3)) return(NULL)

    # need >=2 random-effect levels to be meaningful/stable
    if (length(unique(as.character(d3[[random_effect_col]]))) < 2L) {
      use_re <- FALSE
    } else {
      form <- stats::as.formula(sprintf("y ~ %s + (1|%s)", group_col, random_effect_col))
      fit  <- suppressWarnings(lme4::glmer(form, data = d3, family = stats::binomial(link = "logit")))
      coefs <- try(suppressWarnings(lme4::fixef(fit)), silent = TRUE)
      vcovm <- try(stats::vcov(fit), silent = TRUE)
      if (inherits(coefs, "try-error") || inherits(vcovm, "try-error")) return(NULL)

      tt2 <- setdiff(names(coefs), "(Intercept)")
      if (length(tt2) != 1) return(NULL)
      term <- tt2[1]

      est <- unname(coefs[term])

      # variance index first, then sqrt (and guard)
      se  <- sqrt(diag(vcovm)[term])
      if (!is.finite(se) || se <= 0) return(NULL)

      z   <- stats::qnorm(0.975)
      lo  <- est - z * se
      hi  <- est + z * se
      zval <- est / se
      pval <- 2 * stats::pnorm(-abs(zval))

      return(list(model = "glmer",
                  log_or = est,
                  se = unname(se),
                  p_value = unname(pval),
                  log_or_ci_lower = unname(lo),
                  log_or_ci_upper = unname(hi)))
    }
  }

  # GLM fallback
  form <- stats::as.formula(sprintf("y ~ %s", group_col))
  fit  <- suppressWarnings(stats::glm(form, data = d2, family = stats::binomial(link = "logit")))
  sm   <- summary(fit)$coefficients

  rn <- rownames(sm)
  rn <- rn[rn != "(Intercept)"]
  if (!length(rn)) return(NULL)

  if (length(rn) != 1) return(NULL)
  rr <- rn[1]

  est <- sm[rr, "Estimate"]
  se  <- sm[rr, "Std. Error"]
  pval <- sm[rr, "Pr(>|z|)"]

  if (!is.finite(se) || se <= 0) return(NULL)

  z <- stats::qnorm(0.975)
  lo <- est - z * se
  hi <- est + z * se

  list(model = "glm",
       log_or = unname(est),
       se = unname(se),
       p_value = unname(pval),
       log_or_ci_lower = unname(lo),
       log_or_ci_upper = unname(hi))
}

#' Forest plot for state contrasts
#'
#' Orders contrasts *within each facet* for readability (avoids global factor ordering).
#'
#' @keywords internal
.plot_state_forest <- function(res,
                               out_path,
                               base_family = "Liberation Sans",
                               facet_by = c("scope","state"),
                               point_size = 2.8,
                               line_width = 0.9,
                               dpi = 600) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  facet_by <- match.arg(facet_by)

  # OR scale
  res$odds_ratio  <- exp(res$log_or)
  res$or_ci_lower <- exp(res$log_or_ci_lower)
  res$or_ci_upper <- exp(res$log_or_ci_upper)

  # Build a facet-specific y variable so ordering is per-facet, not global.
  # We keep displayed labels as the original contrast_label.
  facet_id <- if (facet_by == "scope") paste(res$state, res$scope, sep = "||") else as.character(res$state)
  y_key <- paste(res$contrast_label, facet_id, sep = " @ ")

  # Order within each facet by odds_ratio (then name as tiebreak)
  ord_df <- data.frame(y_key = y_key,
                       facet_id = facet_id,
                       odds_ratio = res$odds_ratio,
                       contrast_label = as.character(res$contrast_label),
                       stringsAsFactors = FALSE)
  ord_df <- ord_df[order(ord_df$facet_id, ord_df$odds_ratio, ord_df$contrast_label), , drop = FALSE]

  levs <- character(0)
  for (fid in unique(ord_df$facet_id)) {
    kk <- ord_df[ord_df$facet_id == fid, , drop = FALSE]
    levs <- c(levs, unique(kk$y_key))
  }

  res$y_key <- factor(y_key, levels = rev(levs)) # rev to place largest OR near top visually

  gg <- ggplot2::ggplot(res, ggplot2::aes(x = odds_ratio, y = y_key)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = or_ci_lower, xmax = or_ci_upper),
      height = 0.2,
      linewidth = line_width
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_discrete(labels = function(x) sub(" @ .*$", "", x)) +
    ggplot2::labs(x = "Odds ratio", y = NULL) +
    ggplot2::theme_classic(base_size = 14, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 11)
    )

  if (facet_by == "scope") {
    gg <- gg + ggplot2::facet_grid(state ~ scope, scales = "free_y")
  } else {
    gg <- gg + ggplot2::facet_wrap(~ state, scales = "free_y")
  }

  ext <- tolower(tools::file_ext(out_path))
  if (ext %in% c("png","tiff","jpeg","jpg")) {
    ggplot2::ggsave(out_path, gg, width = 11, height = 8, dpi = dpi)
  } else {
    ggplot2::ggsave(out_path, gg, width = 11, height = 8)
  }
  invisible(NULL)
}

# ---- gene-mode helpers ----

#' Infer group from seqname (current default: capture terminal B/C/D)
#'
#' NOTE: This helper is intentionally conservative. In gene mode, inference failures
#' will STOP with examples rather than silently dropping genes. If you want general
#' haplotype/subgenome handling, prefer group_mode="custom" + gene_group_col.
#'
#' @keywords internal
.infer_group_from_seqname <- function(x) {
  x <- as.character(x)
  sg <- sub("^.*([BCD])$", "\\1", x)
  sg[!grepl("^[BCD]$", sg)] <- NA_character_
  sg
}

#' Build gene rows from pairwise table for query/subject sides
#'
#' Adds:
#'   gene_id, side, seqname, start, end, dNdS_pair, comparison, gene_key
#'
#' @keywords internal
.pairwise_to_gene_rows <- function(d, comp_name, dnds_col, focal_sides = c("query","subject"), random_effect_col = NULL) {
  focal_sides <- match.arg(focal_sides, choices = c("query","subject"), several.ok = TRUE)
  out <- list()

  re_ok <- !is.null(random_effect_col) && nzchar(random_effect_col) && random_effect_col %in% names(d)

  if ("query" %in% focal_sides) {
    need <- c("query_id","q_gff_seqname","q_gff_start","q_gff_end", dnds_col)
    miss <- setdiff(need, names(d))
    if (length(miss)) stop("Missing in pairwise table (query side): ", paste(miss, collapse = ", "))

    tmp <- data.frame(
      gene_id    = as.character(d$query_id),
      side       = "query",
      seqname    = as.character(d$q_gff_seqname),
      start      = suppressWarnings(as.numeric(d$q_gff_start)),
      end        = suppressWarnings(as.numeric(d$q_gff_end)),
      dNdS_pair  = suppressWarnings(as.numeric(d[[dnds_col]])),
      comparison = as.character(comp_name),
      stringsAsFactors = FALSE
    )
    if (re_ok) tmp[[random_effect_col]] <- d[[random_effect_col]]

    tmp$gene_key <- paste(tmp$comparison, tmp$side, tmp$gene_id, sep = "|")
    out[[length(out) + 1]] <- tmp
  }

  if ("subject" %in% focal_sides) {
    need <- c("subject_id","s_gff_seqname","s_gff_start","s_gff_end", dnds_col)
    miss <- setdiff(need, names(d))
    if (length(miss)) stop("Missing in pairwise table (subject side): ", paste(miss, collapse = ", "))

    tmp <- data.frame(
      gene_id    = as.character(d$subject_id),
      side       = "subject",
      seqname    = as.character(d$s_gff_seqname),
      start      = suppressWarnings(as.numeric(d$s_gff_start)),
      end        = suppressWarnings(as.numeric(d$s_gff_end)),
      dNdS_pair  = suppressWarnings(as.numeric(d[[dnds_col]])),
      comparison = as.character(comp_name),
      stringsAsFactors = FALSE
    )
    if (re_ok) tmp[[random_effect_col]] <- d[[random_effect_col]]

    tmp$gene_key <- paste(tmp$comparison, tmp$side, tmp$gene_id, sep = "|")
    out[[length(out) + 1]] <- tmp
  }

  do.call(rbind, out)
}

#' Aggregate gene rows to gene-level table (per gene_key + group)
#'
#' Random effects:
#' - Only retained if random_effect_col is constant within each gene_key.
#' - If it varies within gene_key, it is disabled with a warning.
#'
#' Safety:
#' - In group_mode="subgenome", if any groups cannot be inferred, this function STOPs
#'   with example seqnames to prevent silent gene loss.
#'
#' Adds back `comparison` and `side` parsed from gene_key for debugging/sanity checks.
#'
#' @keywords internal
.aggregate_gene_table <- function(g,
                                  group_mode = c("subgenome","custom"),
                                  group_col = NULL,
                                  agg_fun = c("median","mean"),
                                  min_pairs = 2,
                                  random_effect_col = NULL) {
  group_mode <- match.arg(group_mode)
  agg_fun    <- match.arg(agg_fun)

  g <- g[!is.na(g$gene_key) & nzchar(g$gene_key) &
           !is.na(g$dNdS_pair) & is.finite(g$dNdS_pair), , drop = FALSE]
  if (!nrow(g)) return(g)

  # normalize gene coords
  g$start <- suppressWarnings(as.numeric(g$start))
  g$end   <- suppressWarnings(as.numeric(g$end))
  st0 <- g$start
  en0 <- g$end
  ok_xy <- is.finite(st0) & is.finite(en0)
  g$start[ok_xy] <- pmin(st0[ok_xy], en0[ok_xy])
  g$end[ok_xy]   <- pmax(st0[ok_xy], en0[ok_xy])

  if (group_mode == "subgenome") {
    g$group <- .infer_group_from_seqname(g$seqname)

    if (any(is.na(g$group) | !nzchar(g$group))) {
      bad_seq <- unique(as.character(g$seqname[is.na(g$group) | !nzchar(g$group)]))
      bad_seq <- bad_seq[!is.na(bad_seq) & nzchar(bad_seq)]
      msg <- paste0(
        "Failed to infer 'group' from seqname for ", length(bad_seq), " seqnames in gene mode.\n",
        "This would silently drop genes; stopping for safety.\n\n",
        "Provide group_mode='custom' and gene_group_col (recommended for general haplotypes/subgenomes),\n",
        "or adjust your seqnames to match the inference rule.\n\n",
        "Examples:\n  ", paste(utils::head(bad_seq, 5), collapse = "\n  ")
      )
      stop(msg)
    }
  } else {
    if (is.null(group_col) || !nzchar(group_col)) stop("group_col is required when group_mode='custom'.")
    if (!group_col %in% names(g)) stop("group_col '", group_col, "' not found in gene rows.")
    g$group <- as.character(g[[group_col]])
  }

  g <- g[!is.na(g$group) & nzchar(g$group), , drop = FALSE]
  if (!nrow(g)) return(g)

  fun <- if (agg_fun == "median") stats::median else mean

  # aggregate dNdS
  agg <- stats::aggregate(
    dNdS_pair ~ gene_key + gene_id + group,
    data = g,
    FUN  = function(x) fun(x[is.finite(x)], na.rm = TRUE)
  )
  names(agg)[names(agg) == "dNdS_pair"] <- "dNdS_agg"

  # n_pairs
  n_pairs <- stats::aggregate(
    dNdS_pair ~ gene_key + gene_id + group,
    data = g,
    FUN = function(x) sum(is.finite(x))
  )
  names(n_pairs)[names(n_pairs) == "dNdS_pair"] <- "n_pairs"

  agg <- merge(agg, n_pairs, by = c("gene_key","gene_id","group"), all.x = TRUE)

  # add comparison + side parsed from gene_key: "<comparison>|<side>|<gene_id>"
  parts <- strsplit(as.character(agg$gene_key), "\\|", fixed = FALSE)
  agg$comparison <- vapply(parts, function(p) if (length(p) >= 1) p[1] else NA_character_, character(1))
  agg$side       <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, character(1))

  # attach a random-effect column if requested and valid for gene-level
  if (!is.null(random_effect_col) && nzchar(random_effect_col) && random_effect_col %in% names(g)) {
    re_counts <- tapply(g[[random_effect_col]], g$gene_key, function(x) length(unique(na.omit(x))))
    re_counts[is.na(re_counts)] <- 0L
    bad <- names(re_counts)[re_counts > 1L]

    if (length(bad)) {
      warning("random_effect_col varies within gene_key in gene mode; disabling random effects for gene models.")
    } else {
      re_map <- g[, c("gene_key", random_effect_col), drop = FALSE]
      re_map <- re_map[!is.na(re_map[[random_effect_col]]), , drop = FALSE]
      re_map <- re_map[!duplicated(re_map$gene_key), , drop = FALSE]
      agg <- merge(agg, re_map, by = "gene_key", all.x = TRUE)
    }
  }

  agg <- agg[is.finite(agg$dNdS_agg) & agg$n_pairs >= min_pairs, , drop = FALSE]
  agg
}

# -----------------------------
# Exported function
# -----------------------------

#' dN/dS state contrasts via logistic regression (pairwise or gene-level) + optional forest plot
#'
#' Fits logistic regression models comparing the odds of each dN/dS "state"
#' (Positive / Neutral / Purifying) between groups (e.g., subgenomes B/C/D).
#'
#' Two analysis levels are supported:
#' \itemize{
#'   \item \code{level = "pairwise"}: each ortholog-pair row is an observation.
#'   \item \code{level = "gene"}: converts pairwise rows to gene-side observations
#'     (query/subject), aggregates across partners per gene, then models gene-level states.
#' }
#'
#' Region filtering (optional) treats coordinates as inclusive numeric intervals:
#' a row overlaps a region if \code{start <= region_end AND end >= region_start}.
#' Ensure your \code{regions_file} coordinate system matches your dN/dS tables.
#'
#' @section Core mode:
#' @param level Character. Analysis level: \code{"pairwise"} or \code{"gene"}.
#'
#' @section Input selection (single vs batch):
#' @param dnds_table_file Character. Path to a single input TSV.
#' In pairwise mode, this is a pairwise dN/dS table containing \code{dnds_col} and \code{group_col}.
#' In gene mode, this can be used as a single-file alias for \code{dnds_annot_files}.
#' @param comparison_file Character or data.frame. Batch-mode comparisons table.
#' If a file path, it must contain columns:
#' \code{comparison_name, query_fasta, query_gff, subject_fasta, subject_gff}.
#' If provided, runs once per comparison subdirectory under \code{output_dir}.
#' @param output_dir Character. Base output directory for batch mode (default: \code{getwd()}).
#' @param in_suffix Character or NULL. Batch mode only: relative filename inside each
#' \code{file.path(output_dir, comparison_name)} directory. If NULL, defaults to
#' \code{"<comparison_name>_dnds_annot.tsv"}.
#'
#' @section Core columns:
#' @param dnds_col Character. Column name containing dN/dS values (default \code{"dNdS"}).
#' @param group_col Character. Pairwise mode only: column name defining groups
#' (default \code{"group"}).
#'
#' @section State definitions and filtering:
#' @param pos_threshold Numeric. Positive selection threshold:
#' \code{x > pos_threshold} is Positive. If not finite, Positive is \code{x > neutral_upper}.
#' Default: 1.
#' @param neutral_lower Numeric. Lower bound for Neutral band (default 0.9). Values below are Purifying.
#' @param neutral_upper Numeric. Upper bound for Neutral band (default 1.1). Values above are Positive
#' (unless overridden by \code{pos_threshold}).
#' @param max_dnds Numeric. Maximum allowed dN/dS to retain (rows with \code{dnds_col >= max_dnds}
#' or non-finite values are dropped). Default: 10.
#' @param filter_expr Character or NULL. Optional R expression (as text) evaluated in the table’s
#' column environment; must return a logical vector of length 1 or \code{nrow(d)}.
#'
#' @section Modeling options:
#' @param random_effect_col Character or NULL. Optional column name to use as a random intercept
#' (mixed model via \code{lme4::glmer}). If NULL/unavailable, falls back to \code{glm}.
#' @param min_n_per_group Integer. Minimum number of observations per group to keep that group
#' for modeling within each scope (default: 50).
#' @param fdr_method Character. Multiple-testing correction method applied within each
#' \code{state x scope} block. One of \code{"BH"}, \code{"BY"}, or \code{"none"}.
#'
#' @section Regions and scopes:
#' @param regions_file Character or NULL. Optional BED-like file (expected headerless) specifying regions.
#' Must contain at least 3 columns: seqname, start, end; optional 4th column is region name/label.
#' @param region_seq_col,region_start_col,region_end_col,region_name_col Integer/Character/NULL.
#' Optional selectors for columns in \code{regions_file}. May be 1-based indices (e.g., 1,2,3,4)
#' or names like \code{"V1"}. If NULL, defaults to 1/2/3 and uses column 4 as name if present.
#' @param side Character. Pairwise mode only: which side’s coordinates to use when labeling rows
#' as inside/outside regions. One of \code{"query"} or \code{"subject"}.
#' @param scopes Character vector or NULL. Which scopes to analyze: \code{"global"} and optionally \code{"region"}.
#' If NULL, uses \code{"global"} and adds \code{"region"} when \code{regions_file} is provided.
#'
#' @section Output and plotting:
#' @param make_plots Logical. If TRUE, writes forest plots (PDF + PNG) alongside TSV results.
#' @param base_family Character. Base font family used by ggplot theme (default: \code{"Liberation Sans"}).
#' @param out_prefix Character or NULL. Output filename prefix. If NULL, uses the inferred label
#' (comparison name in batch mode; input basename in single mode).
#'
#' @section Gene-mode only:
#' @param dnds_annot_files Character vector or NULL. Gene mode only: one or more paths to
#' \code{*_dnds_annot.tsv} tables to stack before converting to gene rows.
#' @param focal_sides Character vector. Gene mode only: which sides to include when converting pairwise
#' rows to gene observations. Any of \code{"query"} and/or \code{"subject"}.
#' @param group_mode Character. Gene mode only: how to define groups for gene observations.
#' \code{"subgenome"} infers group from \code{seqname} using a conservative terminal B/C/D rule;
#' \code{"custom"} uses \code{gene_group_col}.
#' @param gene_group_col Character or NULL. Gene mode only: when \code{group_mode="custom"},
#' this is the column in gene rows used as group label.
#' @param agg_fun Character. Gene mode only: how to aggregate per-gene partner dN/dS values.
#' One of \code{"median"} or \code{"mean"}.
#' @param min_pairs Integer. Gene mode only: minimum number of finite partner pairs contributing to a gene
#' before retaining that gene-level observation (default: 2).
#' @param gene_out_dir Character. Batch + gene mode only: directory name under \code{output_dir}
#' where gene-level results are written (default: \code{"gene_level_state_contrast"}).
#'
#' @return Invisibly returns a list. In single mode:
#' \code{list(results = <data.frame>, paths = <character>)}.
#' In batch mode:
#' \code{list(results = <named list>, paths = <character>)}.
#' In gene mode, also includes \code{gene_table} (per-scope gene-level tables) when available.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Pairwise (single)
#' dnds_state_contrast(
#'   level = "pairwise",
#'   dnds_table_file = "BvC_dnds_annot.tsv",
#'   group_col = "subgenome",
#'   dnds_col = "dNdS"
#' )
#'
#' # Gene (single file alias)
#' dnds_state_contrast(
#'   level = "gene",
#'   dnds_table_file = "BvC_dnds_annot.tsv",
#'   group_mode = "subgenome",
#'   focal_sides = c("query","subject")
#' )
#' }
dnds_state_contrast <- function(level            = c("pairwise","gene"),
                                # shared / pairwise inputs
                                dnds_table_file   = NULL,
                                comparison_file   = NULL,
                                output_dir        = getwd(),
                                in_suffix         = NULL,
                                dnds_col          = "dNdS",
                                group_col         = "group",
                                random_effect_col = NULL,
                                pos_threshold     = 1,
                                neutral_lower     = 0.9,
                                neutral_upper     = 1.1,
                                max_dnds          = 10,
                                filter_expr       = NULL,
                                min_n_per_group   = 50,
                                fdr_method        = c("BH","BY","none"),
                                regions_file       = NULL,
                                region_seq_col    = NULL,
                                region_start_col  = NULL,
                                region_end_col    = NULL,
                                region_name_col   = NULL,
                                side              = c("query","subject"),
                                scopes            = NULL,
                                make_plots        = TRUE,
                                base_family       = "Liberation Sans",
                                out_prefix        = NULL,
                                # gene-mode only
                                dnds_annot_files  = NULL,
                                focal_sides       = c("query","subject"),
                                group_mode        = c("subgenome","custom"),
                                gene_group_col    = NULL,
                                agg_fun           = c("median","mean"),
                                min_pairs         = 2,
                                gene_out_dir      = "gene_level_state_contrast") {

  level      <- match.arg(level)
  fdr_method <- match.arg(fdr_method)
  side       <- match.arg(side)
  group_mode <- match.arg(group_mode)
  agg_fun    <- match.arg(agg_fun)
  focal_sides <- match.arg(focal_sides, choices = c("query","subject"), several.ok = TRUE)

  # sanity checks for thresholds
  if (!is.finite(neutral_lower) || !is.finite(neutral_upper) || neutral_lower > neutral_upper) {
    stop("Require neutral_lower <= neutral_upper (both finite).")
  }

  # Was previously a hard stop; now warn (overlap is allowed, Positive overrides Neutral).
  if (is.finite(pos_threshold) && pos_threshold < neutral_upper) {
    warning("pos_threshold (", pos_threshold, ") is < neutral_upper (", neutral_upper, "). ",
            "Positive classification will override Neutral for x > pos_threshold. ",
            "If you intended non-overlapping states, set pos_threshold >= neutral_upper.")
  }

  regions <- .read_regions_file(regions_file,
                                region_seq_col   = region_seq_col,
                                region_start_col = region_start_col,
                                region_end_col   = region_end_col,
                                region_name_col  = region_name_col)

  # determine scopes
  if (is.null(scopes)) {
    scopes <- "global"
    if (!is.null(regions)) scopes <- c(scopes, "region")
  } else {
    scopes <- match.arg(scopes, choices = c("global","region"), several.ok = TRUE)
    if ("region" %in% scopes && is.null(regions)) {
      stop("scopes includes 'region' but regions_file is NULL.")
    }
  }

  .read_tsv <- function(path) {
    d <- utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                           quote = "", comment.char = "", check.names = FALSE)
    .clean_colnames(d)
  }

  # -----------------------------
  # Pairwise runner (table already has group_col)
  # -----------------------------
  .run_one_pairwise_table <- function(d, label, out_dir) {
    if (!group_col %in% names(d)) stop("Missing group_col: ", group_col)
    if (!dnds_col %in% names(d)) stop("Missing dnds_col: ", dnds_col)

    d <- .filter_dnds(d, dnds_col = dnds_col, filter_expr = filter_expr, max_dnds = max_dnds)
    if (!nrow(d)) {
      warning("[dnds_state_contrast/pairwise] No rows after filtering for: ", label)
      return(list(results = NULL, paths = character(0)))
    }

    scopes_pw <- scopes
    if (!is.null(regions) && "region" %in% scopes_pw) {
      req_cols <- if (side == "query") {
        c("q_gff_seqname", "q_gff_start", "q_gff_end")
      } else {
        c("s_gff_seqname", "s_gff_start", "s_gff_end")
      }
      if (!all(req_cols %in% names(d))) {
        warning("[dnds_state_contrast/pairwise] regions_file provided but missing required columns for side='",
                side, "': ", paste(req_cols, collapse = ", "),
                ". Disabling region scope for: ", label)
        d$region_status <- "global"
        scopes_pw <- setdiff(scopes_pw, "region")
      } else {
        d <- .label_region_status_pairwise(d, regions = regions, side = side)
      }
    } else {
      d$region_status <- "global"
    }

    d <- .make_states(d, dnds_col = dnds_col,
                      pos_threshold = pos_threshold,
                      neutral_lower = neutral_lower,
                      neutral_upper = neutral_upper)

    all_res <- list()

    for (scp in scopes_pw) {
      d_sc <- if (scp == "global") d else d[d$region_status == "region", , drop = FALSE]

      tab <- table(d_sc[[group_col]])
      keep_groups <- names(tab)[tab >= min_n_per_group]
      d_sc <- d_sc[d_sc[[group_col]] %in% keep_groups, , drop = FALSE]
      groups <- unique(as.character(d_sc[[group_col]]))

      if (length(groups) < 2) {
        warning("[dnds_state_contrast/pairwise] Not enough groups for scope='", scp, "' in: ", label)
        next
      }

      pairs <- utils::combn(sort(groups), 2, simplify = FALSE)
      states <- levels(d_sc$state)

      for (st in states) {
        d_sc$y <- as.integer(d_sc$state == st)
        for (p in pairs) {
          A <- p[1]; B <- p[2]
          fit <- .fit_pair(d_sc, group_col = group_col, groupA = A, groupB = B,
                           random_effect_col = random_effect_col)
          if (is.null(fit)) next

          all_res[[length(all_res) + 1]] <- data.frame(
            label        = label,
            level        = "pairwise",
            scope        = scp,
            side         = if (scp == "region") side else NA_character_,
            state        = st,
            groupA       = A,
            groupB       = B,
            contrast     = paste0(A, " vs ", B),
            contrast_label = paste0(A, " vs ", B),
            model        = fit$model,
            log_or       = fit$log_or,
            se           = fit$se,
            log_or_ci_lower = fit$log_or_ci_lower,
            log_or_ci_upper = fit$log_or_ci_upper,
            p_value      = fit$p_value,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (!length(all_res)) {
      warning("[dnds_state_contrast/pairwise] No models fit for: ", label)
      return(list(results = NULL, paths = character(0)))
    }

    res <- do.call(rbind, all_res)

    res$p_adj <- NA_real_
    key <- paste(res$state, res$scope, sep = "||")
    for (k in unique(key)) {
      idx <- which(key == k)
      res$p_adj[idx] <- .adjust_p(res$p_value[idx], method = fdr_method)
    }

    res$odds_ratio  <- exp(res$log_or)
    res$or_ci_lower <- exp(res$log_or_ci_lower)
    res$or_ci_upper <- exp(res$log_or_ci_upper)

    res$pos_threshold <- pos_threshold
    res$neutral_lower <- neutral_lower
    res$neutral_upper <- neutral_upper
    res$max_dnds      <- max_dnds
    res$min_n_per_group <- min_n_per_group
    res$fdr_method    <- fdr_method

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    prefix <- if (!is.null(out_prefix) && nzchar(out_prefix)) out_prefix else label

    out_tsv <- file.path(out_dir, sprintf("%s_dnds_state_contrasts.tsv", prefix))
    utils::write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    out_paths <- c(out_tsv)

    if (isTRUE(make_plots)) {
      out_pdf <- file.path(out_dir, sprintf("%s_dnds_state_contrasts_forest.pdf", prefix))
      out_png <- file.path(out_dir, sprintf("%s_dnds_state_contrasts_forest.png", prefix))
      .plot_state_forest(res, out_pdf, base_family = base_family, facet_by = "scope")
      .plot_state_forest(res, out_png, base_family = base_family, facet_by = "scope", dpi = 600)
      out_paths <- c(out_paths, out_pdf, out_png)
    }

    list(results = res, paths = out_paths)
  }

  # -----------------------------
  # Gene-mode runner
  # -----------------------------
  .run_one_gene_from_files <- function(files, label, out_dir) {
    gene_rows <- list()
    used_names <- character(0)

    for (i in seq_along(files)) {
      f <- files[i]
      if (!file.exists(f)) {
        warning("[dnds_state_contrast/gene] Missing file (skipping): ", f)
        next
      }

      comp_name <- sub("_dnds_annot\\.tsv$", "", basename(f))
      comp_name <- sub("\\.tsv$", "", comp_name)

      if (comp_name %in% used_names) {
        comp_name <- make.unique(c(used_names, comp_name))[length(used_names) + 1L]
        warning("[dnds_state_contrast/gene] Duplicate basename detected; using unique comparison id: ", comp_name)
      }
      used_names <- c(used_names, comp_name)

      d <- .read_tsv(f)
      if (!nrow(d)) next
      if (!dnds_col %in% names(d)) stop("Missing dnds_col '", dnds_col, "' in: ", f)

      d <- .filter_dnds(d, dnds_col = dnds_col, filter_expr = filter_expr, max_dnds = max_dnds)
      if (!nrow(d)) next

      g <- .pairwise_to_gene_rows(
        d,
        comp_name = comp_name,
        dnds_col = dnds_col,
        focal_sides = focal_sides,
        random_effect_col = random_effect_col
      )
      if (!nrow(g)) next

      g <- .label_region_status_coords(g, regions = regions, seq_col = "seqname", start_col = "start", end_col = "end")
      gene_rows[[length(gene_rows) + 1]] <- g
    }

    if (!length(gene_rows)) {
      warning("[dnds_state_contrast/gene] No usable rows after filtering for: ", label)
      return(list(results = NULL, paths = character(0), gene_table = NULL))
    }

    g_all <- do.call(rbind, gene_rows)

    all_res <- list()
    gene_tbl_by_scope <- list()

    for (scp in scopes) {
      g_sc <- if (scp == "global") g_all else g_all[g_all$region_status == "region", , drop = FALSE]
      if (!nrow(g_sc)) next

      gene_tbl <- .aggregate_gene_table(
        g_sc,
        group_mode = group_mode,
        group_col  = gene_group_col,
        agg_fun    = agg_fun,
        min_pairs  = min_pairs,
        random_effect_col = random_effect_col
      )
      if (!nrow(gene_tbl)) next

      tab <- table(gene_tbl$group)
      keep_groups <- names(tab)[tab >= min_n_per_group]
      gene_tbl <- gene_tbl[gene_tbl$group %in% keep_groups, , drop = FALSE]
      groups <- unique(as.character(gene_tbl$group))
      if (length(groups) < 2) {
        warning("[dnds_state_contrast/gene] Not enough groups after min_n_per_group in scope='", scp, "' for: ", label)
        next
      }

      gene_tbl <- .make_states(gene_tbl, dnds_col = "dNdS_agg",
                               pos_threshold = pos_threshold,
                               neutral_lower = neutral_lower,
                               neutral_upper = neutral_upper)

      gene_tbl$scope <- scp
      gene_tbl_by_scope[[scp]] <- gene_tbl

      pairs <- utils::combn(sort(groups), 2, simplify = FALSE)
      states <- levels(gene_tbl$state)

      re_for_gene <- NULL
      if (!is.null(random_effect_col) && nzchar(random_effect_col) && random_effect_col %in% names(gene_tbl)) {
        re_for_gene <- random_effect_col
      }

      for (st in states) {
        gene_tbl$y <- as.integer(gene_tbl$state == st)
        for (p in pairs) {
          A <- p[1]; B <- p[2]
          fit <- .fit_pair(gene_tbl, group_col = "group", groupA = A, groupB = B,
                           random_effect_col = re_for_gene)
          if (is.null(fit)) next

          all_res[[length(all_res) + 1]] <- data.frame(
            label        = label,
            level        = "gene",
            scope        = scp,
            side         = NA_character_,
            state        = st,
            groupA       = A,
            groupB       = B,
            contrast     = paste0(A, " vs ", B),
            contrast_label = paste0(A, " vs ", B),
            model        = fit$model,
            log_or       = fit$log_or,
            se           = fit$se,
            log_or_ci_lower = fit$log_or_ci_lower,
            log_or_ci_upper = fit$log_or_ci_upper,
            p_value      = fit$p_value,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (!length(all_res)) {
      warning("[dnds_state_contrast/gene] No models fit for: ", label)
      return(list(results = NULL, paths = character(0), gene_table = gene_tbl_by_scope))
    }

    res <- do.call(rbind, all_res)

    res$p_adj <- NA_real_
    key <- paste(res$state, res$scope, sep = "||")
    for (k in unique(key)) {
      idx <- which(key == k)
      res$p_adj[idx] <- .adjust_p(res$p_value[idx], method = fdr_method)
    }

    res$odds_ratio  <- exp(res$log_or)
    res$or_ci_lower <- exp(res$log_or_ci_lower)
    res$or_ci_upper <- exp(res$log_or_ci_upper)

    res$pos_threshold <- pos_threshold
    res$neutral_lower <- neutral_lower
    res$neutral_upper <- neutral_upper
    res$max_dnds      <- max_dnds
    res$min_n_per_group <- min_n_per_group
    res$fdr_method    <- fdr_method
    res$agg_fun       <- agg_fun
    res$min_pairs     <- min_pairs
    res$group_mode    <- group_mode
    res$focal_sides   <- paste(focal_sides, collapse = ",")

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    prefix <- if (!is.null(out_prefix) && nzchar(out_prefix)) out_prefix else label

    out_tsv <- file.path(out_dir, sprintf("%s_dnds_state_contrasts.tsv", prefix))
    utils::write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    out_paths <- c(out_tsv)

    for (scp in names(gene_tbl_by_scope)) {
      gt <- gene_tbl_by_scope[[scp]]
      out_gt <- file.path(out_dir, sprintf("%s_gene_table_%s.tsv", prefix, scp))
      utils::write.table(gt, out_gt, sep = "\t", quote = FALSE, row.names = FALSE)
      out_paths <- c(out_paths, out_gt)
    }

    if (isTRUE(make_plots)) {
      out_pdf <- file.path(out_dir, sprintf("%s_dnds_state_contrasts_forest.pdf", prefix))
      out_png <- file.path(out_dir, sprintf("%s_dnds_state_contrasts_forest.png", prefix))
      .plot_state_forest(res, out_pdf, base_family = base_family, facet_by = "scope")
      .plot_state_forest(res, out_png, base_family = base_family, facet_by = "scope", dpi = 600)
      out_paths <- c(out_paths, out_pdf, out_png)
    }

    list(results = res, paths = out_paths, gene_table = gene_tbl_by_scope)
  }

  # --------------------
  # Mode dispatch
  # --------------------

  # ---- batch mode ----
  if (!is.null(comparison_file)) {
    cf <- .read_comparisons(comparison_file)
    outs_all  <- list()
    paths_all <- character(0)
    gene_all  <- list()

    for (i in seq_len(nrow(cf))) {
      comp <- as.character(cf$comparison_name[i])
      comp_dir <- file.path(output_dir, comp)

      in_file <- if (!is.null(in_suffix) && nzchar(in_suffix)) {
        file.path(comp_dir, in_suffix)
      } else {
        file.path(comp_dir, sprintf("%s_dnds_annot.tsv", comp))
      }

      if (!file.exists(in_file)) {
        warning("[dnds_state_contrast] Missing input table: ", in_file)
        next
      }

      if (level == "pairwise") {
        d <- .read_tsv(in_file)
        run <- .run_one_pairwise_table(d, label = comp, out_dir = comp_dir)
        outs_all[[comp]] <- run$results
        paths_all <- c(paths_all, run$paths)
      } else {
        out_dir_gene <- file.path(output_dir, gene_out_dir, comp)
        run <- .run_one_gene_from_files(files = c(in_file), label = comp, out_dir = out_dir_gene)
        outs_all[[comp]] <- run$results
        gene_all[[comp]] <- run$gene_table
        paths_all <- c(paths_all, run$paths)
      }
    }

    if (level == "gene") {
      return(invisible(list(results = outs_all, paths = paths_all, gene_table = gene_all)))
    }
    return(invisible(list(results = outs_all, paths = paths_all)))
  }

  # ---- single mode ----
  if (level == "pairwise") {
    if (is.null(dnds_table_file)) stop("Provide dnds_table_file (single) or comparison_file (batch).")
    if (!file.exists(dnds_table_file)) stop("dnds_table_file not found: ", dnds_table_file)

    d <- .read_tsv(dnds_table_file)
    label <- if (!is.null(out_prefix) && nzchar(out_prefix)) out_prefix else sub("\\.tsv$", "", basename(dnds_table_file))
    out_dir <- dirname(dnds_table_file)

    return(invisible(.run_one_pairwise_table(d, label = label, out_dir = out_dir)))
  }

  # gene single mode
  files <- character(0)
  if (!is.null(dnds_annot_files)) {
    files <- as.character(dnds_annot_files)
    if (any(!file.exists(files))) stop("Some dnds_annot_files do not exist.")
  } else if (!is.null(dnds_table_file)) {
    files <- c(dnds_table_file)
    if (!file.exists(files[1])) stop("dnds_table_file not found: ", files[1])
  } else {
    stop("Gene mode requires dnds_annot_files (or dnds_table_file as a single-file alias), or comparison_file for batch.")
  }

  label <- if (!is.null(out_prefix) && nzchar(out_prefix)) out_prefix else {
    sub("_dnds_annot\\.tsv$", "", sub("\\.tsv$", "", basename(files[1])))
  }
  out_dir <- dirname(files[1])

  invisible(.run_one_gene_from_files(files = files, label = label, out_dir = out_dir))
}
