# dnds_state_contrast.R
# dN/dS state contrasts via logistic regression (optional mixed model) + forest plot
#
# Merged "pairwise" + "gene" modes:
#   - level = "pairwise": model states on per-ortholog-pair rows (your modernized version)
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
    ok <- try(eval(parse(text = filter_expr), envir = env), silent = TRUE)
    if (!inherits(ok, "try-error")) {
      ok <- as.logical(ok)
      if (length(ok) == 1L) ok <- rep(ok, nrow(d))
      if (length(ok) != nrow(d)) {
        stop("filter_expr must return a logical vector of length 1 or nrow(d).")
      }
      keep <- keep & ok
    }
  }

  d[keep, , drop = FALSE]
}

#' Read and normalize a regions BED-like file
#'
#' Returns a data.frame with columns: seqname, start, end, region_name
#'
#' Notes:
#' - Tries to infer whether the file has a header. For best results, provide a header row.
#' - Column selectors can be column names (character) or 1-based column indices (numeric).
#'
#' @keywords internal
.read_regions_bed <- function(regions_bed,
                              region_seq_col   = NULL,
                              region_start_col = NULL,
                              region_end_col   = NULL,
                              region_name_col  = NULL) {
  if (is.null(regions_bed)) return(NULL)
  if (!file.exists(regions_bed)) stop("regions_bed not found: ", regions_bed)

  # Try header=TRUE first; if start/end look non-numeric, fall back to header=FALSE
  reg1 <- try(.read_ws(regions_bed, header_try = TRUE), silent = TRUE)
  if (inherits(reg1, "try-error")) {
    reg_raw <- .read_ws(regions_bed, header_try = FALSE)
    names(reg_raw) <- paste0("V", seq_len(ncol(reg_raw)))
  } else {
    reg_raw <- reg1
    if (ncol(reg_raw) >= 3) {
      s2 <- suppressWarnings(as.numeric(reg_raw[[2]]))
      s3 <- suppressWarnings(as.numeric(reg_raw[[3]]))
      frac_bad <- mean(is.na(s2) | is.na(s3))
      if (!is.finite(frac_bad)) frac_bad <- 1
      if (frac_bad > 0.5) {
        reg_raw <- .read_ws(regions_bed, header_try = FALSE)
        names(reg_raw) <- paste0("V", seq_len(ncol(reg_raw)))
      }
    }
  }

  if (ncol(reg_raw) < 3) stop("regions_bed must have at least 3 columns: seqname, start, end.")

  # Resolve column selectors (character name or numeric index)
  .col_get <- function(df, sel, default_idx) {
    if (is.null(sel)) return(df[[default_idx]])
    if (is.numeric(sel)) {
      sel <- as.integer(sel[1])
      if (sel < 1L || sel > ncol(df)) stop("regions_bed column index out of range: ", sel)
      return(df[[sel]])
    }
    sel <- as.character(sel[1])
    if (!sel %in% names(df)) stop("regions_bed column not found: ", sel)
    df[[sel]]
  }

  seqv   <- .col_get(reg_raw, region_seq_col,   1)
  startv <- .col_get(reg_raw, region_start_col, 2)
  endv   <- .col_get(reg_raw, region_end_col,   3)

  # region_name optional
  if (is.null(region_name_col)) {
    namev <- rep("region", nrow(reg_raw))
  } else if (is.numeric(region_name_col)) {
    idx <- as.integer(region_name_col[1])
    if (idx < 1L || idx > ncol(reg_raw)) {
      namev <- rep("region", nrow(reg_raw))
    } else {
      namev <- reg_raw[[idx]]
    }
  } else {
    nm <- as.character(region_name_col[1])
    if (!nm %in% names(reg_raw)) {
      namev <- rep("region", nrow(reg_raw))
    } else {
      namev <- reg_raw[[nm]]
    }
  }

  data.frame(
    seqname     = as.character(seqv),
    start       = as.numeric(startv),
    end         = as.numeric(endv),
    region_name = as.character(namev),
    stringsAsFactors = FALSE
  )
}

#' Label rows as inside/outside regions given explicit coord columns
#'
#' Adds region_status ("region"/"background") based on overlap.
#' Uses data.table::foverlaps if available; otherwise falls back to a loop.
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

  has_dt <- requireNamespace("data.table", quietly = TRUE)
  lab <- rep("background", nrow(d))

  if (has_dt) {
    dt_g <- data.table::data.table(
      idx     = seq_len(nrow(d)),
      seqname = as.character(d[[seq_col]]),
      start   = as.numeric(d[[start_col]]),
      end     = as.numeric(d[[end_col]])
    )
    dt_r <- data.table::data.table(
      seqname = regions$seqname,
      start   = regions$start,
      end     = regions$end
    )
    data.table::setkey(dt_g, seqname, start, end)
    data.table::setkey(dt_r, seqname, start, end)
    ov <- data.table::foverlaps(dt_g, dt_r, nomatch = 0L)
    in_idx <- unique(ov$idx)
    lab[in_idx] <- "region"
  } else {
    for (i in seq_len(nrow(regions))) {
      r <- regions[i, ]
      hits <- which(d[[seq_col]] == r$seqname &
                      d[[start_col]] <= r$end &
                      d[[end_col]]   >= r$start)
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
#' Uses both neutral bounds:
#' - Purifying: x < neutral_lower
#' - Neutral:   neutral_lower <= x <= neutral_upper
#' - Positive:  x > neutral_upper
#'
#' @keywords internal
.make_states <- function(d, dnds_col, pos_threshold, neutral_lower, neutral_upper) {
  x <- d[[dnds_col]]
  d$state <- ifelse(x > neutral_upper, "Positive",
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

  # Separation guard: each group must have both outcomes
  tt <- table(d2[[group_col]], d2$y)
  if (ncol(tt) < 2L) return(NULL)
  if (any(tt[, "0"] == 0L | tt[, "1"] == 0L)) return(NULL)

  d2[[group_col]] <- stats::relevel(factor(d2[[group_col]]), ref = groupB)

  use_re <- !is.null(random_effect_col) &&
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

      # pick the single non-intercept term (since only A vs B present)
      tt2 <- setdiff(names(coefs), "(Intercept)")
      if (length(tt2) != 1) return(NULL)
      term <- tt2[1]

      est <- unname(coefs[term])
      se  <- sqrt(diag(vcovm))[term]
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

  # with only A/B present, there should be one coefficient
  if (length(rn) != 1) return(NULL)
  rr <- rn[1]

  est <- sm[rr, "Estimate"]
  se  <- sm[rr, "Std. Error"]
  pval <- sm[rr, "Pr(>|z|)"]

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

  # ordering for readability within facets
  res$contrast_label <- factor(res$contrast_label, levels = rev(unique(res$contrast_label)))

  gg <- ggplot2::ggplot(res, ggplot2::aes(x = odds_ratio, y = contrast_label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = or_ci_lower, xmax = or_ci_upper),
      height = 0.2,
      linewidth = line_width
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_x_log10() +
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

#' Infer subgenome/group from seqname (default: capture terminal B/C/D)
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
#' If random_effect_col is provided and exists in the pairwise table, it is propagated
#' onto gene rows so gene-mode mixed models can be fit.
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

  if (group_mode == "subgenome") {
    g$group <- .infer_group_from_seqname(g$seqname)
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

  # attach a random-effect column if requested and present in gene rows
  if (!is.null(random_effect_col) && nzchar(random_effect_col) && random_effect_col %in% names(g)) {
    re_map <- g[, c("gene_key", random_effect_col), drop = FALSE]
    re_map <- re_map[!is.na(re_map[[random_effect_col]]), , drop = FALSE]
    re_map <- re_map[!duplicated(re_map$gene_key), , drop = FALSE]
    agg <- merge(agg, re_map, by = "gene_key", all.x = TRUE)
  }

  agg <- agg[is.finite(agg$dNdS_agg) & agg$n_pairs >= min_pairs, , drop = FALSE]
  agg
}

# -----------------------------
# Exported function
# -----------------------------

#' dN/dS state contrasts (pairwise or gene-level) via logistic regression + forest plot
#'
#' Pairwise mode:
#'   - input is a table containing dNdS + group column (or batch reads <comp>_dnds_annot.tsv),
#'   - optionally restrict to regions (regions_bed + side),
#'   - fit pairwise A vs B logistic models for each state within each scope (global/region).
#'
#' Gene mode:
#'   - input is one or more *_dnds_annot.tsv tables (single or batch),
#'   - convert to gene rows on query/subject side(s),
#'   - optionally restrict to regions (regions_bed) before aggregation,
#'   - aggregate across partners per gene (median/mean),
#'   - classify each gene into state, and fit pairwise A vs B models on gene-level observations.
#'
#' Region coordinates are treated as inclusive numeric intervals; ensure your `regions_bed`
#' coordinates are consistent with the coordinate system used in your dN/dS tables.
#'
#' @param level One of c("pairwise","gene").
#'
#' @param dnds_table_file Pairwise mode single input TSV. In gene mode, this can also be
#'   used as a convenience alias for a single *_dnds_annot.tsv file.
#' @param comparison_file Batch mode comparisons file (both levels).
#' @param output_dir Root directory for batch mode.
#' @param in_suffix Filename within each comparison folder:
#'   - pairwise: defaults to "<comp>_dnds_annot.tsv"
#'   - gene:     defaults to "<comp>_dnds_annot.tsv"
#'
#' @param dnds_col Name of dN/dS column (default "dNdS").
#' @param group_col Pairwise mode grouping column name in the table (default "group").
#' @param random_effect_col Optional column for random intercept (e.g., orthogroup) (default NULL).
#'
#' @param pos_threshold,neutral_lower,neutral_upper State thresholds.
#' @param max_dnds Drop rows with dN/dS >= max_dnds or NA.
#' @param filter_expr Optional expression evaluated in raw tables before modeling.
#'
#' @param min_n_per_group Minimum observations per group within a scope.
#'
#' @param fdr_method One of "BH","BY","none". Applied within each (state, scope).
#'
#' @param regions_bed Optional BED-like file. If provided, enables region scope.
#' @param region_seq_col,region_start_col,region_end_col Column selectors in `regions_bed`
#'   giving the chromosome/seqname, start, and end columns. Each selector may be a column
#'   name (character) or 1-based column index (numeric). If NULL (default), the first three
#'   columns of the file are used.
#' @param region_name_col Optional column selector in `regions_bed` for labeling regions
#'   (name or 1-based index). If NULL or not present, a constant label ("region") is used.
#' @param side Pairwise mode region overlap side ("query" or "subject"). Ignored in gene mode
#'   (gene mode uses the gene-row coordinates directly).
#' @param scopes Which scopes to run: subset of c("global","region").
#'
#' @param make_plots If TRUE, write forest plot PDF and PNG.
#' @param base_family Font family for plots.
#' @param out_prefix Output prefix (default derived from file/comp name).
#'
#' ---- gene-mode only ----
#' @param dnds_annot_files Optional vector of explicit *_dnds_annot.tsv paths (gene mode single).
#' @param focal_sides Which genes to treat as focal observations: c("query","subject") (gene mode; default both).
#' @param group_mode "subgenome" infers group from seqname suffix (B/C/D) or "custom".
#' @param gene_group_col If group_mode="custom", name of column in gene rows to use as group
#'   (rare; usually you'd keep group_mode="subgenome"). The column must exist in the gene-row
#'   table prior to aggregation.
#' @param agg_fun "median" or "mean" aggregation across partners.
#' @param min_pairs Minimum pairwise observations per gene required before modeling.
#' @param gene_out_dir Subdirectory (within output_dir) for gene-mode outputs in batch runs.
#'
#' @return Invisibly returns list with:
#'   - results: data.frame (or list of per-comp dfs in batch)
#'   - paths: written TSV + plot paths
#'   - gene_table: (gene mode only) the aggregated gene table (or list in batch)
#'
#' @export
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
                                regions_bed       = NULL,
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

  regions <- .read_regions_bed(regions_bed,
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
      stop("scopes includes 'region' but regions_bed is NULL.")
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

    # filter
    d <- .filter_dnds(d, dnds_col = dnds_col, filter_expr = filter_expr, max_dnds = max_dnds)
    if (!nrow(d)) {
      warning("[dnds_state_contrast/pairwise] No rows after filtering for: ", label)
      return(list(results = NULL, paths = character(0)))
    }

    # region status (pairwise uses side switch)
    d <- .label_region_status_pairwise(d, regions = regions, side = side)

    # states
    d <- .make_states(d, dnds_col = dnds_col,
                      pos_threshold = pos_threshold,
                      neutral_lower = neutral_lower,
                      neutral_upper = neutral_upper)

    all_res <- list()

    for (scp in scopes) {
      d_sc <- if (scp == "global") d else d[d$region_status == "region", , drop = FALSE]

      # power filter per group
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

    # FDR within each (state, scope)
    res$p_adj <- NA_real_
    key <- paste(res$state, res$scope, sep = "||")
    for (k in unique(key)) {
      idx <- which(key == k)
      res$p_adj[idx] <- .adjust_p(res$p_value[idx], method = fdr_method)
    }

    # OR columns
    res$odds_ratio  <- exp(res$log_or)
    res$or_ci_lower <- exp(res$log_or_ci_lower)
    res$or_ci_upper <- exp(res$log_or_ci_upper)

    # context columns
    res$pos_threshold <- pos_threshold
    res$neutral_lower <- neutral_lower
    res$neutral_upper <- neutral_upper
    res$max_dnds      <- max_dnds
    res$min_n_per_group <- min_n_per_group
    res$fdr_method    <- fdr_method

    # write outputs
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
    # stack gene rows from all files
    gene_rows <- list()

    for (i in seq_along(files)) {
      f <- files[i]
      if (!file.exists(f)) next
      comp_name <- sub("_dnds_annot\\.tsv$", "", basename(f))

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

      # label region_status per gene row (uses the gene coords directly)
      g <- .label_region_status_coords(g, regions = regions, seq_col = "seqname", start_col = "start", end_col = "end")

      gene_rows[[length(gene_rows) + 1]] <- g
    }

    if (!length(gene_rows)) {
      warning("[dnds_state_contrast/gene] No usable rows after filtering for: ", label)
      return(list(results = NULL, paths = character(0), gene_table = NULL))
    }

    g_all <- do.call(rbind, gene_rows)

    # compute results per scope
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

      # enforce power at gene level
      tab <- table(gene_tbl$group)
      keep_groups <- names(tab)[tab >= min_n_per_group]
      gene_tbl <- gene_tbl[gene_tbl$group %in% keep_groups, , drop = FALSE]
      groups <- unique(as.character(gene_tbl$group))
      if (length(groups) < 2) {
        warning("[dnds_state_contrast/gene] Not enough groups after min_n_per_group in scope='", scp, "' for: ", label)
        next
      }

      # add state on aggregated dNdS
      gene_tbl <- .make_states(gene_tbl, dnds_col = "dNdS_agg",
                               pos_threshold = pos_threshold,
                               neutral_lower = neutral_lower,
                               neutral_upper = neutral_upper)

      gene_tbl$scope <- scp
      gene_tbl_by_scope[[scp]] <- gene_tbl

      pairs <- utils::combn(sort(groups), 2, simplify = FALSE)
      states <- levels(gene_tbl$state)

      for (st in states) {
        gene_tbl$y <- as.integer(gene_tbl$state == st)
        for (p in pairs) {
          A <- p[1]; B <- p[2]
          fit <- .fit_pair(gene_tbl, group_col = "group", groupA = A, groupB = B,
                           random_effect_col = random_effect_col)
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

    # FDR within each (state, scope)
    res$p_adj <- NA_real_
    key <- paste(res$state, res$scope, sep = "||")
    for (k in unique(key)) {
      idx <- which(key == k)
      res$p_adj[idx] <- .adjust_p(res$p_value[idx], method = fdr_method)
    }

    # OR columns + context
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

    # write outputs
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    prefix <- if (!is.null(out_prefix) && nzchar(out_prefix)) out_prefix else label

    out_tsv <- file.path(out_dir, sprintf("%s_dnds_state_contrasts.tsv", prefix))
    utils::write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    out_paths <- c(out_tsv)

    # also write gene tables (per scope) so you can sanity check
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
        # gene mode: in batch we run per-comp by default (single file), but write into a shared subdir
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
    # allow dnds_table_file as a convenience alias in gene mode too
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
