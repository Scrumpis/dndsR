# contrast.R
# (updated: global-by-default; optional regional restriction via regions_bed)

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

#' Read a contrast file defining pairwise dN/dS comparisons
#'
#' Supports either:
#'   4 columns: contrast_name, compA, compB, side
#'     (side is applied to both compA and compB)
#'   or
#'   5 columns: contrast_name, compA, compA_side, compB, compB_side
#'
#' Returns a data.frame with normalized columns:
#'   contrast_name, compA, compA_side, compB, compB_side
#'
#' Note: when regions_bed is NULL (global mode), these side columns are accepted
#' but ignored.
#'
#' @keywords internal
.read_contrasts <- function(x) {
  req5 <- c("contrast_name", "compA", "compA_side", "compB", "compB_side")
  req4 <- c("contrast_name", "compA", "compB", "side")

  normalize_df <- function(df) {
    if (all(req5 %in% names(df))) return(df[, req5, drop = FALSE])
    if (all(req4 %in% names(df))) {
      return(data.frame(
        contrast_name = df$contrast_name,
        compA         = df$compA,
        compA_side    = df$side,
        compB         = df$compB,
        compB_side    = df$side,
        stringsAsFactors = FALSE
      ))
    }
    stop("contrast_file must have either 4 cols (",
         paste(req4, collapse = ", "),
         ") or 5 cols (",
         paste(req5, collapse = ", "),
         "), or a header with those names.")
  }

  if (is.data.frame(x)) return(normalize_df(x))

  df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
  if (!inherits(df1, "try-error")) {
    if (all(req5 %in% names(df1)) || all(req4 %in% names(df1))) {
      return(normalize_df(df1))
    }
  }

  df2 <- .read_ws(x, header_try = FALSE)
  if (ncol(df2) == 5) {
    names(df2)[1:5] <- req5
  } else if (ncol(df2) == 4) {
    names(df2)[1:4] <- req4
  } else {
    stop("contrast_file must have 4 or 5 columns if no header is present.")
  }
  normalize_df(df2)
}

#' Filter dNdS annotation table by NA / max_dnds / optional logical expression
#'
#' @keywords internal
.filter_dnds <- function(d, filter_expr = NULL, max_dnds = 10) {
  keep <- !is.na(d$dNdS) & d$dNdS < max_dnds
  if (!is.null(filter_expr) && nzchar(filter_expr)) {
    ok <- try(eval(parse(text = filter_expr),
                   envir  = d,
                   enclos = parent.frame()),
              silent = TRUE)
    if (!inherits(ok, "try-error")) {
      keep <- keep & isTRUE(as.vector(ok))
    }
  }
  d[keep, , drop = FALSE]
}

#' Read and normalize a regions BED-like file
#'
#' Returns a data.frame with columns: seqname, start, end, region_name
#'
#' @keywords internal
.read_regions_bed <- function(regions_bed,
                              region_seq_col   = NULL,
                              region_start_col = NULL,
                              region_end_col   = NULL,
                              region_name_col  = NULL) {
  if (missing(regions_bed) || is.null(regions_bed) || !nzchar(regions_bed)) {
    stop("regions_bed is required.")
  }
  if (!file.exists(regions_bed)) stop("regions_bed not found: ", regions_bed)

  reg_raw <- .read_ws(regions_bed, header_try = TRUE)
  if (ncol(reg_raw) < 3) stop("regions_bed must have at least 3 columns: seqname, start, end.")

  if (is.null(region_seq_col))   region_seq_col   <- names(reg_raw)[1]
  if (is.null(region_start_col)) region_start_col <- names(reg_raw)[2]
  if (is.null(region_end_col))   region_end_col   <- names(reg_raw)[3]
  if (is.null(region_name_col) || !region_name_col %in% names(reg_raw)) {
    reg_raw$region_name <- "region"
    region_name_col <- "region_name"
  }

  data.frame(
    seqname     = as.character(reg_raw[[region_seq_col]]),
    start       = as.numeric(reg_raw[[region_start_col]]),
    end         = as.numeric(reg_raw[[region_end_col]]),
    region_name = as.character(reg_raw[[region_name_col]]),
    stringsAsFactors = FALSE
  )
}

#' Label dNdS rows as inside/outside regions for a given side
#'
#' Adds a column q_region_status or s_region_status with values "region"/"background".
#' Uses data.table::foverlaps if available; otherwise falls back to a simple loop.
#'
#' @keywords internal
.label_region_side <- function(d, regions, side = c("query", "subject")) {
  side   <- match.arg(side)
  has_dt <- requireNamespace("data.table", quietly = TRUE)

  if (side == "query") {
    seq_col    <- "q_gff_seqname"
    start_col  <- "q_gff_start"
    end_col    <- "q_gff_end"
    status_col <- "q_region_status"
  } else {
    seq_col    <- "s_gff_seqname"
    start_col  <- "s_gff_start"
    end_col    <- "s_gff_end"
    status_col <- "s_region_status"
  }

  if (!all(c(seq_col, start_col, end_col) %in% names(d))) {
    stop("Missing required columns for side='", side, "': ",
         paste(c(seq_col, start_col, end_col), collapse = ", "))
  }

  if (has_dt) {
    dt_g <- data.table::data.table(
      idx     = seq_len(nrow(d)),
      seqname = as.character(d[[seq_col]]),
      start   = as.numeric(d[[start_col]]),
      end     = as.numeric(d[[end_col]])
    )
    dt_r <- data.table::data.table(
      seqname     = regions$seqname,
      start       = regions$start,
      end         = regions$end,
      region_name = regions$region_name
    )
    data.table::setkey(dt_g, seqname, start, end)
    data.table::setkey(dt_r, seqname, start, end)
    ov     <- data.table::foverlaps(dt_g, dt_r, nomatch = 0L)
    in_idx <- unique(ov$idx)
    lab <- rep("background", nrow(d))
    lab[in_idx] <- "region"
    d[[status_col]] <- lab
  } else {
    lab <- rep("background", nrow(d))
    for (i in seq_len(nrow(regions))) {
      r <- regions[i, ]
      hits <- which(d[[seq_col]] == r$seqname &
                      d[[start_col]] <= r$end &
                      d[[end_col]]   >= r$start)
      if (length(hits)) lab[hits] <- "region"
    }
    d[[status_col]] <- lab
  }
  d
}

#' Compute mean/median/SD/SE/CI for a numeric vector
#'
#' @keywords internal
.calc_mean_ci <- function(x, ci_method = c("normal", "bootstrap"), n_boot = 1000) {
  ci_method <- match.arg(ci_method)
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_))
  m   <- mean(x)
  med <- stats::median(x)
  sd  <- stats::sd(x)
  se  <- sd / sqrt(max(1L, n))
  if (ci_method == "normal") {
    ci <- m + c(-1.96, 1.96) * se
  } else {
    means <- replicate(n_boot, mean(sample(x, size = n, replace = TRUE)))
    ci <- stats::quantile(means, probs = c(0.025, 0.975), names = FALSE)
  }
  c(m, med, sd, se, ci[1], ci[2])
}

#' Clean column names
#'
#' Strips CRLF artifacts and trims leading/trailing whitespace from column names.
#'
#' @param d A data.frame-like object.
#' @return `d` with cleaned column names.
#' @keywords internal
.clean_colnames <- function(d) {
  nn <- names(d)
  nn <- sub("\r$", "", nn)
  nn <- trimws(nn)
  names(d) <- nn
  d
}

# -----------------------------
# dN/dS contrasts function
# -----------------------------

#' Pairwise dN/dS contrasts (genome-wide by default; optional region restriction)
#'
#' Compare dN/dS between pairs of dNdS annotation files using paired nonparametric
#' tests. For each contrast (compA vs compB), the function:
#'   - filters each table by max_dnds and optional filter_expr,
#'   - (optional) restricts to genes in regions_bed (side-aware),
#'   - matches rows across comparisons,
#'   - computes delta = dNdS_A - dNdS_B per matched gene,
#'   - runs a Wilcoxon signed-rank test on delta,
#'   - runs Fisher's exact test for enrichment of dN/dS > pos_threshold in A vs B,
#'   - writes a TSV summary + optional plots.
#'
#' Default behavior is genome-wide contrasts (regions_bed = NULL).
#'
#' Matching behavior:
#'   - Global (regions_bed = NULL): merge on c("query_id","subject_id") by default
#'     (recommended). You may override with merge_cols.
#'   - Regional (regions_bed provided): merge on side-aware focal IDs by default
#'     (query_id if side == "query", else subject_id). You may override with merge_cols.
#'
#' Modes:
#'   1) Single-contrast: supply dnds_annot_file_a and dnds_annot_file_b.
#'   2) Batch + explicit contrasts: comparison_file + contrast_file + output_dir.
#'   3) Batch + auto-all-pairs: comparison_file + contrast_file = NULL.
#'
#' @param dnds_annot_file_a Path to first <compA>_dnds_annot.tsv (single-contrast).
#' @param dnds_annot_file_b Path to second <compB>_dnds_annot.tsv (single-contrast).
#'
#' @param comparison_file Path to whitespace-delimited file with columns:
#'   comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#' @param contrast_file Optional file defining contrasts (4-col or 5-col layout).
#'   Note: side columns are ignored in global mode (regions_bed = NULL).
#'
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#'
#' @param regions_bed Optional BED-like file of regions. If provided, restrict to
#'   genes overlapping regions (regional mode). If NULL, run genome-wide (default).
#' @param region_seq_col,region_start_col,region_end_col,region_name_col
#'   Column names in regions_bed for seqname/start/end/label.
#'
#' @param merge_cols Optional character vector of column names used to match rows
#'   across comparisons. If NULL, defaults are used (see Matching behavior above).
#'
#' @param sides Character vector among c("query","subject") indicating which side(s)
#'   to evaluate in auto-all-pairs mode and single-contrast mode. Ignored in global
#'   mode (regions_bed = NULL).
#'
#' @param filter_expr Optional character with a logical expression evaluated in
#'   each dNdS table prior to region labeling and merging.
#' @param max_dnds Numeric. Drop rows with dNdS >= max_dnds or NA dNdS (default 10).
#'
#' @param ci_method One of "normal" (mean +/- 1.96*SE) or "bootstrap" (percentile).
#' @param n_boot Integer; number of bootstrap resamples if ci_method = "bootstrap".
#'
#' @param make_plots Logical; if TRUE, write paired scatter, delta histogram,
#'   and delta ECDF per contrast.
#'
#' @param pos_threshold Numeric threshold for calling "positive selection"
#'   (genes with dN/dS > pos_threshold are counted in the enrichment test; default 1).
#'
#' @return Invisibly, a character vector of summary TSV paths (one per contrast x side_tag).
#' @export
contrast <- function(dnds_annot_file_a = NULL,
                           dnds_annot_file_b = NULL,
                           comparison_file   = NULL,
                           contrast_file     = NULL,
                           output_dir        = getwd(),
                           regions_bed       = NULL,
                           region_seq_col    = NULL,
                           region_start_col  = NULL,
                           region_end_col    = NULL,
                           region_name_col   = NULL,
                           merge_cols        = NULL,
                           sides             = c("query", "subject"),
                           filter_expr       = NULL,
                           max_dnds          = 10,
                           ci_method         = c("normal", "bootstrap"),
                           n_boot            = 1000,
                           make_plots        = TRUE,
                           pos_threshold     = 1) {

  ci_method <- match.arg(ci_method)
  sides     <- match.arg(sides, choices = c("query", "subject"), several.ok = TRUE)

  do_regions <- !is.null(regions_bed) && nzchar(regions_bed)

  regions <- NULL
  if (do_regions) {
    regions <- .read_regions_bed(regions_bed,
                                 region_seq_col   = region_seq_col,
                                 region_start_col = region_start_col,
                                 region_end_col   = region_end_col,
                                 region_name_col  = region_name_col)
  }

  .wilcox_signed <- function(delta) {
    delta <- delta[is.finite(delta) & delta != 0]
    if (length(delta) < 2L) return(NA_real_)
    stats::wilcox.test(delta, mu = 0, alternative = "two.sided")$p.value
  }

  .plot_contrast <- function(merged, contrast_name, compA_name, compB_name, out_dir, side_tag, mode_label) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
    if (!nrow(merged)) return(invisible(NULL))

    # Scatter: B vs A (log1p both)
    gg <- ggplot2::ggplot(merged, ggplot2::aes(x = dNdS_B, y = dNdS_A)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      ggplot2::geom_point(alpha = 0.4, size = 1) +
      ggplot2::scale_x_continuous(trans = "log1p") +
      ggplot2::scale_y_continuous(trans = "log1p") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = paste0(contrast_name, " (", mode_label, "; ", side_tag, ")"),
        x     = paste0("dN/dS (", compB_name, ")"),
        y     = paste0("dN/dS (", compA_name, ")")
      )

    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(contrast_name, "_", side_tag, "_dnds_scatter.pdf")),
      plot     = gg,
      width    = 5,
      height   = 4
    )

    # Delta histogram
    gg2 <- ggplot2::ggplot(merged, ggplot2::aes(x = delta)) +
      ggplot2::geom_histogram(bins = 40, alpha = 0.7) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        x     = "delta dN/dS (A - B)",
        y     = "Count",
        title = paste0(contrast_name, " — delta histogram (", side_tag, ")")
      )

    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(contrast_name, "_", side_tag, "_delta_hist.pdf")),
      plot     = gg2,
      width    = 5,
      height   = 4
    )

    # Delta ECDF
    gg3 <- ggplot2::ggplot(merged, ggplot2::aes(x = delta)) +
      ggplot2::stat_ecdf() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        x     = "delta dN/dS (A - B)",
        y     = "ECDF",
        title = paste0(contrast_name, " — delta ECDF (", side_tag, ")")
      )

    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(contrast_name, "_", side_tag, "_delta_ecdf.pdf")),
      plot     = gg3,
      width    = 5,
      height   = 4
    )

    invisible(NULL)
  }

  .require_cols <- function(d, cols, context) {
    miss <- setdiff(cols, names(d))
    if (length(miss)) stop("Missing required column(s) in ", context, ": ", paste(miss, collapse = ", "))
    invisible(TRUE)
  }

  .merge_dnds <- function(dA, dB, by_cols, suffixes = c("_A","_B")) {
    # drop empty/NA keys
    for (cc in by_cols) {
      dA <- dA[!is.na(dA[[cc]]) & nzchar(as.character(dA[[cc]])), , drop = FALSE]
      dB <- dB[!is.na(dB[[cc]]) & nzchar(as.character(dB[[cc]])), , drop = FALSE]
    }

    # deduplicate by key (keep first)
    keyA <- do.call(paste, c(dA[by_cols], sep = "\r"))
    keyB <- do.call(paste, c(dB[by_cols], sep = "\r"))
    dA <- dA[!duplicated(keyA), , drop = FALSE]
    dB <- dB[!duplicated(keyB), , drop = FALSE]

    merge(dA, dB, by = by_cols, suffixes = suffixes)
  }

  .run_one_contrast <- function(contrast_name,
                                compA_name,
                                compB_name,
                                sideA,
                                sideB,
                                fileA,
                                fileB,
                                out_dir) {

    if (!file.exists(fileA)) stop("dNdS file A not found: ", fileA)
    if (!file.exists(fileB)) stop("dNdS file B not found: ", fileB)

    dA_raw <- utils::read.table(fileA, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                               quote = "", comment.char = "", check.names = FALSE)
    dA_raw <- .clean_colnames(dA_raw)

    dB_raw <- utils::read.table(fileB, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                               quote = "", comment.char = "", check.names = FALSE)
    dB_raw <- .clean_colnames(dB_raw)

    if (!nrow(dA_raw) || !nrow(dB_raw)) {
      warning("Empty dNdS table(s) for contrast ", contrast_name)
      return(NA_character_)
    }

    .require_cols(dA_raw, c("dNdS"), context = fileA)
    .require_cols(dB_raw, c("dNdS"), context = fileB)

    # Apply filtering
    dA <- .filter_dnds(dA_raw, filter_expr = filter_expr, max_dnds = max_dnds)
    dB <- .filter_dnds(dB_raw, filter_expr = filter_expr, max_dnds = max_dnds)

    if (!nrow(dA) || !nrow(dB)) {
      warning("No rows after filtering for contrast ", contrast_name)
      return(NA_character_)
    }

    mode_label <- if (do_regions) "regional (region-only)" else "genome-wide"

    # ----------------------------
    # Regional restriction (optional)
    # ----------------------------
    if (do_regions) {
      if (!sideA %in% c("query","subject")) stop("sideA must be 'query' or 'subject' in regional mode.")
      if (!sideB %in% c("query","subject")) stop("sideB must be 'query' or 'subject' in regional mode.")

      dA <- .label_region_side(dA, regions = regions, side = sideA)
      dB <- .label_region_side(dB, regions = regions, side = sideB)

      status_col_A <- if (sideA == "query") "q_region_status" else "s_region_status"
      status_col_B <- if (sideB == "query") "q_region_status" else "s_region_status"

      dA <- dA[dA[[status_col_A]] == "region", , drop = FALSE]
      dB <- dB[dB[[status_col_B]] == "region", , drop = FALSE]

      if (!nrow(dA) || !nrow(dB)) {
        warning("No region rows for contrast ", contrast_name)
        return(NA_character_)
      }
    }

    # ----------------------------
    # Matching / merge
    # ----------------------------
    if (is.null(merge_cols)) {
      if (!do_regions) {
        # Global default: ortholog-pair key
        merge_cols_use <- c("query_id","subject_id")
        .require_cols(dA, merge_cols_use, context = paste0("table A (", compA_name, ")"))
        .require_cols(dB, merge_cols_use, context = paste0("table B (", compB_name, ")"))

        dA_key <- data.frame(
          query_id   = as.character(dA$query_id),
          subject_id = as.character(dA$subject_id),
          dNdS       = as.numeric(dA$dNdS),
          stringsAsFactors = FALSE
        )
        dB_key <- data.frame(
          query_id   = as.character(dB$query_id),
          subject_id = as.character(dB$subject_id),
          dNdS       = as.numeric(dB$dNdS),
          stringsAsFactors = FALSE
        )

        merged <- .merge_dnds(dA_key, dB_key, by_cols = merge_cols_use)
        side_tag <- "global"

      } else {
        # Regional default: side-aware focal matching
        .require_cols(dA, c("query_id","subject_id"), context = paste0("table A (", compA_name, ")"))
        .require_cols(dB, c("query_id","subject_id"), context = paste0("table B (", compB_name, ")"))

        focal_id_A <- if (sideA == "query") dA$query_id else dA$subject_id
        focal_id_B <- if (sideB == "query") dB$query_id else dB$subject_id

        dA_key <- data.frame(
          focal_id = as.character(focal_id_A),
          dNdS     = as.numeric(dA$dNdS),
          stringsAsFactors = FALSE
        )
        dB_key <- data.frame(
          focal_id = as.character(focal_id_B),
          dNdS     = as.numeric(dB$dNdS),
          stringsAsFactors = FALSE
        )

        merged <- .merge_dnds(dA_key, dB_key, by_cols = "focal_id")
        side_tag <- paste0("A_", sideA, "__B_", sideB)
      }

    } else {
      # User-specified merge columns (applies in either mode)
      merge_cols_use <- as.character(merge_cols)
      .require_cols(dA, merge_cols_use, context = paste0("table A (", compA_name, ")"))
      .require_cols(dB, merge_cols_use, context = paste0("table B (", compB_name, ")"))

      # keep only merge cols + dNdS
      dA_key <- dA[, unique(c(merge_cols_use, "dNdS")), drop = FALSE]
      dB_key <- dB[, unique(c(merge_cols_use, "dNdS")), drop = FALSE]

      # coerce merge cols to character for consistency
      for (cc in merge_cols_use) {
        dA_key[[cc]] <- as.character(dA_key[[cc]])
        dB_key[[cc]] <- as.character(dB_key[[cc]])
      }
      dA_key$dNdS <- as.numeric(dA_key$dNdS)
      dB_key$dNdS <- as.numeric(dB_key$dNdS)

      merged <- .merge_dnds(dA_key, dB_key, by_cols = merge_cols_use)
      side_tag <- if (do_regions) paste0("custommerge__A_", sideA, "__B_", sideB) else "custommerge__global"
    }

    if (!nrow(merged)) {
      warning("No overlapping matched rows for contrast ", contrast_name, " (", mode_label, ")")
      return(NA_character_)
    }

    # standardize names for downstream code
    if (!("dNdS_A" %in% names(merged))) {
      # merge() uses suffixes; ensure exactly dNdS_A / dNdS_B exist
      # If user supplied merge_cols, dNdS columns should have been suffixed.
      # Defensive fallback:
      dcols <- grep("^dNdS", names(merged), value = TRUE)
      if (length(dcols) == 2) {
        names(merged)[match(dcols[1], names(merged))] <- "dNdS_A"
        names(merged)[match(dcols[2], names(merged))] <- "dNdS_B"
      }
    }

    if (!all(c("dNdS_A","dNdS_B") %in% names(merged))) {
      stop("Internal error: merged table is missing dNdS_A and/or dNdS_B after merging.")
    }

    merged$delta <- merged$dNdS_A - merged$dNdS_B

    # ----------------------------
    # Stats / tests
    # ----------------------------
    stats_delta  <- .calc_mean_ci(merged$delta, ci_method = ci_method, n_boot = n_boot)
    mean_delta   <- stats_delta[1]
    median_delta <- stats_delta[2]
    delta_sd     <- stats_delta[3]
    delta_se     <- stats_delta[4]
    delta_ci_lo  <- stats_delta[5]
    delta_ci_hi  <- stats_delta[6]

    mean_A <- mean(merged$dNdS_A, na.rm = TRUE)
    med_A  <- stats::median(merged$dNdS_A, na.rm = TRUE)
    mean_B <- mean(merged$dNdS_B, na.rm = TRUE)
    med_B  <- stats::median(merged$dNdS_B, na.rm = TRUE)
    n      <- sum(is.finite(merged$delta))

    p_wilcox <- .wilcox_signed(merged$delta)

    # Positive selection enrichment: dN/dS > pos_threshold
    is_finite_A <- is.finite(merged$dNdS_A)
    is_finite_B <- is.finite(merged$dNdS_B)
    is_pos_A    <- is_finite_A & merged$dNdS_A > pos_threshold
    is_pos_B    <- is_finite_B & merged$dNdS_B > pos_threshold

    nA         <- sum(is_finite_A)
    nB         <- sum(is_finite_B)
    n_pos_A    <- sum(is_pos_A)
    n_pos_B    <- sum(is_pos_B)
    n_non_A    <- nA - n_pos_A
    n_non_B    <- nB - n_pos_B
    frac_pos_A <- if (nA > 0) n_pos_A / nA else NA_real_
    frac_pos_B <- if (nB > 0) n_pos_B / nB else NA_real_

    fisher_p <- NA_real_
    if (nA > 0 && nB > 0) {
      tab <- matrix(c(n_pos_A, n_non_A,
                      n_pos_B, n_non_B),
                    nrow = 2, byrow = TRUE)
      fisher_p <- tryCatch(
        stats::fisher.test(tab, alternative = "greater")$p.value,
        error = function(e) NA_real_
      )
    }

    summary_df <- data.frame(
      contrast_name   = contrast_name,
      compA           = compA_name,
      compB           = compB_name,
      mode            = if (do_regions) "regional" else "global",
      sideA           = if (do_regions) sideA else "global",
      sideB           = if (do_regions) sideB else "global",
      side_tag        = side_tag,
      n               = n,
      max_dnds        = max_dnds,
      mean_dnds_A     = mean_A,
      median_dnds_A   = med_A,
      mean_dnds_B     = mean_B,
      median_dnds_B   = med_B,
      mean_delta      = mean_delta,
      median_delta    = median_delta,
      delta_sd        = delta_sd,
      delta_se        = delta_se,
      delta_ci_lower  = delta_ci_lo,
      delta_ci_upper  = delta_ci_hi,
      wilcox_p_value  = p_wilcox,
      n_pos_A         = n_pos_A,
      n_nonpos_A      = n_non_A,
      frac_pos_A      = frac_pos_A,
      n_pos_B         = n_pos_B,
      n_nonpos_B      = n_non_B,
      frac_pos_B      = frac_pos_B,
      pos_threshold   = pos_threshold,
      fisher_p_A_gt_B_dnds_gt_pos_threshold = fisher_p,
      stringsAsFactors = FALSE
    )

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    out_tsv <- file.path(out_dir, paste0(contrast_name, "_", side_tag, "_dnds_contrast.tsv"))
    utils::write.table(summary_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    if (make_plots) {
      .plot_contrast(merged, contrast_name, compA_name, compB_name, out_dir, side_tag, mode_label)
    }

    out_tsv
  }

  # ------------------ mode dispatch ------------------

  out_paths <- character(0)

  # Batch modes using comparison_file
  if (!is.null(comparison_file)) {
    comps      <- .read_comparisons(comparison_file)
    comp_names <- comps$comparison_name

    contrast_dir <- file.path(output_dir, "contrast")
    dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

    if (!is.null(contrast_file)) {
      # Explicit contrasts (side fields ignored if global mode)
      contr <- .read_contrasts(contrast_file)

      for (i in seq_len(nrow(contr))) {
        cn    <- contr$contrast_name[i]
        cA    <- contr$compA[i]
        cB    <- contr$compB[i]
        sideA <- contr$compA_side[i]
        sideB <- contr$compB_side[i]

        if (!cA %in% comp_names) stop("compA '", cA, "' not found in comparison_file.")
        if (!cB %in% comp_names) stop("compB '", cB, "' not found in comparison_file.")

        # Only validate sides if doing regional mode
        if (do_regions) {
          if (!sideA %in% c("query", "subject")) stop("compA_side must be 'query' or 'subject' in contrast_file.")
          if (!sideB %in% c("query", "subject")) stop("compB_side must be 'query' or 'subject' in contrast_file.")
        } else {
          # keep values but they won't be used
          if (is.na(sideA) || !nzchar(sideA)) sideA <- "query"
          if (is.na(sideB) || !nzchar(sideB)) sideB <- "query"
        }

        fileA <- file.path(output_dir, cA, paste0(cA, "_dnds_annot.tsv"))
        fileB <- file.path(output_dir, cB, paste0(cB, "_dnds_annot.tsv"))

        out_paths <- c(out_paths,
                       .run_one_contrast(cn, cA, cB, sideA, sideB, fileA, fileB, contrast_dir))
      }

      message("dN/dS contrasts complete (explicit contrast_file).")
      return(invisible(out_paths))
    }

    # Auto all-pairs
    pairs <- utils::combn(comp_names, 2, simplify = FALSE)
    for (p in pairs) {
      cA    <- p[1]
      cB    <- p[2]
      fileA <- file.path(output_dir, cA, paste0(cA, "_dnds_annot.tsv"))
      fileB <- file.path(output_dir, cB, paste0(cB, "_dnds_annot.tsv"))

      if (!do_regions) {
        # Global: run once per pair
        cn <- paste0(cA, "_vs_", cB)
        out_paths <- c(out_paths,
                       .run_one_contrast(cn, cA, cB, "query", "query", fileA, fileB, contrast_dir))
      } else {
        # Regional: run for each side (same side for A and B)
        for (sd in sides) {
          cn <- paste0(cA, "_vs_", cB)
          out_paths <- c(out_paths,
                         .run_one_contrast(cn, cA, cB, sd, sd, fileA, fileB, contrast_dir))
        }
      }
    }

    message("dN/dS contrasts complete (auto all-pairs).")
    return(invisible(out_paths))
  }

  # Single-contrast mode
  if (!is.null(dnds_annot_file_a) && !is.null(dnds_annot_file_b)) {
    compA_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file_a))
    compB_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file_b))

    if (!file.exists(dnds_annot_file_a)) stop("dnds_annot_file_a not found: ", dnds_annot_file_a)
    if (!file.exists(dnds_annot_file_b)) stop("dnds_annot_file_b not found: ", dnds_annot_file_b)

    out_dir <- dirname(dnds_annot_file_a)

    if (!do_regions) {
      # Global: run once
      cn <- paste0(compA_name, "_vs_", compB_name)
      out_paths <- c(out_paths,
                     .run_one_contrast(cn, compA_name, compB_name,
                                       "query", "query",
                                       dnds_annot_file_a, dnds_annot_file_b,
                                       out_dir))
    } else {
      # Regional: run per side (same for A and B)
      for (sd in sides) {
        cn <- paste0(compA_name, "_vs_", compB_name)
        out_paths <- c(out_paths,
                       .run_one_contrast(cn, compA_name, compB_name, sd, sd,
                                         dnds_annot_file_a, dnds_annot_file_b,
                                         out_dir))
      }
    }

    message("dN/dS contrast(s) complete for: ", compA_name, " vs ", compB_name)
    return(invisible(out_paths))
  }

  stop("Either provide (dnds_annot_file_a & dnds_annot_file_b) for single-contrast ",
       "OR comparison_file (with optional contrast_file) for batch mode.")
}
