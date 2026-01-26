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

#' Read a contrast file defining pairwise regional dN/dS comparisons
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
#' @keywords internal
.read_contrasts <- function(x) {
  req5 <- c("contrast_name", "compA", "compA_side", "compB", "compB_side")
  req4 <- c("contrast_name", "compA", "compB", "side")

  normalize_df <- function(df) {
    # If already 5-column style
    if (all(req5 %in% names(df))) {
      return(df[, req5, drop = FALSE])
    }
    # If 4-column style, expand side -> compA_side / compB_side
    if (all(req4 %in% names(df))) {
      df_out <- data.frame(
        contrast_name = df$contrast_name,
        compA         = df$compA,
        compA_side    = df$side,
        compB         = df$compB,
        compB_side    = df$side,
        stringsAsFactors = FALSE
      )
      return(df_out)
    }
    stop("contrast_file must have either 4 cols (",
         paste(req4, collapse = ", "),
         ") or 5 cols (",
         paste(req5, collapse = ", "),
         "), or a header with those names.")
  }

  if (is.data.frame(x)) {
    return(normalize_df(x))
  }

  # Try with header (only accept if header matches expected names)
  df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
  if (!inherits(df1, "try-error")) {
    if (all(req5 %in% names(df1)) || all(req4 %in% names(df1))) {
      return(normalize_df(df1))
    }
    # else fall through to no-header parsing
  }

  # Fallback: no header, assign by position
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
                   envir = d,
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
  if (missing(regions_bed)) {
    stop("regions_bed is required.")
  }
  if (!file.exists(regions_bed)) {
    stop("regions_bed not found: ", regions_bed)
  }
  reg_raw <- .read_ws(regions_bed, header_try = TRUE)
  if (ncol(reg_raw) < 3) {
    stop("regions_bed must have at least 3 columns: seqname, start, end.")
  }
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
    ov    <- data.table::foverlaps(dt_g, dt_r, nomatch = 0L)
    in_idx <- unique(ov$idx)
    lab   <- rep("background", nrow(d))
    lab[in_idx] <- "region"
    d[[status_col]] <- lab
  } else {
    lab <- rep("background", nrow(d))
    for (i in seq_len(nrow(regions))) {
      r    <- regions[i, ]
      hits <- which(d[[seq_col]] == r$seqname &
                    d[[start_col]] <= r$end   &
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
  if (!n) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_))
  }
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
  nn <- sub("\r$", "", nn)   # strip CRLF artifacts
  nn <- trimws(nn)          # strip leading/trailing spaces
  names(d) <- nn
  d
}

# -----------------------------
# 1) Regional summary function
# -----------------------------

#' dN/dS summaries (global + optional region vs background tests)
#'
#' Summarize dN/dS distributions for each dNdS annotation file.
#'
#' - Always computes a **global** summary across all retained rows.
#' - If `regions_bed` is provided, also labels rows as inside/outside regions
#'   (per side) and computes region/background summaries + a Wilcoxon rank-sum
#'   test (region vs background) for each side.
#'
#' Works in single mode (`dnds_annot_file`) or batch mode (`comparison_file` +
#' `output_dir`), mirroring other dndsR commands.
#'
#' @param dnds_annot_file Path to a single <comp>_dnds_annot.tsv (single mode).
#' @param comparison_file Path to whitespace-delimited file (tabs/spaces; header or not)
#'   with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#'   If provided, batch mode reads:
#'   file.path(output_dir, comparison_name, paste0(comparison_name, "_dnds_annot.tsv")).
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#'
#' @param regions_bed Optional. Path to BED-like file of regions with at least 3 columns:
#'   seqname, start, end, and optionally a name/label column. If NULL, only global
#'   summaries (and global plots) are produced.
#' @param region_seq_col,region_start_col,region_end_col,region_name_col
#'   Column names in regions_bed for seqname/start/end/label (default: first three +
#'   optional "region_name").
#'
#' @param sides Character vector among c("query","subject") indicating which genome
#'   side to use for overlap when `regions_bed` is provided (default both).
#'   Coordinates are taken from q_gff_* or s_gff_* columns.
#'
#' @param filter_expr Optional character with a logical expression evaluated in the data
#'   (e.g., "q_gff_seqname == s_gff_seqname & dNdS < 5").
#' @param max_dnds Numeric. Drop rows with dNdS >= max_dnds or NA dNdS (default 10).
#' @param make_plots Logical; if TRUE, write global plots always and regional plots
#'   when `regions_bed` is provided (default TRUE).
#'
#' @param ci_method One of "normal" (mean +/- 1.96*SE) or "bootstrap" for 95% CI.
#' @param n_boot Number of bootstrap resamples if ci_method = "bootstrap" (default 1000).
#'
#' @return Invisibly, in single mode: path to summary TSV.
#'         In batch mode: vector of summary TSV paths.
#'         The TSV contains:
#'           - one row with scope="global" (side=NA, region_status="all")
#'           - plus (if regions_bed) two rows per side with scope="regional"
#'             (region/background) and a Wilcoxon rank-sum p_value for region vs background.
#' @export
dnds_summary <- function(dnds_annot_file = NULL,
                         comparison_file = NULL,
                         output_dir      = getwd(),
                         regions_bed     = NULL,
                         region_seq_col   = NULL,
                         region_start_col = NULL,
                         region_end_col   = NULL,
                         region_name_col  = NULL,
                         sides            = c("query", "subject"),
                         filter_expr      = NULL,
                         max_dnds         = 10,
                         make_plots       = TRUE,
                         ci_method        = c("normal", "bootstrap"),
                         n_boot           = 1000) {

  ci_method <- match.arg(ci_method)
  sides     <- match.arg(sides, choices = c("query", "subject"), several.ok = TRUE)

  has_regions <- !is.null(regions_bed)
  regions <- NULL
  if (has_regions) {
    regions <- .read_regions_bed(regions_bed,
                                 region_seq_col   = region_seq_col,
                                 region_start_col = region_start_col,
                                 region_end_col   = region_end_col,
                                 region_name_col  = region_name_col)
  }

  .wilcox_region <- function(x_region, x_bg) {
    x_region <- x_region[is.finite(x_region)]
    x_bg     <- x_bg[is.finite(x_bg)]
    if (length(x_region) < 2L || length(x_bg) < 2L) return(NA_real_)
    stats::wilcox.test(x_region, x_bg, alternative = "two.sided")$p.value
  }

  .plot_global <- function(d, comp_name, comp_dir) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
    if (!nrow(d) || !"dNdS" %in% names(d)) return(invisible(NULL))

    dfp <- d[is.finite(d$dNdS), , drop = FALSE]
    if (!nrow(dfp)) return(invisible(NULL))

    # Global histogram (log1p x)
    gg_hist <- ggplot2::ggplot(dfp, ggplot2::aes(x = dNdS)) +
      ggplot2::geom_histogram(bins = 60, alpha = 0.7) +
      ggplot2::scale_x_continuous(trans = "log1p") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = paste0(comp_name, " — global dN/dS distribution"),
        x     = "dN/dS (log1p scale)",
        y     = "Count"
      )

    ggplot2::ggsave(
      filename = file.path(comp_dir, paste0(comp_name, "_dnds_global_hist.pdf")),
      plot     = gg_hist,
      width    = 6,
      height   = 4
    )

    # Global ECDF (log1p x)
    gg_ecdf <- ggplot2::ggplot(dfp, ggplot2::aes(x = dNdS)) +
      ggplot2::stat_ecdf() +
      ggplot2::scale_x_continuous(trans = "log1p") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = paste0(comp_name, " — global dN/dS ECDF"),
        x     = "dN/dS (log1p scale)",
        y     = "ECDF"
      )

    ggplot2::ggsave(
      filename = file.path(comp_dir, paste0(comp_name, "_dnds_global_ecdf.pdf")),
      plot     = gg_ecdf,
      width    = 6,
      height   = 4
    )

    invisible(NULL)
  }

  .plot_regional <- function(d, comp_name, comp_dir, sides) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
    if (!nrow(d) || !"dNdS" %in% names(d)) return(invisible(NULL))

    for (sd in sides) {
      status_col <- if (sd == "query") "q_region_status" else "s_region_status"
      if (!status_col %in% names(d)) next

      df_plot <- d[is.finite(d$dNdS) & !is.na(d[[status_col]]), ]
      if (!nrow(df_plot)) next

      # Violin + box (quick region vs background comparison)
      gg <- ggplot2::ggplot(df_plot,
                            ggplot2::aes_string(x = status_col, y = "dNdS")) +
        ggplot2::geom_violin(trim = TRUE, alpha = 0.7) +
        ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.5) +
        ggplot2::scale_y_continuous(trans = "log1p") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(
          x     = paste0(sd, " region status"),
          y     = "dN/dS (log1p scale)",
          title = paste0(comp_name, " — region vs background (", sd, " side)")
        )

      ggplot2::ggsave(
        filename = file.path(comp_dir, paste0(comp_name, "_regional_dnds_", sd, "_violin.pdf")),
        plot     = gg,
        width    = 5,
        height   = 4
      )

      # Overlaid histograms (often easier to interpret than violins)
      gg2 <- ggplot2::ggplot(df_plot, ggplot2::aes_string(x = "dNdS", fill = status_col)) +
        ggplot2::geom_histogram(position = "identity", bins = 60, alpha = 0.45) +
        ggplot2::scale_x_continuous(trans = "log1p") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(
          title = paste0(comp_name, " — region/background overlay (", sd, " side)"),
          x     = "dN/dS (log1p scale)",
          y     = "Count"
        )

      ggplot2::ggsave(
        filename = file.path(comp_dir, paste0(comp_name, "_regional_dnds_", sd, "_hist_overlay.pdf")),
        plot     = gg2,
        width    = 6,
        height   = 4
      )
    }

    invisible(NULL)
  }

  .summarize_one <- function(comp_name, comp_dir) {
    in_file <- file.path(comp_dir, paste0(comp_name, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) stop("Annotated dN/dS file not found: ", in_file)

    d_raw <- utils::read.table(in_file,
                               sep = "\t",
                               header = TRUE,
                               stringsAsFactors = FALSE,
                               quote = "",
                               comment.char = "",
                               check.names = FALSE)
    d_raw <- .clean_colnames(d_raw)

    if (!nrow(d_raw)) {
      warning("No rows in file for ", comp_name)
      return(character(0))
    }
    if (!"dNdS" %in% names(d_raw)) {
      stop("Missing required column 'dNdS' in: ", in_file)
    }

    # Bookkeeping before filtering
    n_total <- nrow(d_raw)
    n_removed_na     <- sum(is.na(d_raw$dNdS))
    n_removed_ge_max <- sum(is.finite(d_raw$dNdS) & d_raw$dNdS >= max_dnds)

    # Apply filtering (was previously missing in the original function)
    d <- .filter_dnds(d_raw, filter_expr = filter_expr, max_dnds = max_dnds)

    if (!nrow(d)) {
      warning("No rows after filtering for ", comp_name)
      return(character(0))
    }

    # Always plot global patterns (if requested)
    if (make_plots) .plot_global(d, comp_name, comp_dir)

    res_list <- list()

    # ---------- global summary ----------
    vals_all  <- d$dNdS
    stats_all <- .calc_mean_ci(vals_all, ci_method = ci_method, n_boot = n_boot)

    res_list[["global"]] <- data.frame(
      comparison        = comp_name,
      scope             = "global",
      side              = NA_character_,
      region_status     = "all",
      n_total           = n_total,
      n_kept            = nrow(d),
      n_removed_na      = n_removed_na,
      n_removed_ge_max  = n_removed_ge_max,
      max_dnds          = max_dnds,
      n                 = sum(is.finite(vals_all)),
      mean_dnds         = stats_all[1],
      median_dnds       = stats_all[2],
      sd_dnds           = stats_all[3],
      se_dnds           = stats_all[4],
      ci_lower          = stats_all[5],
      ci_upper          = stats_all[6],
      p_value           = NA_real_,
      stringsAsFactors  = FALSE
    )

    # ---------- regional summaries (optional) ----------
    if (has_regions) {
      # Label region status for each requested side
      for (sd in sides) {
        d <- tryCatch(
          .label_region_side(d, regions = regions, side = sd),
          error = function(e) {
            stop(
              sprintf(
                "dnds_summary failed for comparison '%s' (file: %s) while labeling side='%s':\n%s",
                comp_name, in_file, sd, conditionMessage(e)
              ),
              call. = FALSE
            )
          }
        )
      }

      # Regional plots (if requested)
      if (make_plots) .plot_regional(d, comp_name, comp_dir, sides)

      for (sd in sides) {
        status_col <- if (sd == "query") "q_region_status" else "s_region_status"
        if (!status_col %in% names(d)) next

        # Stats for region/background
        for (st in c("region", "background")) {
          vals  <- d$dNdS[d[[status_col]] == st]
          stats <- .calc_mean_ci(vals, ci_method = ci_method, n_boot = n_boot)
          res_list[[paste0("regional__", sd, "__", st)]] <- data.frame(
            comparison        = comp_name,
            scope             = "regional",
            side              = sd,
            region_status     = st,
            n_total           = n_total,
            n_kept            = nrow(d),
            n_removed_na      = n_removed_na,
            n_removed_ge_max  = n_removed_ge_max,
            max_dnds          = max_dnds,
            n                 = sum(is.finite(vals)),
            mean_dnds         = stats[1],
            median_dnds       = stats[2],
            sd_dnds           = stats[3],
            se_dnds           = stats[4],
            ci_lower          = stats[5],
            ci_upper          = stats[6],
            p_value           = NA_real_,  # filled below for both rows
            stringsAsFactors  = FALSE
          )
        }

        # Wilcoxon rank-sum: region vs background (within comparison/side)
        x_reg <- d$dNdS[d[[status_col]] == "region"]
        x_bg  <- d$dNdS[d[[status_col]] == "background"]
        p     <- .wilcox_region(x_reg, x_bg)

        for (nm in names(res_list)) {
          if (startsWith(nm, paste0("regional__", sd, "__"))) {
            res_list[[nm]]$p_value <- p
          }
        }
      }
    }

    summary_df <- do.call(rbind, res_list)

    out_tsv <- file.path(comp_dir, paste0(comp_name, "_dnds_summary.tsv"))
    utils::write.table(summary_df, file = out_tsv, sep = "\t",
                       quote = FALSE, row.names = FALSE)

    out_tsv
  }

  # ---- batch vs single ----
  if (!is.null(comparison_file)) {
    df   <- .read_comparisons(comparison_file)
    outs <- character(0)

    for (i in seq_len(nrow(df))) {
      comp     <- df$comparison_name[i]
      comp_dir <- file.path(output_dir, comp)
      dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

      one <- tryCatch(
        .summarize_one(comp, comp_dir),
        error = function(e) {
          warning(
            sprintf(
              "Skipping comparison '%s' (dir: %s): %s",
              comp, comp_dir, conditionMessage(e)
            ),
            call. = FALSE
          )
          character(0)
        }
      )

      outs <- c(outs, one)
    }

    message("dN/dS summaries complete.")
    return(invisible(outs))
  }

  if (is.null(dnds_annot_file)) {
    stop("Provide either comparison_file (batch) OR dnds_annot_file (single).")
  }
  if (!file.exists(dnds_annot_file)) {
    stop("dnds_annot_file not found: ", dnds_annot_file)
  }

  comp_dir  <- dirname(dnds_annot_file)
  comp_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file))
  out       <- .summarize_one(comp_name, comp_dir)
  message("dN/dS summary complete for: ", dnds_annot_file)
  invisible(out)
}
