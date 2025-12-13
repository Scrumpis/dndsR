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
#' Expects columns: contrast_name, compA, compB, side
#'
#' @keywords internal
.read_contrasts <- function(x) {
  req <- c("contrast_name", "compA", "compB", "side")
  if (is.data.frame(x)) {
    stopifnot(all(req %in% names(x)))
    return(x[, req, drop = FALSE])
  }
  df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
  if (!inherits(df1, "try-error") && all(req %in% names(df1))) {
    return(df1[, req, drop = FALSE])
  }
  df2 <- .read_ws(x, header_try = FALSE)
  if (ncol(df2) < 4) {
    stop("contrast_file must have 4 columns or a header with: ",
         paste(req, collapse = ", "))
  }
  names(df2)[1:4] <- req
  df2[, req, drop = FALSE]
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
    seq_col    <- "q_seqname"
    start_col  <- "q_start"
    end_col    <- "q_end"
    status_col <- "q_region_status"
  } else {
    seq_col    <- "s_seqname"
    start_col  <- "s_start"
    end_col    <- "s_end"
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

# -----------------------------
# 1) Regional summary function
# -----------------------------

#' Regional dN/dS summaries and region vs background tests
#'
#' For each dNdS annotation file, label rows as inside/outside a set of genomic
#' regions and summarize dN/dS distributions by region_status. Optionally tests
#' whether dN/dS differs between genes inside vs outside the regions using
#' Wilcoxon rank-sum tests.
#'
#' Works in single mode (dnds_annot_file) or batch mode (comparison_file +
#' output_dir), mirroring other dndsR commands.
#'
#' @param dnds_annot_file Path to a single <comp>_dnds_annot.tsv (single mode).
#' @param comparison_file Path to whitespace-delimited file (tabs/spaces; header or not)
#'   with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#'   If provided, batch mode reads:
#'   file.path(output_dir, comparison_name, paste0(comparison_name, "_dnds_annot.tsv")).
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#'
#' @param regions_bed Path to BED-like file of regions with at least 3 columns:
#'   seqname, start, end, and optionally a name/label column.
#' @param region_seq_col,region_start_col,region_end_col,region_name_col
#'   Column names in regions_bed for seqname/start/end/label (default: first three +
#'   optional "region_name").
#'
#' @param sides Character vector among c("query","subject") indicating which genome
#'   side to use for overlap (default both). Coordinates are taken from
#'   q_seqname/q_start/q_end or s_seqname/s_start/s_end.
#'
#' @param filter_expr Optional character with a logical expression evaluated in the data
#'   (e.g., "q_seqname == s_seqname & dNdS < 5").
#' @param max_dnds Numeric. Drop rows with dNdS >= max_dnds or NA dNdS (default 10).
#' @param make_plots Logical; if TRUE, write box/violin plots per comparison (default TRUE).
#'
#' @param ci_method One of "normal" (mean ± 1.96*SE) or "bootstrap" for 95% CI.
#' @param n_boot Number of bootstrap resamples if ci_method = "bootstrap" (default 1000).
#'
#' @return Invisibly, in single mode: path to summary TSV.
#'         In batch mode: vector of summary TSV paths.
#'         Each TSV contains, per comparison × side × region_status:
#'           comparison, side, region_status, n, mean_dnds, median_dnds,
#'           sd_dnds, se_dnds, ci_lower, ci_upper, p_value (region vs background).
#' @export
regional_dnds_summary <- function(dnds_annot_file = NULL,
                                  comparison_file = NULL,
                                  output_dir      = getwd(),
                                  regions_bed,
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

  regions <- .read_regions_bed(regions_bed,
                               region_seq_col   = region_seq_col,
                               region_start_col = region_start_col,
                               region_end_col   = region_end_col,
                               region_name_col  = region_name_col)

  .wilcox_region <- function(x_region, x_bg) {
    x_region <- x_region[is.finite(x_region)]
    x_bg     <- x_bg[is.finite(x_bg)]
    if (length(x_region) < 2L || length(x_bg) < 2L) return(NA_real_)
    stats::wilcox.test(x_region, x_bg, alternative = "two.sided")$p.value
  }

  .summarize_one <- function(comp_name, comp_dir) {
    in_file <- file.path(comp_dir, paste0(comp_name, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) {
      stop("Annotated dN/dS file not found: ", in_file)
    }
    d <- utils::read.table(in_file,
                           sep       = "\t",
                           header    = TRUE,
                           stringsAsFactors = FALSE,
                           quote     = "",
                           comment.char = "")
    d <- .filter_dnds(d, filter_expr = filter_expr, max_dnds = max_dnds)
    if (!nrow(d)) {
      warning("No rows after filter for ", comp_name)
      return(character(0))
    }

    for (sd in sides) {
      d <- .label_region_side(d, regions = regions, side = sd)
    }

    res_list <- list()

    for (sd in sides) {
      status_col <- if (sd == "query") "q_region_status" else "s_region_status"
      if (!status_col %in% names(d)) next

      for (st in c("region", "background")) {
        vals  <- d$dNdS[d[[status_col]] == st]
        stats <- .calc_mean_ci(vals, ci_method = ci_method, n_boot = n_boot)
        res_list[[paste(sd, st, sep = "_")]] <- data.frame(
          comparison    = comp_name,
          side          = sd,
          region_status = st,
          n             = sum(is.finite(vals)),
          mean_dnds     = stats[1],
          median_dnds   = stats[2],
          sd_dnds       = stats[3],
          se_dnds       = stats[4],
          ci_lower      = stats[5],
          ci_upper      = stats[6],
          stringsAsFactors = FALSE
        )
      }

      x_reg <- d$dNdS[d[[status_col]] == "region"]
      x_bg  <- d$dNdS[d[[status_col]] == "background"]
      p     <- .wilcox_region(x_reg, x_bg)
      for (nm in names(res_list)) {
        if (startsWith(nm, paste0(sd, "_"))) {
          res_list[[nm]]$p_value <- p
        }
      }
    }

    if (!length(res_list)) {
      warning("No summaries computed for ", comp_name)
      return(character(0))
    }
    summary_df <- do.call(rbind, res_list)

    out_tsv <- file.path(comp_dir, paste0(comp_name, "_regional_dnds_summary.tsv"))
    utils::write.table(summary_df, file = out_tsv, sep = "\t",
                       quote = FALSE, row.names = FALSE)

    if (make_plots && requireNamespace("ggplot2", quietly = TRUE)) {
      for (sd in sides) {
        status_col <- if (sd == "query") "q_region_status" else "s_region_status"
        if (!status_col %in% names(d)) next

        df_plot <- d[is.finite(d$dNdS) & !is.na(d[[status_col]]), ]
        if (!nrow(df_plot)) next

        gg <- ggplot2::ggplot(df_plot,
                              ggplot2::aes_string(x = status_col, y = "dNdS")) +
          ggplot2::geom_violin(trim = TRUE, alpha = 0.7) +
          ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.5) +
          ggplot2::scale_y_continuous(trans = "log1p") +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::labs(x     = paste0(sd, " region status"),
                        y     = "dN/dS",
                        title = paste0(comp_name, " (", sd, " side)"))

        out_pdf <- file.path(comp_dir,
                             paste0(comp_name, "_regional_dnds_", sd, "_violin.pdf"))
        ggplot2::ggsave(out_pdf, gg, width = 5, height = 4)
      }
    }

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
      outs <- c(outs, .summarize_one(comp, comp_dir))
    }
    message("Regional dN/dS summaries complete.")
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
  message("Regional dN/dS summary complete for: ", dnds_annot_file)
  invisible(out)
}

# -----------------------------
# 2) Regional contrasts function
# -----------------------------

#' Pairwise regional dN/dS contrasts (Wilcoxon signed-rank)
#'
#' Compare dN/dS between pairs of dNdS annotation files within specified genomic
#' regions using paired nonparametric tests. For each contrast (compA vs compB),
#' the function:
#'   - restricts to genes in the provided regions (query or subject side),
#'   - matches rows across comparisons by user-specified ID columns,
#'   - computes delta = dNdS_A - dNdS_B per matched gene,
#'   - runs a Wilcoxon signed-rank test on delta,
#'   - summarizes statistics and writes a TSV (and optional plots) per contrast.
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
#'   Used in batch mode together with contrast_file or auto-all-pairs.
#' @param contrast_file Optional path to a whitespace-delimited file defining pairwise
#'   contrasts. Required columns:
#'     contrast_name, compA, compB, side
#'   where side is "query" or "subject" indicating which genome's coordinates
#'   should be overlapped with regions. If contrast_file is NULL, all unique
#'   unordered pairs of comparisons are tested for each side in 'sides'.
#'
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#'
#' @param regions_bed Path to BED-like file of regions with at least 3 columns:
#'   seqname, start, end, and optionally a name/label column.
#' @param region_seq_col,region_start_col,region_end_col,region_name_col
#'   Column names in regions_bed for seqname/start/end/label.
#'
#' @param merge_cols Character vector of column names used to match rows across
#'   comparisons (e.g. "query_id", "subject_id", or an orthogroup ID). If NULL
#'   (default), the function uses "query_id" when side = "query" and "subject_id"
#'   when side = "subject".
#'
#' @param sides Character vector among c("query","subject") indicating which side(s)
#'   to evaluate in auto-all-pairs mode and single-contrast mode. Ignored when
#'   contrast_file explicitly specifies side.
#'
#' @param filter_expr Optional character with a logical expression evaluated in
#'   each dNdS table before merging (e.g., "q_seqname == s_seqname & dNdS < 5").
#' @param max_dnds Numeric. Drop rows with dNdS >= max_dnds or NA dNdS (default 10).
#'
#' @param ci_method One of "normal" (mean ± 1.96*SE) or "bootstrap" (percentile).
#' @param n_boot Integer; number of bootstrap resamples if ci_method = "bootstrap".
#'
#' @param make_plots Logical; if TRUE, write a paired scatter and delta histogram
#'   per contrast and side.
#'
#' @return Invisibly, a character vector of summary TSV paths (one per contrast × side).
#'         Each TSV has columns:
#'           contrast_name, compA, compB, side, n,
#'           mean_dnds_A, median_dnds_A,
#'           mean_dnds_B, median_dnds_B,
#'           mean_delta, median_delta,
#'           delta_sd, delta_se, delta_ci_lower, delta_ci_upper,
#'           wilcox_p_value.
#' @export
regional_dnds_contrasts <- function(dnds_annot_file_a = NULL,
                                    dnds_annot_file_b = NULL,
                                    comparison_file   = NULL,
                                    contrast_file     = NULL,
                                    output_dir        = getwd(),
                                    regions_bed,
                                    region_seq_col   = NULL,
                                    region_start_col = NULL,
                                    region_end_col   = NULL,
                                    region_name_col  = NULL,
                                    merge_cols       = NULL,
                                    sides            = c("query", "subject"),
                                    filter_expr      = NULL,
                                    max_dnds         = 10,
                                    ci_method        = c("normal", "bootstrap"),
                                    n_boot           = 1000,
                                    make_plots       = TRUE) {
  ci_method <- match.arg(ci_method)
  sides     <- match.arg(sides, choices = c("query", "subject"), several.ok = TRUE)

  regions <- .read_regions_bed(regions_bed,
                               region_seq_col   = region_seq_col,
                               region_start_col = region_start_col,
                               region_end_col   = region_end_col,
                               region_name_col  = region_name_col)

  .wilcox_signed <- function(delta) {
    delta <- delta[is.finite(delta) & delta != 0]
    if (length(delta) < 2L) return(NA_real_)
    stats::wilcox.test(delta, mu = 0, alternative = "two.sided")$p.value
  }

  .run_one_contrast <- function(contrast_name,
                                compA_name,
                                compB_name,
                                side,
                                fileA,
                                fileB,
                                out_dir) {

    if (!file.exists(fileA)) stop("dNdS file A not found: ", fileA)
    if (!file.exists(fileB)) stop("dNdS file B not found: ", fileB)

    dA <- utils::read.table(fileA,
                            sep       = "\t",
                            header    = TRUE,
                            stringsAsFactors = FALSE,
                            quote     = "",
                            comment.char = "")
    dB <- utils::read.table(fileB,
                            sep       = "\t",
                            header    = TRUE,
                            stringsAsFactors = FALSE,
                            quote     = "",
                            comment.char = "")

    dA <- .filter_dnds(dA, filter_expr = filter_expr, max_dnds = max_dnds)
    dB <- .filter_dnds(dB, filter_expr = filter_expr, max_dnds = max_dnds)

    if (!nrow(dA) || !nrow(dB)) {
      warning("No rows after filtering for contrast ", contrast_name)
      return(NA_character_)
    }

    dA <- .label_region_side(dA, regions = regions, side = side)
    dB <- .label_region_side(dB, regions = regions, side = side)

    status_col <- if (side == "query") "q_region_status" else "s_region_status"

    dA_reg <- dA[dA[[status_col]] == "region", , drop = FALSE]
    dB_reg <- dB[dB[[status_col]] == "region", , drop = FALSE]

    if (!nrow(dA_reg) || !nrow(dB_reg)) {
      warning("No region rows for contrast ", contrast_name)
      return(NA_character_)
    }

    # Decide which columns to use for matching rows across comparisons
    merge_cols_eff <- merge_cols
    if (is.null(merge_cols_eff) || identical(merge_cols_eff, "auto")) {
      merge_cols_eff <- if (side == "query") "query_id" else "subject_id"
    }

    missing_colsA <- setdiff(merge_cols_eff, names(dA_reg))
    missing_colsB <- setdiff(merge_cols_eff, names(dB_reg))
    if (length(missing_colsA) || length(missing_colsB)) {
      stop("merge_cols must be present in both dNdS tables. Missing: A=",
           paste(missing_colsA, collapse = ","), " B=",
           paste(missing_colsB, collapse = ","))
    }

    merged <- merge(
      dA_reg[, c(merge_cols_eff, "dNdS"), drop = FALSE],
      dB_reg[, c(merge_cols_eff, "dNdS"), drop = FALSE],
      by       = merge_cols_eff,
      suffixes = c("_A", "_B")
    )

    if (!nrow(merged)) {
      warning("No overlapping IDs between comps A and B in regions for contrast ",
              contrast_name)
      return(NA_character_)
    }

    merged$delta   <- merged$dNdS_A - merged$dNdS_B
    stats_delta    <- .calc_mean_ci(merged$delta, ci_method = ci_method, n_boot = n_boot)
    mean_delta     <- stats_delta[1]
    median_delta   <- stats_delta[2]
    delta_sd       <- stats_delta[3]
    delta_se       <- stats_delta[4]
    delta_ci_lo    <- stats_delta[5]
    delta_ci_hi    <- stats_delta[6]

    mean_A <- mean(merged$dNdS_A, na.rm = TRUE)
    med_A  <- stats::median(merged$dNdS_A, na.rm = TRUE)
    mean_B <- mean(merged$dNdS_B, na.rm = TRUE)
    med_B  <- stats::median(merged$dNdS_B, na.rm = TRUE)
    n      <- sum(is.finite(merged$delta))

    p_wilcox <- .wilcox_signed(merged$delta)

    summary_df <- data.frame(
      contrast_name   = contrast_name,
      compA           = compA_name,
      compB           = compB_name,
      side            = side,
      n               = n,
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
      stringsAsFactors = FALSE
    )

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    out_tsv <- file.path(out_dir,
                         paste0(contrast_name, "_", side, "_regional_dnds_contrast.tsv"))
    utils::write.table(summary_df, file = out_tsv, sep = "\t",
                       quote = FALSE, row.names = FALSE)

    if (make_plots && requireNamespace("ggplot2", quietly = TRUE)) {
      # Scatter plot dNdS_B vs dNdS_A
      gg <- ggplot2::ggplot(merged,
                            ggplot2::aes(x = dNdS_B, y = dNdS_A)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        ggplot2::geom_point(alpha = 0.4, size = 1) +
        ggplot2::scale_x_continuous(trans = "log1p") +
        ggplot2::scale_y_continuous(trans = "log1p") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(
          title = paste0(contrast_name, " (", side, " side, region-only)"),
          x     = paste0("dN/dS (", compB_name, ")"),
          y     = paste0("dN/dS (", compA_name, ")")
        )

      out_pdf <- file.path(out_dir,
                           paste0(contrast_name, "_", side, "_regional_dnds_scatter.pdf"))
      ggplot2::ggsave(out_pdf, gg, width = 5, height = 4)

      # Delta histogram
      gg2 <- ggplot2::ggplot(merged,
                             ggplot2::aes(x = delta)) +
        ggplot2::geom_histogram(bins = 40, alpha = 0.7) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(
          x     = "delta dN/dS (A - B)",
          y     = "Count",
          title = paste0(contrast_name, " — delta distribution (region-only, ", side, ")")
        )
      out_pdf2 <- file.path(out_dir,
                            paste0(contrast_name, "_", side, "_regional_dnds_delta_hist.pdf"))
      ggplot2::ggsave(out_pdf2, gg2, width = 5, height = 4)
    }

    out_tsv
  }

  # ------------------ mode dispatch ------------------

  out_paths <- character(0)

  # Batch modes using comparison_file
  if (!is.null(comparison_file)) {
    comps      <- .read_comparisons(comparison_file)
    comp_names <- comps$comparison_name

    if (!is.null(contrast_file)) {
      # Explicit contrasts
      contr <- .read_contrasts(contrast_file)
      for (i in seq_len(nrow(contr))) {
        cn   <- contr$contrast_name[i]
        cA   <- contr$compA[i]
        cB   <- contr$compB[i]
        side <- contr$side[i]

        if (!cA %in% comp_names) stop("compA '", cA, "' not found in comparison_file.")
        if (!cB %in% comp_names) stop("compB '", cB, "' not found in comparison_file.")
        if (!side %in% c("query", "subject")) {
          stop("side must be 'query' or 'subject' in contrast_file.")
        }

        fileA <- file.path(output_dir, cA, paste0(cA, "_dnds_annot.tsv"))
        fileB <- file.path(output_dir, cB, paste0(cB, "_dnds_annot.tsv"))

        contrast_dir <- file.path(output_dir, "regional_contrasts")
        out_paths <- c(out_paths,
                       .run_one_contrast(cn, cA, cB, side, fileA, fileB, contrast_dir))
      }
      message("Regional dN/dS contrasts complete (explicit contrast_file).")
      return(invisible(out_paths))
    }

    # Auto all-pairs mode
    pairs        <- utils::combn(comp_names, 2, simplify = FALSE)
    contrast_dir <- file.path(output_dir, "regional_contrasts")
    for (p in pairs) {
      cA    <- p[1]
      cB    <- p[2]
      fileA <- file.path(output_dir, cA, paste0(cA, "_dnds_annot.tsv"))
      fileB <- file.path(output_dir, cB, paste0(cB, "_dnds_annot.tsv"))
      for (sd in sides) {
        cn <- paste0(cA, "_vs_", cB)
        out_paths <- c(out_paths,
                       .run_one_contrast(cn, cA, cB, sd, fileA, fileB, contrast_dir))
      }
    }
    message("Regional dN/dS contrasts complete (auto all-pairs).")
    return(invisible(out_paths))
  }

  # Single-contrast mode
  if (!is.null(dnds_annot_file_a) && !is.null(dnds_annot_file_b)) {
    compA_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file_a))
    compB_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file_b))

    if (!file.exists(dnds_annot_file_a)) {
      stop("dnds_annot_file_a not found: ", dnds_annot_file_a)
    }
    if (!file.exists(dnds_annot_file_b)) {
      stop("dnds_annot_file_b not found: ", dnds_annot_file_b)
    }

    out_dir <- dirname(dnds_annot_file_a)
    for (sd in sides) {
      cn <- paste0(compA_name, "_vs_", compB_name)
      out_paths <- c(out_paths,
                     .run_one_contrast(cn, compA_name, compB_name, sd,
                                       dnds_annot_file_a, dnds_annot_file_b,
                                       out_dir))
    }
    message("Regional dN/dS contrast(s) complete for: ",
            compA_name, " vs ", compB_name)
    return(invisible(out_paths))
  }

  stop("Either provide (dnds_annot_file_a & dnds_annot_file_b) for single-contrast ",
       "OR comparison_file (with optional contrast_file) for batch mode.")
}
