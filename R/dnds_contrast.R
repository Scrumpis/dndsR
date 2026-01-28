# dnds_contrast.R
#
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
        contrast_name    = df$contrast_name,
        compA            = df$compA,
        compA_side       = df$side,
        compB            = df$compB,
        compB_side       = df$side,
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

#' Filter dNdS annotation table by NA / max_dnds
#'
#' Keeps rows with finite dNdS and dNdS < max_dnds. This is the only filtering
#' performed (by design).
#'
#' @keywords internal
.filter_dnds <- function(d, max_dnds = 10) {
  keep <- is.finite(d$dNdS) & d$dNdS < max_dnds
  d[keep, , drop = FALSE]
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

#' Read and normalize a regions BED-like file
#'
#' Returns a data.frame with columns: seqname, start, end, region_name
#'
#' Coordinate conventions:
#'   - If regions_coord = "bed0", interpret as BED 0-based, half-open [start, end)
#'     and convert to 1-based, closed [start+1, end] for overlap with GFF-like coordinates.
#'   - If regions_coord = "gff1", interpret as 1-based, closed [start, end].
#'
#' Robust to headerless BED files: tries header=TRUE, then falls back to header=FALSE.
#'
#' @keywords internal
.read_regions_bed <- function(regions_bed,
                              region_seq_col   = NULL,
                              region_start_col = NULL,
                              region_end_col   = NULL,
                              region_name_col  = NULL,
                              regions_coord    = c("bed0", "gff1"),
                              stop_on_invalid  = TRUE,
                              ...) {
  regions_coord <- match.arg(regions_coord)

  if (missing(regions_bed) || is.null(regions_bed) || !nzchar(regions_bed)) {
    stop("regions_bed is required.")
  }
  if (!file.exists(regions_bed)) stop("regions_bed not found: ", regions_bed)

  reg_raw <- try(.read_ws(regions_bed, header_try = TRUE), silent = TRUE)
  if (inherits(reg_raw, "try-error") || ncol(reg_raw) < 3) {
    reg_raw <- .read_ws(regions_bed, header_try = FALSE)
  }
  if (ncol(reg_raw) < 3) stop("regions_bed must have at least 3 columns: seqname, start, end.")

  # If headerless, assign standard names for first 3 columns.
  if (is.null(names(reg_raw)) || any(!nzchar(names(reg_raw)[1:3]))) {
    names(reg_raw)[1:3] <- c("seqname", "start", "end")
  }

  if (is.null(region_seq_col))   region_seq_col   <- names(reg_raw)[1]
  if (is.null(region_start_col)) region_start_col <- names(reg_raw)[2]
  if (is.null(region_end_col))   region_end_col   <- names(reg_raw)[3]
  if (is.null(region_name_col) || !region_name_col %in% names(reg_raw)) {
    reg_raw$region_name <- "region"
    region_name_col <- "region_name"
  }

  seqname <- as.character(reg_raw[[region_seq_col]])
  start0  <- suppressWarnings(as.numeric(reg_raw[[region_start_col]]))
  end0    <- suppressWarnings(as.numeric(reg_raw[[region_end_col]]))
  rname   <- as.character(reg_raw[[region_name_col]])

  bad <- is.na(seqname) | !nzchar(seqname) | is.na(start0) | is.na(end0)
  if (any(bad)) {
    msg <- paste0(
      "regions_bed has invalid rows (missing/NA seqname or non-numeric start/end). ",
      "First few bad row indices: ",
      paste(utils::head(which(bad), 10), collapse = ", ")
    )
    if (stop_on_invalid) stop(msg) else warning(msg)
  }

  # Convert coordinate convention to 1-based closed intervals for downstream overlap
  if (regions_coord == "bed0") {
    # BED: 0-based, half-open [start, end) => 1-based closed [start+1, end]
    start <- start0 + 1
    end   <- end0
  } else {
    start <- start0
    end   <- end0
  }

  # Basic sanity
  bad2 <- !is.finite(start) | !is.finite(end) | start > end
  if (any(bad2)) {
    msg <- paste0(
      "regions_bed has invalid intervals (start > end or non-finite). ",
      "First few bad row indices: ",
      paste(utils::head(which(bad2), 10), collapse = ", ")
    )
    if (stop_on_invalid) stop(msg) else warning(msg)
  }

  data.frame(
    seqname     = seqname,
    start       = start,
    end         = end,
    region_name = rname,
    stringsAsFactors = FALSE
  )
}

#' Label dNdS rows as inside/outside regions for a given side
#'
#' Adds a column q_region_status or s_region_status with values "region"/"background".
#' Uses data.table::foverlaps if available; otherwise falls back to a simple loop.
#'
#' NOTE: Assumes regions and gene coords are in the same convention (1-based closed).
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

  # coerce coords numeric & validate
  d_seq   <- as.character(d[[seq_col]])
  d_start <- suppressWarnings(as.numeric(d[[start_col]]))
  d_end   <- suppressWarnings(as.numeric(d[[end_col]]))

  bad <- is.na(d_seq) | !nzchar(d_seq) | is.na(d_start) | is.na(d_end) | d_start > d_end
  if (any(bad)) {
    stop("Gene coordinate columns contain invalid rows for side='", side, "'. ",
         "First few bad row indices: ", paste(utils::head(which(bad), 10), collapse = ", "))
  }

  if (has_dt) {
    dt_g <- data.table::data.table(
      idx     = seq_len(nrow(d)),
      seqname = d_seq,
      start   = d_start,
      end     = d_end
    )
    dt_r <- data.table::data.table(
      seqname     = as.character(regions$seqname),
      start       = as.numeric(regions$start),
      end         = as.numeric(regions$end),
      region_name = as.character(regions$region_name)
    )

    # Validate regions too
    badr <- is.na(dt_r$seqname) | !nzchar(dt_r$seqname) |
      is.na(dt_r$start) | is.na(dt_r$end) | dt_r$start > dt_r$end
    if (any(badr)) {
      stop("Regions contain invalid rows. First few bad region indices: ",
           paste(utils::head(which(badr), 10), collapse = ", "))
    }

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
      hits <- which(d_seq == r$seqname &
                      d_start <= r$end &
                      d_end   >= r$start)
      if (length(hits)) lab[hits] <- "region"
    }
    d[[status_col]] <- lab
  }
  d
}

#' Compute mean/median/SD/SE/CI for a numeric vector
#'
#' @keywords internal
.calc_mean_ci <- function(x,
                          ci_method = c("normal", "bootstrap"),
                          n_boot = 1000,
                          seed = NULL,
                          ...) {
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
    if (!is.null(seed)) set.seed(seed)
    means <- replicate(n_boot, mean(sample(x, size = n, replace = TRUE)))
    ci <- stats::quantile(means, probs = c(0.025, 0.975), names = FALSE)
  }
  c(m, med, sd, se, ci[1], ci[2])
}

# -----------------------------
# dN/dS contrasts function
# -----------------------------

#' Pairwise dN/dS contrasts (genome-wide by default; optional region restriction)
#'
#' Compare dN/dS between pairs of dNdS annotation files using paired nonparametric
#' tests. For each contrast (compA vs compB), the function:
#'   - filters each table by max_dnds,
#'   - (optional) restricts to genes in regions_bed (side-aware),
#'   - matches rows across comparisons,
#'   - computes delta = dNdS_A - dNdS_B per matched gene,
#'   - runs a Wilcoxon signed-rank test on delta,
#'   - runs McNemar's test (paired) for enrichment of dN/dS > pos_threshold in A vs B,
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
#' @param regions_coord Coordinate convention for regions_bed: "bed0" (BED 0-based half-open)
#'   or "gff1" (1-based closed). Default "bed0".
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
#' @param max_dnds Numeric. Drop rows with dNdS >= max_dnds or non-finite dNdS (default 10).
#'
#' @param ci_method One of "normal" (mean +/- 1.96*SE) or "bootstrap" (percentile).
#' @param n_boot Integer; number of bootstrap resamples if ci_method = "bootstrap".
#' @param ci_seed Optional integer seed for bootstrap reproducibility (default NULL).
#'
#' @param make_plots Logical; if TRUE, write paired scatter, delta histogram,
#'   and delta ECDF per contrast.
#'
#' @param pos_threshold Numeric threshold for calling "positive selection"
#'   (genes with dN/dS > pos_threshold are counted in the paired enrichment test; default 1).
#'
#' @param drop_zero_deltas Logical; if TRUE, drop delta==0 before Wilcoxon (legacy behavior).
#'   Default FALSE (keep zeros; ties handled by wilcox.test).
#'
#' @param dedup_keys What to do if merge keys are duplicated within A or B:
#'   "error" (default) stops with a message; "first" keeps first occurrence.
#'
#' @return Invisibly, a character vector of summary TSV paths (one per contrast x side_tag).
#' @export
dnds_contrast <- function(dnds_annot_file_a = NULL,
                          dnds_annot_file_b = NULL,
                          comparison_file   = NULL,
                          contrast_file     = NULL,
                          output_dir        = getwd(),
                          regions_bed       = NULL,
                          regions_coord     = c("bed0", "gff1"),
                          region_seq_col    = NULL,
                          region_start_col  = NULL,
                          region_end_col    = NULL,
                          region_name_col   = NULL,
                          merge_cols        = NULL,
                          sides             = c("query", "subject"),
                          max_dnds          = 10,
                          ci_method         = c("normal", "bootstrap"),
                          n_boot            = 1000,
                          ci_seed           = NULL,
                          make_plots        = TRUE,
                          pos_threshold     = 1,
                          drop_zero_deltas  = FALSE,
                          dedup_keys        = c("error", "first")) {

  ci_method     <- match.arg(ci_method)
  sides         <- match.arg(sides, choices = c("query", "subject"), several.ok = TRUE)
  regions_coord <- match.arg(regions_coord)
  dedup_keys    <- match.arg(dedup_keys)

  do_regions <- !is.null(regions_bed) && nzchar(regions_bed)

  regions <- NULL
  if (do_regions) {
    regions <- .read_regions_bed(regions_bed,
                                 region_seq_col   = region_seq_col,
                                 region_start_col = region_start_col,
                                 region_end_col   = region_end_col,
                                 region_name_col  = region_name_col,
                                 regions_coord    = regions_coord,
                                 stop_on_invalid  = TRUE)
  }

  .wilcox_signed <- function(delta, drop_zero = FALSE) {
    delta <- delta[is.finite(delta)]
    if (drop_zero) delta <- delta[delta != 0]
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
        title = paste0(contrast_name, " - delta histogram (", side_tag, ")")
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
        title = paste0(contrast_name, " - delta ECDF (", side_tag, ")")
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

  .check_dups_or_handle <- function(d, by_cols, label, action = c("error", "first")) {
    action <- match.arg(action)

    # Use an ASCII separator to avoid non-ASCII/non-printable characters
    key <- do.call(paste, c(d[by_cols], sep = "||"))
    dup <- duplicated(key)
    if (any(dup)) {
      # include a few example keys for debugging
      ex_keys <- unique(key[dup])
      msg <- paste0(
        "Duplicate merge keys detected in ", label, " (", sum(dup), " duplicated rows; ",
        length(ex_keys), " duplicated key(s)). ",
        "Example duplicated key(s): ",
        paste(utils::head(ex_keys, 5), collapse = " | "),
        ".\nThis can bias results if handled silently. Resolve upstream (recommended) ",
        "or set dedup_keys='first' to keep the first occurrence."
      )
      if (action == "error") stop(msg)
      # else keep first
      d <- d[!dup, , drop = FALSE]
    }
    d
  }

  .merge_dnds <- function(dA, dB, by_cols, suffixes = c("_A", "_B")) {
    by_cols <- as.character(by_cols)

    # drop empty/NA keys
    for (cc in by_cols) {
      dA <- dA[!is.na(dA[[cc]]) & nzchar(as.character(dA[[cc]])), , drop = FALSE]
      dB <- dB[!is.na(dB[[cc]]) & nzchar(as.character(dB[[cc]])), , drop = FALSE]
    }

    # check duplicates (hard-fail by default)
    dA <- .check_dups_or_handle(dA, by_cols, label = "table A", action = dedup_keys)
    dB <- .check_dups_or_handle(dB, by_cols, label = "table B", action = dedup_keys)

    merge(dA, dB, by = by_cols, suffixes = suffixes)
  }

  .paired_pos_tests <- function(dNdS_A, dNdS_B, pos_threshold = 1) {
    # Paired binary outcomes; only keep pairs where BOTH values are finite
    ok <- is.finite(dNdS_A) & is.finite(dNdS_B)
    if (!any(ok)) {
      return(list(
        n_pairs_both_finite = 0L,
        n_pos_A = NA_integer_, n_pos_B = NA_integer_,
        b_Apos_Bnon = NA_integer_, c_Anon_Bpos = NA_integer_,
        mcnemar_p_two_sided = NA_real_,
        binom_p_A_gt_B = NA_real_
      ))
    }

    Apos <- dNdS_A[ok] > pos_threshold
    Bpos <- dNdS_B[ok] > pos_threshold

    # 2x2 for McNemar:
    #           Bpos   !Bpos
    # Apos        a      b
    # !Apos       c      d
    a <- sum(Apos & Bpos)
    b <- sum(Apos & !Bpos)  # discordant in A direction
    c <- sum(!Apos & Bpos)  # discordant in B direction
    d <- sum(!Apos & !Bpos)

    # McNemar (paired), 2-sided
    mcn_p <- tryCatch(
      stats::mcnemar.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value,
      error = function(e) NA_real_
    )

    # One-sided exact sign/binomial test on discordant pairs only:
    # H0: b and c equally likely; H1: b > c (more Apos when Bnonpos)
    bin_p <- NA_real_
    if ((b + c) > 0) {
      bin_p <- stats::binom.test(x = b, n = b + c, p = 0.5, alternative = "greater")$p.value
    }

    list(
      n_pairs_both_finite = sum(ok),
      n_pos_A = sum(Apos),
      n_pos_B = sum(Bpos),
      b_Apos_Bnon = b,
      c_Anon_Bpos = c,
      mcnemar_p_two_sided = mcn_p,
      binom_p_A_gt_B = bin_p
    )
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

    # Apply filtering (finite dNdS and below threshold)
    dA <- .filter_dnds(dA_raw, max_dnds = max_dnds)
    dB <- .filter_dnds(dB_raw, max_dnds = max_dnds)

    if (!nrow(dA) || !nrow(dB)) {
      warning("No rows after filtering for contrast ", contrast_name)
      return(NA_character_)
    }

    mode_label <- if (do_regions) paste0("regional (region-only; regions_coord=", regions_coord, ")") else "genome-wide"

    # ----------------------------
    # Regional restriction (optional)
    # ----------------------------
    if (do_regions) {
      if (!sideA %in% c("query", "subject")) stop("sideA must be 'query' or 'subject' in regional mode.")
      if (!sideB %in% c("query", "subject")) stop("sideB must be 'query' or 'subject' in regional mode.")

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
        merge_cols_use <- c("query_id", "subject_id")
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
        .require_cols(dA, c("query_id", "subject_id"), context = paste0("table A (", compA_name, ")"))
        .require_cols(dB, c("query_id", "subject_id"), context = paste0("table B (", compB_name, ")"))

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
      dcols <- grep("^dNdS", names(merged), value = TRUE)
      if (length(dcols) == 2) {
        names(merged)[match(dcols[1], names(merged))] <- "dNdS_A"
        names(merged)[match(dcols[2], names(merged))] <- "dNdS_B"
      }
    }

    if (!all(c("dNdS_A", "dNdS_B") %in% names(merged))) {
      stop("Internal error: merged table is missing dNdS_A and/or dNdS_B after merging.")
    }

    merged$delta <- merged$dNdS_A - merged$dNdS_B

    # ----------------------------
    # Stats / tests
    # ----------------------------
    stats_delta  <- .calc_mean_ci(merged$delta, ci_method = ci_method, n_boot = n_boot, seed = ci_seed)
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

    p_wilcox <- .wilcox_signed(merged$delta, drop_zero = drop_zero_deltas)

    # Paired positive-selection enrichment (McNemar + one-sided exact binomial on discordant pairs)
    pos_tests <- .paired_pos_tests(merged$dNdS_A, merged$dNdS_B, pos_threshold = pos_threshold)

    summary_df <- data.frame(
      contrast_name   = contrast_name,
      compA           = compA_name,
      compB           = compB_name,
      mode            = if (do_regions) "regional" else "global",
      regions_coord   = if (do_regions) regions_coord else NA_character_,
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
      wilcox_drop_zero_deltas = drop_zero_deltas,
      pos_threshold   = pos_threshold,
      n_pairs_both_finite_for_pos = pos_tests$n_pairs_both_finite,
      n_pos_A         = pos_tests$n_pos_A,
      n_pos_B         = pos_tests$n_pos_B,
      discordant_Apos_Bnon = pos_tests$b_Apos_Bnon,
      discordant_Anon_Bpos = pos_tests$c_Anon_Bpos,
      mcnemar_p_two_sided = pos_tests$mcnemar_p_two_sided,
      binom_p_A_gt_B_on_discordant = pos_tests$binom_p_A_gt_B,
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

    # NOTE: Keep output subdir name as "contrast" to avoid breaking existing paths/workflows.
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
