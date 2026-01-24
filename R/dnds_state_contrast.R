# dnds_state_contrast.R
# dN/dS state contrasts via logistic regression (optional mixed model) + forest plot
#
# Modernized to match dnds_summary / dnds_contrasts style:
# - consistent helpers (.read_ws, .clean_colnames, .filter_dnds)
# - single + batch modes
# - global-by-default; optional regional filtering via regions_bed + side
# - explicit pairwise contrasts (A vs B), no fragile coefficient-term parsing
# - optional forest plot generation built-in
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
    ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)
    if (!inherits(ok, "try-error")) keep <- keep & isTRUE(as.vector(ok))
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
  if (is.null(regions_bed)) return(NULL)
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

#' Label rows as inside/outside regions for a given side
#'
#' Adds region_status ("region"/"background") based on overlap.
#' Uses data.table::foverlaps if available; otherwise falls back to a loop.
#'
#' Requires columns like q_gff_seqname/q_gff_start/q_gff_end or s_* depending on side.
#'
#' @keywords internal
.label_region_status <- function(d, regions, side = c("query","subject")) {
  side <- match.arg(side)
  if (is.null(regions)) {
    d$region_status <- "global"
    return(d)
  }

  has_dt <- requireNamespace("data.table", quietly = TRUE)
  if (side == "query") {
    seq_col   <- "q_gff_seqname"
    start_col <- "q_gff_start"
    end_col   <- "q_gff_end"
  } else {
    seq_col   <- "s_gff_seqname"
    start_col <- "s_gff_start"
    end_col   <- "s_gff_end"
  }

  if (!all(c(seq_col, start_col, end_col) %in% names(d))) {
    stop("Missing required columns for region labeling on side='", side, "': ",
         paste(c(seq_col, start_col, end_col), collapse = ", "))
  }

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

#' Make state labels from numeric dN/dS
#'
#' @keywords internal
.make_states <- function(d, dnds_col, pos_threshold, neutral_lower, neutral_upper) {
  x <- d[[dnds_col]]
  d$state <- ifelse(x > pos_threshold, "Positive",
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
#' @keywords internal
.fit_pair <- function(d, group_col, groupA, groupB,
                      random_effect_col = NULL) {

  # keep only A and B
  d2 <- d[d[[group_col]] %in% c(groupA, groupB), , drop = FALSE]
  if (!nrow(d2)) return(NULL)

  # binary response already in d2$y
  d2[[group_col]] <- stats::relevel(factor(d2[[group_col]]), ref = groupB)

  use_re <- !is.null(random_effect_col) &&
    random_effect_col %in% names(d2) &&
    any(!is.na(d2[[random_effect_col]]))

  # model formula
  if (use_re) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("random_effect_col set but package 'lme4' is not installed.")
    }
    d3 <- d2[!is.na(d2[[random_effect_col]]), , drop = FALSE]
    if (!nrow(d3)) return(NULL)

    form <- stats::as.formula(sprintf("y ~ %s + (1|%s)", group_col, random_effect_col))
    fit  <- suppressWarnings(lme4::glmer(form, data = d3, family = stats::binomial(link = "logit")))

    coefs <- try(suppressWarnings(lme4::fixef(fit)), silent = TRUE)
    vcovm <- try(stats::vcov(fit), silent = TRUE)
    if (inherits(coefs, "try-error") || inherits(vcovm, "try-error")) return(NULL)

    term <- paste0(group_col, groupA)
    if (!term %in% names(coefs)) {
      # sometimes factor naming differs; try to find the one non-intercept term
      tt <- setdiff(names(coefs), "(Intercept)")
      if (length(tt) != 1) return(NULL)
      term <- tt
    }

    est <- unname(coefs[term])
    se  <- sqrt(diag(vcovm))[term]
    z   <- stats::qnorm(0.975)
    lo  <- est - z * se
    hi  <- est + z * se
    zval <- est / se
    pval <- 2 * stats::pnorm(-abs(zval))

    return(list(model = "glmer", log_or = est, se = unname(se),
                p_value = unname(pval), log_or_ci_lower = unname(lo), log_or_ci_upper = unname(hi)))
  }

  # GLM
  form <- stats::as.formula(sprintf("y ~ %s", group_col))
  fit  <- stats::glm(form, data = d2, family = stats::binomial(link = "logit"))
  sm   <- summary(fit)$coefficients

  # coefficient row name usually like group_colA
  rn <- rownames(sm)
  rn <- rn[rn != "(Intercept)"]
  if (!length(rn)) return(NULL)

  # ideally pick the A term; else if only one term, take it
  pick <- paste0(group_col, groupA)
  if (pick %in% rn) {
    rr <- pick
  } else if (length(rn) == 1) {
    rr <- rn[1]
  } else {
    return(NULL)
  }

  est <- sm[rr, "Estimate"]
  se  <- sm[rr, "Std. Error"]
  pval <- sm[rr, "Pr(>|z|)"]

  z <- stats::qnorm(0.975)
  lo <- est - z * se
  hi <- est + z * se

  list(model = "glm", log_or = unname(est), se = unname(se),
       p_value = unname(pval), log_or_ci_lower = unname(lo), log_or_ci_upper = unname(hi))
}

#' Forest plot for state contrasts
#'
#' @keywords internal
.plot_state_forest <- function(res,
                               out_path,
                               base_family = "Liberation Sans",
                               facet_by = c("state","scope"),
                               point_size = 2.8,
                               line_width = 0.9) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  facet_by <- match.arg(facet_by)

  # OR scale
  res$odds_ratio  <- exp(res$log_or)
  res$or_ci_lower <- exp(res$log_or_ci_lower)
  res$or_ci_upper <- exp(res$log_or_ci_upper)

  # ordering for readability within facets
  # Use the displayed label; keep in factor order of appearance
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

  if (facet_by == "state") {
    gg <- gg + ggplot2::facet_wrap(~ state, scales = "free_y")
  } else {
    gg <- gg + ggplot2::facet_grid(state ~ scope, scales = "free_y")
  }

  ggplot2::ggsave(out_path, gg, width = 11, height = 8)
  invisible(NULL)
}

# -----------------------------
# Exported function
# -----------------------------

#' dN/dS state contrasts (logistic regression; optional mixed model) + forest plot
#'
#' Classifies genes into three dN/dS "states" (Positive/Neutral/Purifying) using
#' user thresholds, then tests whether the probability of each state differs
#' between groups using pairwise logistic regression. Optionally includes a
#' random intercept via glmer (mixed logistic regression).
#'
#' Global-by-default. If `regions_bed` is provided, the analysis can be restricted
#' to genes overlapping regions on the chosen `side`. In that case, results are
#' produced for both scopes: "global" and "region" unless `scopes` is restricted.
#'
#' @param dnds_table_file Path to TSV containing dNdS and a grouping column (single mode).
#' @param comparison_file Optional batch file; if provided, reads per-comparison
#'   tables from file.path(output_dir, comparison_name, in_suffix) or default <comp>_dnds_annot.tsv.
#' @param output_dir Root directory for batch mode.
#' @param in_suffix Filename within each comparison folder. If NULL, defaults to "<comp>_dnds_annot.tsv".
#'
#' @param dnds_col Name of dN/dS column (default "dNdS").
#' @param group_col Name of grouping column (default "group").
#' @param random_effect_col Optional column for random intercept (e.g., "orthogroup") (default NULL).
#'
#' @param pos_threshold Numeric; dN/dS > this is "Positive" (default 1).
#' @param neutral_lower Lower bound for Neutral band (default 0.9).
#' @param neutral_upper Upper bound for Neutral band (default 1.1).
#'
#' @param max_dnds Drop rows with dN/dS >= max_dnds or NA (default 10).
#' @param filter_expr Optional expression evaluated in the table to filter rows.
#' @param min_n_per_group Minimum rows required per group within a scope to run comparisons (default 50).
#'
#' @param fdr_method One of "BH","BY","none" (default "BH"). Applied within each (state, scope).
#'
#' @param regions_bed Optional BED-like file. If provided, enables regional scope.
#' @param side Which side to use for region overlap if regions_bed provided ("query" or "subject").
#' @param scopes Which scopes to run: subset of c("global","region"). Default runs global always,
#'   plus region if regions_bed is provided.
#'
#' @param make_plots If TRUE, write forest plot PDF and PNG (default TRUE).
#' @param base_family Font family for plots (default "Liberation Sans").
#' @param out_prefix Output prefix (default derived from file/comp name).
#'
#' @return Invisibly returns list with:
#'   - results: data.frame of pairwise contrasts across states/scopes
#'   - paths: written TSV + plot paths
#' @export
dnds_state_contrast <- function(dnds_table_file   = NULL,
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
                                out_prefix        = NULL) {

  fdr_method <- match.arg(fdr_method)
  side       <- match.arg(side)

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

  .run_one_table <- function(d, label, out_dir) {
    if (!group_col %in% names(d)) stop("Missing group_col: ", group_col)
    if (!dnds_col %in% names(d)) stop("Missing dnds_col: ", dnds_col)

    # global filtering first (consistent with other commands)
    d <- .filter_dnds(d, dnds_col = dnds_col, filter_expr = filter_expr, max_dnds = max_dnds)
    if (!nrow(d)) {
      warning("[dnds_state_contrast] No rows after filtering for: ", label)
      return(list(results = NULL, paths = character(0)))
    }

    # add region status if possible
    d <- .label_region_status(d, regions = regions, side = side)
    d <- .make_states(d, dnds_col = dnds_col,
                      pos_threshold = pos_threshold,
                      neutral_lower = neutral_lower,
                      neutral_upper = neutral_upper)

    # compute within each scope (global, region)
    all_res <- list()

    for (scp in scopes) {
      if (scp == "global") {
        d_sc <- d
      } else {
        d_sc <- d[d$region_status == "region", , drop = FALSE]
      }

      # require enough per group
      tab <- table(d_sc[[group_col]])
      keep_groups <- names(tab)[tab >= min_n_per_group]
      d_sc <- d_sc[d_sc[[group_col]] %in% keep_groups, , drop = FALSE]

      groups <- unique(as.character(d_sc[[group_col]]))
      if (length(groups) < 2) {
        warning("[dnds_state_contrast] Not enough groups for scope='", scp, "' in: ", label)
        next
      }

      pairs <- utils::combn(sort(groups), 2, simplify = FALSE)
      states <- levels(d_sc$state)

      for (st in states) {
        # binary response for this state
        d_sc$y <- as.integer(d_sc$state == st)

        for (p in pairs) {
          A <- p[1]
          B <- p[2]

          fit <- .fit_pair(d_sc, group_col = group_col, groupA = A, groupB = B,
                           random_effect_col = random_effect_col)
          if (is.null(fit)) next

          # label: "A vs B" means OR > 1 favors A relative to B
          all_res[[length(all_res) + 1]] <- data.frame(
            label        = label,
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
      warning("[dnds_state_contrast] No models fit for: ", label)
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

    # derived OR columns
    res$odds_ratio  <- exp(res$log_or)
    res$or_ci_lower <- exp(res$log_or_ci_lower)
    res$or_ci_upper <- exp(res$log_or_ci_upper)

    # write outputs
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    prefix <- if (!is.null(out_prefix) && nzchar(out_prefix)) out_prefix else label

    out_tsv <- file.path(out_dir, sprintf("%s_dnds_state_contrasts.tsv", prefix))
    utils::write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    out_paths <- c(out_tsv)

    if (isTRUE(make_plots)) {
      # split into scope plots if both exist, or one combined faceted by scope
      out_pdf <- file.path(out_dir, sprintf("%s_dnds_state_contrasts_forest.pdf", prefix))
      out_png <- file.path(out_dir, sprintf("%s_dnds_state_contrasts_forest.png", prefix))

      # PDF
      .plot_state_forest(res, out_pdf, base_family = base_family, facet_by = "scope")
      # PNG (higher dpi)
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        # re-use plot function by saving again using last plot isn't reliable; regenerate quickly:
        .plot_state_forest(res, out_png, base_family = base_family, facet_by = "scope")
        # overwrite with better dpi if png
        # ggsave inside .plot_state_forest uses default; adjust here:
        ggplot2::ggsave(out_png, ggplot2::last_plot(), width = 11, height = 8, dpi = 600)
      }

      out_paths <- c(out_paths, out_pdf, out_png)
    }

    list(results = res, paths = out_paths)
  }

  # -------------------- Batch mode --------------------
  if (!is.null(comparison_file)) {
    cf <- .read_comparisons(comparison_file)
    outs_all <- list()
    paths_all <- character(0)

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

      d <- .read_tsv(in_file)
      run <- .run_one_table(d, label = comp, out_dir = comp_dir)
      outs_all[[comp]] <- run$results
      paths_all <- c(paths_all, run$paths)
    }

    return(invisible(list(results = outs_all, paths = paths_all)))
  }

  # -------------------- Single mode --------------------
  if (is.null(dnds_table_file)) {
    stop("Provide either dnds_table_file (single) or comparison_file (batch).")
  }
  if (!file.exists(dnds_table_file)) stop("dnds_table_file not found: ", dnds_table_file)

  d <- .read_tsv(dnds_table_file)
  label <- if (!is.null(out_prefix) && nzchar(out_prefix)) out_prefix else {
    sub("\\.tsv$", "", basename(dnds_table_file))
  }
  out_dir <- dirname(dnds_table_file)

  invisible(.run_one_table(d, label = label, out_dir = out_dir))
}
