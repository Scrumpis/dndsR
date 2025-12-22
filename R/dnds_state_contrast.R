#' dN/dS state contrasts via (mixed) logistic regression
#'
#' Compares the frequency of dN/dS "states" across groups using one-vs-rest
#' logistic regression, optionally with a random intercept (mixed model).
#'
#' States are:
#'   - Positive:  dN/dS > pos_threshold
#'   - Neutral:   neutral_lower <= dN/dS <= neutral_upper
#'   - Purifying: dN/dS < neutral_lower
#'
#' Runs all pairwise contrasts (A vs B) and returns odds ratios with CIs and FDR.
#' Produces a forest plot faceted by state and reference group.
#'
#' @param dnds_table_file Path to a TSV with dN/dS + group column (single mode).
#' @param comparison_file Optional batch file; if provided, reads per-comparison
#'   tables from file.path(output_dir, comparison_name, in_suffix).
#'   Must have a column 'comparison_name' (header or first column).
#' @param output_dir Root dir for batch outputs (default getwd()).
#' @param in_suffix Filename within each comparison folder (default "<comp>_dnds.tsv").
#'
#' @param dnds_col Name of dN/dS column (default "dNdS").
#' @param group_col Name of grouping column (e.g., "subgenome", "species", "region") (default "group").
#' @param random_effect_col Optional column for random intercept (e.g., "orthogroup") (default NULL).
#'
#' @param pos_threshold Numeric; dN/dS > this is "Positive" (default 1).
#' @param neutral_lower Numeric; lower bound for "Neutral" (default 0.9).
#' @param neutral_upper Numeric; upper bound for "Neutral" (default 1.1).
#' @param max_dnds Numeric; drop rows with dN/dS >= max_dnds or NA (default 10).
#' @param filter_expr Optional expression evaluated in the table to filter rows.
#' @param min_n_per_group Minimum rows required per group to run comparisons (default 50).
#'
#' @param fdr_method One of "BH","BY","none" (default "BH").
#' @param make_plots Logical; if TRUE, writes forest plot (default TRUE).
#' @param out_prefix Output prefix (default derived from file basename).
#'
#' @return Invisibly returns list with:
#'   - results: data.frame of contrasts
#'   - paths: written TSV + plot paths
#' @export
dnds_state_contrast <- function(dnds_table_file = NULL,
                                comparison_file = NULL,
                                output_dir = getwd(),
                                in_suffix = NULL,
                                dnds_col = "dNdS",
                                group_col = "group",
                                random_effect_col = NULL,
                                pos_threshold = 1,
                                neutral_lower = 0.9,
                                neutral_upper = 1.1,
                                max_dnds = 10,
                                filter_expr = NULL,
                                min_n_per_group = 50,
                                fdr_method = c("BH","BY","none"),
                                make_plots = TRUE,
                                out_prefix = NULL) {

  fdr_method <- match.arg(fdr_method)

  .read_ws <- function(path, header_try = TRUE) {
    utils::read.table(path, header = header_try, sep = "", quote = "\"",
                      stringsAsFactors = FALSE, comment.char = "",
                      strip.white = TRUE, blank.lines.skip = TRUE, check.names = FALSE)
  }

  .read_tsv <- function(path) {
    utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                      quote = "", comment.char = "", check.names = FALSE)
  }

  .apply_filter <- function(d) {
    if (!dnds_col %in% names(d)) stop("Missing dnds_col: ", dnds_col)
    if (!group_col %in% names(d)) stop("Missing group_col: ", group_col)

    d[[dnds_col]] <- suppressWarnings(as.numeric(d[[dnds_col]]))

    keep <- !is.na(d[[dnds_col]]) & d[[dnds_col]] < max_dnds
    if (!is.null(filter_expr) && nzchar(filter_expr)) {
      ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)
      if (!inherits(ok, "try-error")) keep <- keep & isTRUE(as.vector(ok))
    }
    d <- d[keep, , drop = FALSE]
    d
  }

  .make_states <- function(d) {
    x <- d[[dnds_col]]
    d$state <- ifelse(x > pos_threshold, "Positive",
                      ifelse(x < neutral_lower, "Purifying", "Neutral"))
    d
  }

  .pairwise_refs <- function(groups) {
    # Returns list of reference levels to relevel through (so we can extract all A vs B contrasts cleanly)
    unique(as.character(groups))
  }

  .fit_one_vs_rest <- function(d, state_name, ref_level) {
    # Build binary response
    d$y <- as.integer(d$state == state_name)
    d[[group_col]] <- stats::relevel(factor(d[[group_col]]), ref = ref_level)

    # Require lme4 only if random_effect_col provided
    use_re <- !is.null(random_effect_col) && random_effect_col %in% names(d) && any(!is.na(d[[random_effect_col]]))

    if (use_re) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("random_effect_col set but package 'lme4' is not installed.")
      }
      # drop NA in random effect for modeling
      d2 <- d[!is.na(d[[random_effect_col]]), , drop = FALSE]
      if (!nrow(d2)) return(NULL)

      form <- stats::as.formula(sprintf("y ~ %s + (1|%s)", group_col, random_effect_col))
      fit  <- suppressWarnings(lme4::glmer(form, data = d2, family = stats::binomial(link = "logit")))
      coefs <- try(suppressWarnings(lme4::fixef(fit)), silent = TRUE)
      vcovm <- try(stats::vcov(fit), silent = TRUE)

      if (inherits(coefs, "try-error") || inherits(vcovm, "try-error")) return(NULL)

      # tidy manually: group terms are like group_colLEVEL
      terms <- names(coefs)
      terms <- terms[terms != "(Intercept)"]
      if (!length(terms)) return(NULL)

      se <- sqrt(diag(vcovm))[terms]
      est <- unname(coefs[terms])

      # Wald CI
      z <- stats::qnorm(0.975)
      lo <- est - z * se
      hi <- est + z * se

      # p-values from Wald z
      zval <- est / se
      pval <- 2 * stats::pnorm(-abs(zval))

      data.frame(
        model = "glmer",
        state = state_name,
        reference = ref_level,
        term = terms,
        log_or = est,
        se = se,
        p_value = pval,
        log_or_ci_lower = lo,
        log_or_ci_upper = hi,
        stringsAsFactors = FALSE
      )
    } else {
      # simple logistic
      form <- stats::as.formula(sprintf("y ~ %s", group_col))
      fit  <- stats::glm(form, data = d, family = stats::binomial(link = "logit"))

      sm <- summary(fit)$coefficients
      sm <- sm[rownames(sm) != "(Intercept)", , drop = FALSE]
      if (!nrow(sm)) return(NULL)

      est <- sm[, "Estimate"]
      se  <- sm[, "Std. Error"]
      pval <- sm[, "Pr(>|z|)"]

      z <- stats::qnorm(0.975)
      lo <- est - z * se
      hi <- est + z * se

      data.frame(
        model = "glm",
        state = state_name,
        reference = ref_level,
        term = rownames(sm),
        log_or = unname(est),
        se = unname(se),
        p_value = unname(pval),
        log_or_ci_lower = unname(lo),
        log_or_ci_upper = unname(hi),
        stringsAsFactors = FALSE
      )
    }
  }

  .decode_term_to_contrast <- function(term, ref_level) {
    # term looks like group_colLEVEL (glm) or factor coding from glmer
    # We'll try to strip the group_col prefix.
    lvl <- sub(sprintf("^%s", group_col), "", term)
    lvl <- sub("^", "", lvl)
    lvl <- sub("^\\)", "", lvl)
    # Fallback: if nothing stripped, just keep as-is.
    if (!nzchar(lvl) || identical(lvl, term)) lvl <- term
    sprintf("%s vs %s", lvl, ref_level)
  }

  .adjust <- function(p) {
    if (fdr_method == "none") return(p)
    stats::p.adjust(p, method = fdr_method)
  }

  .plot_forest <- function(res, out_plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

    # plot on OR scale (log10 axis looks nice for OR)
    res$or <- exp(res$log_or)
    res$or_lo <- exp(res$log_or_ci_lower)
    res$or_hi <- exp(res$log_or_ci_upper)

    gg <- ggplot2::ggplot(res,
                          ggplot2::aes(x = or, y = contrast)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
      ggplot2::geom_point() +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = or_lo, xmax = or_hi), height = 0.2) +
      ggplot2::scale_x_log10() +
      ggplot2::facet_grid(state ~ reference, scales = "free_y") +
      ggplot2::labs(x = "Odds ratio (one-vs-rest)", y = NULL) +
      ggplot2::theme_bw(base_size = 12)

    ggplot2::ggsave(out_plot, gg, width = 11, height = 8)
    invisible(NULL)
  }

  .run_one <- function(d, label, out_dir) {
    d <- .apply_filter(d)
    d <- .make_states(d)

    # require enough data per group
    tab <- table(d[[group_col]])
    keep_groups <- names(tab)[tab >= min_n_per_group]
    d <- d[d[[group_col]] %in% keep_groups, , drop = FALSE]
    if (length(unique(d[[group_col]])) < 2) {
      warning("[dnds_state_contrast] Not enough groups after filtering for: ", label)
      return(list(results = NULL, paths = character(0)))
    }

    refs <- .pairwise_refs(d[[group_col]])
    states <- c("Positive","Neutral","Purifying")

    all <- list()
    for (st in states) {
      for (r in refs) {
        tmp <- .fit_one_vs_rest(d, st, r)
        if (!is.null(tmp) && nrow(tmp)) all[[length(all) + 1]] <- tmp
      }
    }
    if (!length(all)) return(list(results = NULL, paths = character(0)))

    res <- do.call(rbind, all)

    # add derived fields
    res$contrast <- vapply(seq_len(nrow(res)),
                          function(i) .decode_term_to_contrast(res$term[i], res$reference[i]),
                          character(1))
    # drop self-comparisons if any odd coding slips in
    res <- res[!grepl(sprintf("^%s vs %s$", res$reference, res$reference), res$contrast), , drop = FALSE]

    # adjust within each (state, reference) block -- makes interpretation clean
    res$p_adj <- NA_real_
    key <- paste(res$state, res$reference, sep = "||")
    for (k in unique(key)) {
      idx <- which(key == k)
      res$p_adj[idx] <- .adjust(res$p_value[idx])
    }

    # OR scale columns
    res$odds_ratio <- exp(res$log_or)
    res$or_ci_lower <- exp(res$log_or_ci_lower)
    res$or_ci_upper <- exp(res$log_or_ci_upper)

    # write outputs
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    if (is.null(out_prefix) || !nzchar(out_prefix)) out_prefix <- label

    out_tsv  <- file.path(out_dir, sprintf("%s_dnds_state_contrasts.tsv", out_prefix))
    utils::write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    out_plot <- file.path(out_dir, sprintf("%s_dnds_state_contrasts_forest.png", out_prefix))
    if (isTRUE(make_plots)) {
      # order contrasts for readability
      res2 <- res
      res2$contrast <- factor(res2$contrast, levels = rev(unique(res2$contrast)))
      .plot_forest(res2, out_plot)
    }

    list(results = res, paths = c(out_tsv, if (make_plots) out_plot))
  }

  # -------------------- Batch mode --------------------
  if (!is.null(comparison_file)) {
    cf <- try(.read_ws(comparison_file, header_try = TRUE), silent = TRUE)
    if (inherits(cf, "try-error") || !"comparison_name" %in% names(cf)) {
      cf <- .read_ws(comparison_file, header_try = FALSE)
      names(cf)[1] <- "comparison_name"
    }
    outs_all <- list()
    paths_all <- character(0)

    for (i in seq_len(nrow(cf))) {
      comp <- as.character(cf$comparison_name[i])
      comp_dir <- file.path(output_dir, comp)
      in_file <- if (!is.null(in_suffix) && nzchar(in_suffix)) {
        file.path(comp_dir, in_suffix)
      } else {
        file.path(comp_dir, sprintf("%s_dnds.tsv", comp))
      }
      if (!file.exists(in_file)) {
        warning("Missing dN/dS table: ", in_file)
        next
      }
      d <- .read_tsv(in_file)
      out_dir <- comp_dir
      out_prefix <<- comp  # per-comp naming
      run <- .run_one(d, label = comp, out_dir = out_dir)
      outs_all[[comp]] <- run$results
      paths_all <- c(paths_all, run$paths)
    }

    return(invisible(list(results = outs_all, paths = paths_all)))
  }

  # -------------------- Single mode --------------------
  if (is.null(dnds_table_file)) stop("Provide either dnds_table_file (single) or comparison_file (batch).")
  if (!file.exists(dnds_table_file)) stop("dnds_table_file not found: ", dnds_table_file)

  d <- .read_tsv(dnds_table_file)
  label <- if (!is.null(out_prefix) && nzchar(out_prefix)) out_prefix else {
    sub("\\.tsv$", "", basename(dnds_table_file))
  }
  out_dir <- dirname(dnds_table_file)

  invisible(.run_one(d, label = label, out_dir = out_dir))
}
