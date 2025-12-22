#' Gene-level dN/dS state contrasts across groups (aggregate across comparisons)
#'
#' This is the "gene-centric" version of dnds_state_contrast().
#' Instead of modeling per *pairwise row* (one ortholog pair), it:
#'   1) stacks one or more <comp>_dnds_annot.tsv tables,
#'   2) converts them into *gene-level* observations by taking each gene on the
#'      query and/or subject side and attaching the pairwise dN/dS value,
#'   3) aggregates across all pairwise partners per gene (median/mean),
#'   4) bins each gene into a state (Positive / Neutral / Purifying),
#'   5) runs one-vs-rest (mixed) logistic regression across subgenomes/groups.
#'
#' For Chenopodium BBCCDD this is what you want if you're trying to say:
#'   "Genes on subgenome D are more likely to be in dN/dS>1 than B or C"
#' because each gene's state is determined from its distribution across partners.
#'
#' Requires standard dndsR columns:
#'   query_id, subject_id, dNdS, q_gff_seqname, s_gff_seqname
#'
#' @param dnds_annot_files Optional character vector of one or more *_dnds_annot.tsv files.
#'   If supplied, runs in "explicit files" mode.
#' @param comparison_file Optional path (batch mode) used throughout dndsR.
#'   If provided, reads <output_dir>/<comparison_name>/<comparison_name>_dnds_annot.tsv.
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#'
#' @param focal_sides Which genes to treat as focal observations: c("query","subject") (default both).
#'   Use both for B/C/D gene-level aggregation across all comps.
#'
#' @param group_mode "subgenome" derives group from q_gff_seqname/s_gff_seqname like "Chr01C" -> "C",
#'   or "custom" uses group_col from the data you provide (after aggregation).
#' @param group_col If group_mode="custom", the group column name to use (must exist after aggregation).
#'
#' @param dnds_col Column name holding pairwise dN/dS (default "dNdS").
#' @param agg_fun Aggregation function across partners per gene: "median" (default) or "mean".
#' @param min_pairs Minimum number of pairwise observations required per gene (default 2).
#'
#' @param pos_threshold Numeric; dN/dS > pos_threshold defines "Positive" (default 1).
#' @param neutral_lower Numeric; lower bound for "Neutral" (default 0.9).
#' @param neutral_upper Numeric; upper bound for "Neutral" (default 1.1).
#'
#' @param filter_expr Optional character expression evaluated in the raw pairwise tables.
#' @param max_dnds Drop rows with dN/dS >= max_dnds or NA (default 10).
#'
#' @param random_effect_col Optional column for random intercept (e.g., gene family / orthogroup).
#'   NOTE: this column must exist in the *aggregated gene table*. If you want this,
#'   you typically merge it in beforehand or add it via your pipeline outputs.
#'
#' @param min_n_per_group Minimum number of genes per group to keep (default 200).
#' @param fdr_scope Adjust p-values within each (state, reference) block ("block")
#'   or across all tests ("global"). Default "block".
#' @param fdr_method One of "BH","BY","none". Default "BH".
#'
#' @param make_plots Logical; if TRUE, write a forest plot (default TRUE).
#' @param plot_format One of "pdf","png". Default "pdf".
#'
#' @return Invisibly:
#'   list(gene_table = ..., results_df = ..., out_tsv = ..., out_plot = ...)
#'   In batch mode (comparison_file), writes outputs to file.path(output_dir, "gene_level_state_contrast").
#' @export
gene_level_state_contrast <- function(dnds_annot_files   = NULL,
                                      comparison_file    = NULL,
                                      output_dir         = getwd(),
                                      focal_sides        = c("query", "subject"),
                                      group_mode         = c("subgenome", "custom"),
                                      group_col          = NULL,
                                      dnds_col           = "dNdS",
                                      agg_fun            = c("median", "mean"),
                                      min_pairs          = 2,
                                      pos_threshold      = 1,
                                      neutral_lower      = 0.9,
                                      neutral_upper      = 1.1,
                                      filter_expr        = NULL,
                                      max_dnds           = 10,
                                      random_effect_col  = NULL,
                                      min_n_per_group    = 200,
                                      fdr_scope          = c("block", "global"),
                                      fdr_method         = c("BH", "BY", "none"),
                                      make_plots         = TRUE,
                                      plot_format        = c("pdf", "png")) {

  focal_sides <- match.arg(focal_sides, choices = c("query", "subject"), several.ok = TRUE)
  group_mode  <- match.arg(group_mode)
  agg_fun     <- match.arg(agg_fun)
  fdr_scope   <- match.arg(fdr_scope)
  fdr_method  <- match.arg(fdr_method)
  plot_format <- match.arg(plot_format)

  # -----------------------------
  # Internal helpers (expects your standard dndsR helpers exist in the namespace):
  #   .read_comparisons, .clean_colnames, .filter_dnds
  # -----------------------------

  .safe_num <- function(x) suppressWarnings(as.numeric(x))

  .adjust <- function(p) {
    if (fdr_method == "none") return(p)
    stats::p.adjust(p, method = fdr_method)
  }

  .make_state <- function(x) {
    ifelse(x > pos_threshold, "Positive",
           ifelse(x < neutral_lower, "Purifying", "Neutral"))
  }

  .subgenome_from_seqname <- function(x) {
    x <- as.character(x)
    sg <- sub("^.*([BCD])$", "\\1", x)
    sg[!grepl("^[BCD]$", sg)] <- NA_character_
    sg
  }

  .read_one_pairwise <- function(path) {
    d <- utils::read.table(path,
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           quote = "",
                           comment.char = "",
                           check.names = FALSE)
    d <- .clean_colnames(d)

    if (!dnds_col %in% names(d)) stop("Missing dnds_col '", dnds_col, "' in: ", path)
    d$dNdS <- .safe_num(d[[dnds_col]])

    d <- .filter_dnds(d, filter_expr = filter_expr, max_dnds = max_dnds)
    d
  }

  .pairwise_to_gene_rows <- function(d, comp_name = NA_character_) {
    out <- list()

    if ("query" %in% focal_sides) {
      need <- c("query_id", "q_gff_seqname", "dNdS")
      miss <- setdiff(need, names(d))
      if (length(miss)) stop("Missing in pairwise table (query side): ", paste(miss, collapse = ", "))
      out[[length(out) + 1]] <- data.frame(
        gene_id   = as.character(d$query_id),
        side      = "query",
        seqname   = as.character(d$q_gff_seqname),
        dNdS_pair = .safe_num(d$dNdS),
        comparison = comp_name,
        stringsAsFactors = FALSE
      )
    }

    if ("subject" %in% focal_sides) {
      need <- c("subject_id", "s_gff_seqname", "dNdS")
      miss <- setdiff(need, names(d))
      if (length(miss)) stop("Missing in pairwise table (subject side): ", paste(miss, collapse = ", "))
      out[[length(out) + 1]] <- data.frame(
        gene_id   = as.character(d$subject_id),
        side      = "subject",
        seqname   = as.character(d$s_gff_seqname),
        dNdS_pair = .safe_num(d$dNdS),
        comparison = comp_name,
        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, out)
  }

  .aggregate_gene_table <- function(g) {
    g <- g[!is.na(g$gene_id) & nzchar(g$gene_id) & is.finite(g$dNdS_pair), , drop = FALSE]
    if (!nrow(g)) return(g)

    if (group_mode == "subgenome") {
      g$group <- .subgenome_from_seqname(g$seqname)
    } else {
      # group_mode == "custom": user must provide group_col in the *gene rows* (rare)
      if (is.null(group_col) || !nzchar(group_col)) stop("group_col is required when group_mode='custom'.")
      if (!group_col %in% names(g)) stop("group_col '", group_col, "' not found in gene rows.")
      g$group <- as.character(g[[group_col]])
    }

    g <- g[!is.na(g$group) & nzchar(g$group), , drop = FALSE]
    if (!nrow(g)) return(g)

    # Aggregate per gene (and group). Group should be stable; keep first non-NA.
    fun <- if (agg_fun == "median") stats::median else mean

    agg <- stats::aggregate(
      dNdS_pair ~ gene_id + group,
      data = g,
      FUN  = function(x) fun(x[is.finite(x)], na.rm = TRUE)
    )
    names(agg)[names(agg) == "dNdS_pair"] <- "dNdS_agg"

    n_pairs <- stats::aggregate(
      dNdS_pair ~ gene_id + group,
      data = g,
      FUN  = function(x) sum(is.finite(x))
    )
    names(n_pairs)[names(n_pairs) == "dNdS_pair"] <- "n_pairs"

    agg <- merge(agg, n_pairs, by = c("gene_id", "group"), all.x = TRUE)

    # Add optional random effect column if it exists in g (e.g., pre-attached family)
    if (!is.null(random_effect_col) && random_effect_col %in% names(g)) {
      re_map <- g[, c("gene_id", random_effect_col), drop = FALSE]
      re_map <- re_map[!is.na(re_map[[random_effect_col]]), , drop = FALSE]
      re_map <- re_map[!duplicated(re_map$gene_id), , drop = FALSE]
      agg <- merge(agg, re_map, by = "gene_id", all.x = TRUE)
    }

    agg <- agg[agg$n_pairs >= min_pairs & is.finite(agg$dNdS_agg), , drop = FALSE]
    if (!nrow(agg)) return(agg)

    agg$state <- .make_state(agg$dNdS_agg)
    agg
  }

  .glm_or_glmer <- function(d, ref_level, state_name) {
    d$y <- as.integer(d$state == state_name)
    d$group <- stats::relevel(factor(d$group), ref = ref_level)

    use_re <- !is.null(random_effect_col) &&
      random_effect_col %in% names(d) &&
      any(!is.na(d[[random_effect_col]])) &&
      length(unique(na.omit(d[[random_effect_col]]))) >= 2L

    if (use_re) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("random_effect_col provided but 'lme4' is not installed.")
      }
      d2 <- d[!is.na(d[[random_effect_col]]), , drop = FALSE]
      if (!nrow(d2)) return(NULL)

      form <- stats::as.formula(sprintf("y ~ group + (1|%s)", random_effect_col))
      fit  <- suppressWarnings(lme4::glmer(form, data = d2, family = stats::binomial()))

      coefs <- try(lme4::fixef(fit), silent = TRUE)
      vcovm <- try(stats::vcov(fit), silent = TRUE)
      if (inherits(coefs, "try-error") || inherits(vcovm, "try-error")) return(NULL)

      terms <- names(coefs)
      terms <- terms[terms != "(Intercept)"]
      if (!length(terms)) return(NULL)

      se  <- sqrt(diag(vcovm))[terms]
      est <- unname(coefs[terms])

      z  <- stats::qnorm(0.975)
      lo <- est - z * se
      hi <- est + z * se

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
      fit <- stats::glm(y ~ group, data = d, family = stats::binomial())
      sm  <- summary(fit)$coefficients
      sm  <- sm[rownames(sm) != "(Intercept)", , drop = FALSE]
      if (!nrow(sm)) return(NULL)

      est  <- sm[, "Estimate"]
      se   <- sm[, "Std. Error"]
      pval <- sm[, "Pr(>|z|)"]

      z  <- stats::qnorm(0.975)
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

  .decode_contrast <- function(term, ref_level) {
    # term like "groupC" etc.
    lvl <- sub("^group", "", term)
    if (!nzchar(lvl) || identical(lvl, term)) lvl <- term
    sprintf("%s vs %s", lvl, ref_level)
  }

  .forest_plot <- function(res, out_plot, plot_format) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

    res$odds_ratio  <- exp(res$log_or)
    res$or_ci_lower <- exp(res$log_or_ci_lower)
    res$or_ci_upper <- exp(res$log_or_ci_upper)

    res$contrast <- factor(res$contrast, levels = rev(unique(res$contrast)))

    gg <- ggplot2::ggplot(res, ggplot2::aes(x = odds_ratio, y = contrast)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
      ggplot2::geom_point() +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = or_ci_lower, xmax = or_ci_upper), height = 0.2) +
      ggplot2::scale_x_log10() +
      ggplot2::facet_grid(state ~ reference, scales = "free_y") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(x = "Odds ratio (gene-level, one-vs-rest)", y = NULL)

    if (plot_format == "pdf") {
      ggplot2::ggsave(out_plot, gg, width = 12, height = 8)
    } else {
      ggplot2::ggsave(out_plot, gg, width = 12, height = 8, dpi = 300)
    }
    invisible(NULL)
  }

  # -----------------------------
  # Resolve input files
  # -----------------------------
  files <- character(0)
  comp_names <- character(0)

  if (!is.null(dnds_annot_files)) {
    files <- as.character(dnds_annot_files)
    if (any(!file.exists(files))) stop("Some dnds_annot_files do not exist.")
    comp_names <- sub("_dnds_annot\\.tsv$", "", basename(files))
  } else if (!is.null(comparison_file)) {
    comps <- .read_comparisons(comparison_file)
    comp_names <- as.character(comps$comparison_name)
    files <- file.path(output_dir, comp_names, paste0(comp_names, "_dnds_annot.tsv"))
    miss <- files[!file.exists(files)]
    if (length(miss)) {
      warning("Some comparison files missing; they will be skipped:\n", paste(miss, collapse = "\n"), call. = FALSE)
    }
  } else {
    stop("Provide either dnds_annot_files OR comparison_file.")
  }

  # -----------------------------
  # Build stacked gene rows
  # -----------------------------
  gene_rows <- list()
  used_files <- character(0)

  for (i in seq_along(files)) {
    f <- files[i]
    if (!file.exists(f)) next
    cn <- if (length(comp_names) >= i) comp_names[i] else sub("_dnds_annot\\.tsv$", "", basename(f))

    d <- .read_one_pairwise(f)
    if (!nrow(d)) next

    g <- .pairwise_to_gene_rows(d, comp_name = cn)
    if (!is.null(g) && nrow(g)) {
      gene_rows[[length(gene_rows) + 1]] <- g
      used_files <- c(used_files, f)
    }
  }

  if (!length(gene_rows)) stop("No usable rows found across input files after filtering.")

  g_all <- do.call(rbind, gene_rows)

  # -----------------------------
  # Aggregate to gene-level table
  # -----------------------------
  gene_tbl <- .aggregate_gene_table(g_all)
  if (!nrow(gene_tbl)) stop("No genes remain after aggregation/min_pairs filtering.")

  # Drop underpowered groups
  tabg <- table(gene_tbl$group)
  keep_groups <- names(tabg)[tabg >= min_n_per_group]
  gene_tbl <- gene_tbl[gene_tbl$group %in% keep_groups, , drop = FALSE]
  if (length(unique(gene_tbl$group)) < 2) stop("Not enough groups after min_n_per_group filtering.")

  # -----------------------------
  # Fit one-vs-rest models with each reference
  # -----------------------------
  refs   <- unique(as.character(gene_tbl$group))
  states <- c("Positive", "Neutral", "Purifying")

  all <- list()
  for (st in states) {
    for (r in refs) {
      tmp <- .glm_or_glmer(gene_tbl, ref_level = r, state_name = st)
      if (!is.null(tmp) && nrow(tmp)) all[[length(all) + 1]] <- tmp
    }
  }
  if (!length(all)) stop("No model results produced (unexpected).")

  res <- do.call(rbind, all)

  # Decode contrasts + OR
  res$contrast <- vapply(
    seq_len(nrow(res)),
    function(i) .decode_contrast(res$term[i], res$reference[i]),
    character(1)
  )
  res$odds_ratio  <- exp(res$log_or)
  res$or_ci_lower <- exp(res$log_or_ci_lower)
  res$or_ci_upper <- exp(res$log_or_ci_upper)

  # Context columns
  res$agg_fun       <- agg_fun
  res$min_pairs     <- min_pairs
  res$pos_threshold <- pos_threshold
  res$neutral_lower <- neutral_lower
  res$neutral_upper <- neutral_upper
  res$max_dnds      <- max_dnds
  res$min_n_per_group <- min_n_per_group
  res$group_mode    <- group_mode
  res$focal_sides   <- paste(focal_sides, collapse = ",")
  res$random_effect_col <- if (!is.null(random_effect_col)) random_effect_col else NA_character_

  # P-adjust
  if (fdr_scope == "global") {
    res$p_adj <- .adjust(res$p_value)
  } else {
    res$p_adj <- NA_real_
    key <- paste(res$state, res$reference, sep = "||")
    for (k in unique(key)) {
      idx <- which(key == k)
      res$p_adj[idx] <- .adjust(res$p_value[idx])
    }
  }

  # -----------------------------
  # Write outputs
  # -----------------------------
  out_dir <- file.path(output_dir, "gene_level_state_contrast")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  tag <- paste0(
    "gene_level__",
    "sides-", gsub(",", "-", res$focal_sides[1]),
    "__agg-", agg_fun,
    "__minpairs-", min_pairs
  )

  out_tsv  <- file.path(out_dir, paste0(tag, "_results.tsv"))
  out_gene <- file.path(out_dir, paste0(tag, "_gene_table.tsv"))
  utils::write.table(res,      file = out_tsv,  sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(gene_tbl, file = out_gene, sep = "\t", quote = FALSE, row.names = FALSE)

  out_plot <- file.path(out_dir, paste0(tag, "_forest.", plot_format))
  if (isTRUE(make_plots)) .forest_plot(res, out_plot, plot_format)

  invisible(list(
    gene_table = gene_tbl,
    results_df = res,
    out_tsv    = out_tsv,
    out_gene_table_tsv = out_gene,
    out_plot   = if (make_plots) out_plot else NA_character_,
    used_files = used_files
  ))
}
