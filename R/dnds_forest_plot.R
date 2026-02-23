#' Forest plot from dnds_contrast summary TSVs
#'
#' Reads one or more *_dnds_contrast.tsv summary files (produced by dnds_contrast)
#' and generates a forest plot of mean_delta with confidence intervals.
#'
#' Typical use: run after stats are already computed.
#'
#' @param input Path to a directory OR a character vector of TSV file paths.
#'   If a directory, files are discovered recursively (optional) by pattern.
#' @param pattern Regex used to match summary TSV filenames.
#' @param recursive Logical; if TRUE and input is a directory, search recursively.
#' @param output_dir Where to write outputs. Default: if input is a directory, that directory;
#'   otherwise dirname of the first file.
#' @param file_prefix Prefix for output files (PDF/TSV).
#' @param facet One of "auto", "none", "mode", "side_tag", "both".
#' @param order_by Ordering of rows in the plot: "mean_delta" or "contrast_name".
#' @return Invisibly returns a list with paths: plot_pdf, table_tsv, n_rows.
#' @export
dnds_forest_plot <- function(input,
                             pattern = "_dnds_contrast\\.tsv$",
                             recursive = TRUE,
                             output_dir = NULL,
                             file_prefix = "contrast_forest_plot",
                             facet = c("auto", "none", "mode", "side_tag", "both"),
                             order_by = c("mean_delta", "contrast_name")) {

  facet <- match.arg(facet)
  order_by <- match.arg(order_by)

  if (missing(input) || is.null(input)) stop("input is required (directory or TSV path(s)).")

  # Resolve file list
  paths <- character(0)

  if (length(input) == 1L && dir.exists(input)) {
    paths <- list.files(input,
                        pattern = pattern,
                        recursive = recursive,
                        full.names = TRUE)
    if (is.null(output_dir) || is.na(output_dir) || !nzchar(output_dir)) output_dir <- input
  } else {
    paths <- as.character(input)
    paths <- paths[file.exists(paths)]
    paths <- paths[grepl(pattern, basename(paths))]
    if (is.null(output_dir) || is.na(output_dir) || !nzchar(output_dir)) {
      output_dir <- dirname(paths[1])
    }
  }

  if (!length(paths)) stop("No summary TSVs found matching pattern: ", pattern)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  # Read and combine
  dfs <- lapply(paths, function(p) {
    d <- utils::read.table(p, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                          quote = "", comment.char = "", check.names = FALSE)
    # clean colnames (CRLF artifacts etc.)
    nn <- names(d)
    nn <- sub("\r$", "", nn)
    nn <- trimws(nn)
    names(d) <- nn
    d$.source_tsv <- p
    d
  })
  dd <- do.call(rbind, dfs)

  # Validate required columns
  need <- c("contrast_name", "mean_delta", "delta_ci_lower", "delta_ci_upper", "n")
  miss <- setdiff(need, names(dd))
  if (length(miss)) stop("Missing required column(s) in summary TSV(s): ", paste(miss, collapse = ", "))

  # Optional columns used for nicer labeling/facets
  if (!("side_tag" %in% names(dd))) dd$side_tag <- "NA"
  if (!("mode" %in% names(dd))) dd$mode <- "NA"

  # Coerce numerics defensively
  dd$mean_delta     <- suppressWarnings(as.numeric(dd$mean_delta))
  dd$delta_ci_lower <- suppressWarnings(as.numeric(dd$delta_ci_lower))
  dd$delta_ci_upper <- suppressWarnings(as.numeric(dd$delta_ci_upper))
  dd$n              <- suppressWarnings(as.integer(dd$n))

  dd <- dd[is.finite(dd$mean_delta) & is.finite(dd$delta_ci_lower) & is.finite(dd$delta_ci_upper), , drop = FALSE]
  if (!nrow(dd)) stop("No finite rows available after numeric coercion/filtering.")

  # Build label: contrast + side_tag + n (keeps duplicates distinct)
  dd$label <- paste0(dd$contrast_name, " (", dd$side_tag, ", n=", dd$n, ")")

  # Order
  if (order_by == "mean_delta") {
    dd <- dd[order(dd$mean_delta, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  } else {
    dd <- dd[order(dd$contrast_name, dd$side_tag, na.last = TRUE), , drop = FALSE]
  }
  dd$label <- factor(dd$label, levels = dd$label)

  # Decide facetting
  facet_mode <- length(unique(dd$mode)) > 1
  facet_side <- length(unique(dd$side_tag)) > 1

  if (facet == "none") { facet_mode <- FALSE; facet_side <- FALSE }
  if (facet == "mode") { facet_side <- FALSE }
  if (facet == "side_tag") { facet_mode <- FALSE }
  if (facet == "both") { facet_mode <- TRUE; facet_side <- TRUE }
  # facet == "auto": keep computed booleans

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Write combined table
  out_tsv <- file.path(output_dir, paste0(file_prefix, ".tsv"))
  utils::write.table(dd, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  # Plot
  gg <- ggplot2::ggplot(dd, ggplot2::aes(y = label, x = mean_delta)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = delta_ci_lower, xmax = delta_ci_upper), height = 0) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Mean delta dN/dS (A - B) with CI", y = NULL)

  if (facet_mode && facet_side) {
    gg <- gg + ggplot2::facet_grid(mode ~ side_tag, scales = "free_y", space = "free_y")
  } else if (facet_mode) {
    gg <- gg + ggplot2::facet_wrap(~ mode, scales = "free_y")
  } else if (facet_side) {
    gg <- gg + ggplot2::facet_wrap(~ side_tag, scales = "free_y")
  }

  out_pdf <- file.path(output_dir, paste0(file_prefix, ".pdf"))
  ggplot2::ggsave(filename = out_pdf,
                  plot = gg,
                  width = 9,
                  height = max(4, 0.22 * nrow(dd)))

  invisible(list(plot_pdf = out_pdf, table_tsv = out_tsv, n_rows = nrow(dd)))
}
