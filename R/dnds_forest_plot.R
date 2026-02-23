#' Forest plot from dnds_contrast summary TSVs
#'
#' Reads one or more *_dnds_contrast.tsv summary files (produced by dnds_contrast)
#' and generates a forest plot of mean_delta with confidence intervals.
#'
#' Typical use: run after stats are already computed.
#'
#' @param input Path to a directory OR a character vector of TSV file paths.
#'   If a directory, files are discovered by pattern.
#' @param pattern Regex used to match summary TSV filenames.
#' @param recursive Logical; if TRUE and input is a directory, search recursively.
#' @param output_dir Where to write outputs. Default: if input is a directory, that directory;
#'   otherwise dirname of the first file.
#' @param file_prefix Prefix for output files (PDF/TSV/PNG).
#'
#' @param filter_mode Optional character vector. If provided, keep only rows where mode is in filter_mode.
#'   Example: "global" or "regional".
#' @param filter_side_tag Optional character vector. If provided, keep only rows where side_tag is in filter_side_tag.
#'   Example: "global_q_vs_q".
#'
#' @param flip_sign Logical. If TRUE, multiply delta and CI by -1 (and swap CI bounds).
#'   Use this if you want "positive = divergence in the direction B - A" for aesthetics.
#' @param x_label Label for x-axis. If NULL, chosen automatically based on flip_sign.
#'
#' @param base_family Font family name (container-safe default).
#' @param base_size Base font size.
#' @param point_size Point size.
#' @param dodge_width Dodge width for multiple comparisons per species row.
#'
#' @param species_from "compA" or "compB": which column to infer species from.
#' @param species_map Named character vector mapping raw prefixes to display names.
#' @param species_levels Optional explicit ordering of species factor (top-to-bottom).
#'
#' @param comparison_from Column to use for color grouping. Default "contrast_name".
#' @param comparison_levels Optional explicit ordering of comparisons in legend.
#'
#' @param write_png Logical; if TRUE, also write a high-DPI PNG.
#'
#' @return Invisibly returns a list with paths: plot_pdf, table_tsv, plot_png (optional), n_rows.
#' @export
dnds_forest_plot <- function(input,
                             pattern = "_dnds_contrast\\.tsv$",
                             recursive = TRUE,
                             output_dir = NULL,
                             file_prefix = "contrast_forest_plot",
                             filter_mode = NULL,
                             filter_side_tag = NULL,
                             flip_sign = FALSE,
                             x_label = NULL,
                             base_family = "Liberation Sans",
                             base_size = 15,
                             point_size = 4.5,
                             dodge_width = 0.6,
                             species_from = c("compA", "compB"),
                             species_map = c(
                               "CalbumUk" = "C. album (UK)",
                               "CalbumUK" = "C. album (UK)",
                               "CalbumUS" = "C. album (US)",
                               "Calbum"   = "C. album (US)",
                               "Cformosanum" = "C. formosanum"
                             ),
                             species_levels = NULL,
                             comparison_from = c("contrast_name"),
                             comparison_levels = NULL,
                             write_png = TRUE) {

  species_from <- match.arg(species_from)
  comparison_from <- match.arg(comparison_from)

  if (missing(input) || is.null(input)) stop("input is required (directory or TSV path(s)).")

  # ---- file discovery ----
  paths <- character(0)
  if (length(input) == 1L && dir.exists(input)) {
    paths <- list.files(input, pattern = pattern, recursive = recursive, full.names = TRUE)
    if (is.null(output_dir) || is.na(output_dir) || !nzchar(output_dir)) output_dir <- input
  } else {
    paths <- as.character(input)
    paths <- paths[file.exists(paths)]
    paths <- paths[grepl(pattern, basename(paths))]
    if (!length(paths)) stop("No existing TSVs in 'input' matching pattern: ", pattern)
    if (is.null(output_dir) || is.na(output_dir) || !nzchar(output_dir)) output_dir <- dirname(paths[1])
  }
  if (!length(paths)) stop("No summary TSVs found matching pattern: ", pattern)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  # ---- read + bind ----
  clean_colnames <- function(d) {
    nn <- names(d)
    nn <- sub("\r$", "", nn)
    nn <- trimws(nn)
    names(d) <- nn
    d
  }

  dfs <- lapply(paths, function(p) {
    d <- utils::read.table(p, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                          quote = "", comment.char = "", check.names = FALSE)
    d <- clean_colnames(d)
    d$.source_tsv <- p
    d
  })
  dd <- do.call(rbind, dfs)

  # ---- required columns ----
  need <- c("contrast_name", "compA", "compB", "mean_delta", "delta_ci_lower", "delta_ci_upper", "n")
  miss <- setdiff(need, names(dd))
  if (length(miss)) stop("Missing required column(s) in summary TSV(s): ", paste(miss, collapse = ", "))

  # optional columns
  if (!("mode" %in% names(dd))) dd$mode <- "NA"
  if (!("side_tag" %in% names(dd))) dd$side_tag <- "NA"

  # ---- filters (prevents mixing global/regional etc.) ----
  if (!is.null(filter_mode)) {
    dd <- dd[dd$mode %in% filter_mode, , drop = FALSE]
  }
  if (!is.null(filter_side_tag)) {
    dd <- dd[dd$side_tag %in% filter_side_tag, , drop = FALSE]
  }
  if (!nrow(dd)) stop("No rows left after filtering (filter_mode/filter_side_tag).")

  # ---- numeric coercion ----
  dd$mean_delta     <- suppressWarnings(as.numeric(dd$mean_delta))
  dd$delta_ci_lower <- suppressWarnings(as.numeric(dd$delta_ci_lower))
  dd$delta_ci_upper <- suppressWarnings(as.numeric(dd$delta_ci_upper))
  dd$n              <- suppressWarnings(as.integer(dd$n))

  dd <- dd[is.finite(dd$mean_delta) & is.finite(dd$delta_ci_lower) & is.finite(dd$delta_ci_upper), , drop = FALSE]
  if (!nrow(dd)) stop("No finite rows available after numeric coercion/filtering.")

  # ---- infer species from compA/compB prefix ----
  pick <- if (species_from == "compA") dd$compA else dd$compB
  species_raw <- sub("_.*$", "", as.character(pick))
  dd$species <- species_raw
  # apply mapping
  if (length(species_map)) {
    hit <- species_raw %in% names(species_map)
    dd$species[hit] <- unname(species_map[species_raw[hit]])
  }
  dd$species <- as.character(dd$species)

  # ---- comparison grouping for color ----
  dd$comparison <- as.character(dd[[comparison_from]])

  # ---- sign flipping (for prettier “positive divergence”) ----
  if (flip_sign) {
    m  <- -dd$mean_delta
    lo <- -dd$delta_ci_upper  # swap after flipping
    hi <- -dd$delta_ci_lower
    dd$mean_delta     <- m
    dd$delta_ci_lower <- lo
    dd$delta_ci_upper <- hi
  }

  # ---- factor ordering ----
  if (!is.null(species_levels)) {
    dd$species <- factor(dd$species, levels = species_levels)
  } else {
    # default: preserve the order as it appears (but stable)
    dd$species <- factor(dd$species, levels = rev(unique(dd$species)))
  }

  if (!is.null(comparison_levels)) {
    dd$comparison <- factor(dd$comparison, levels = comparison_levels)
  } else {
    dd$comparison <- factor(dd$comparison)
  }

  # ---- axis label ----
  if (is.null(x_label)) {
    x_label <- if (flip_sign) {
      # IMPORTANT: delta is now inverted; make that explicit
      expression(bold("Inverted " * Delta * " dN/dS"))
    } else {
      expression(bold(Delta * " dN/dS"))
    }
  }

  # ---- output dir + combined table ----
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_tsv <- file.path(output_dir, paste0(file_prefix, ".tsv"))
  utils::write.table(dd, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- plot ----
  pd <- ggplot2::position_dodge(width = dodge_width)

  p <- ggplot2::ggplot(dd, ggplot2::aes(x = mean_delta, y = species, color = comparison)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = delta_ci_lower, xmax = delta_ci_upper),
      position    = pd,
      width       = 0.25,
      linewidth   = 1.0,
      orientation = "y"
    ) +
    ggplot2::geom_point(position = pd, size = point_size, stroke = 1.1) +
    ggplot2::labs(x = x_label, y = NULL, color = "Comparison") +
    ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(size = base_size + 2, margin = ggplot2::margin(t = 10)),
      axis.text.x  = ggplot2::element_text(size = base_size - 1),
      axis.text.y  = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text  = ggplot2::element_text(size = base_size - 1),
      legend.position = "right"
    )

  # ---- save (PDF + optional PNG) ----
  out_pdf <- file.path(output_dir, paste0(file_prefix, ".pdf"))
  # cairo_pdf only if available; fallback otherwise
  dev_ok <- capabilities("cairo")
  if (dev_ok) {
    ggplot2::ggsave(out_pdf, p, width = 7.5, height = 4.2, device = grDevices::cairo_pdf)
  } else {
    ggplot2::ggsave(out_pdf, p, width = 7.5, height = 4.2)
  }

  out_png <- NULL
  if (isTRUE(write_png)) {
    out_png <- file.path(output_dir, paste0(file_prefix, ".png"))
    ggplot2::ggsave(out_png, p, width = 7.5, height = 4.2, dpi = 600)
  }

  message("Forest plot written: ", out_pdf)
  message("Combined table written: ", out_tsv)

  invisible(list(plot_pdf = out_pdf, table_tsv = out_tsv, plot_png = out_png, n_rows = nrow(dd)))
}
