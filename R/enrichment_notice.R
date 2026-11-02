#' Calibration notice for enrichment commands that do not adjust for length or dS
#'
#' The rates come from simulated evolution, with omega independent of annotation, at one public plant
#' pair: `docs/validation/enrichment-confounding.md`. They are not predictions for other data.
#'
#' @param command `"go_enrichment"` or `"term_enrichment"`.
#' @return One string.
#' @keywords internal
.enrichment_notice <- function(command) {
  rate <- c(
    go_enrichment = "at default settings reported at least one false term in up to 43% of runs",
    term_enrichment = "Fisher's test, which this command uses, reported at least one false term in up to about 70% of runs"
  )[[command]]
  paste0("[", command, "] Calibration notice: this test does not adjust for gene length or dS. In simulated evolution ",
         "with annotation unrelated to omega, it ", rate, " at one plant pair (A. thaliana against A. lyrata), and ",
         "much less at others. Treat its terms as exploratory. covariate_enrichment adjusts for length and dS, and ",
         "has its own limits. Evidence: docs/validation/enrichment-confounding.md.")
}

#' Record the calibration notice beside an enrichment command's results
#'
#' Tables keep their columns and content. The notice and settings go in a separate file, with no
#' timestamp, so identical runs write identical files.
#'
#' @param comp_dir Folder of the results.
#' @param comp_name Comparison name.
#' @param command `"go_enrichment"` or `"term_enrichment"`.
#' @param settings Named list of the settings that define positive genes and the tests.
#' @return The path written, invisibly.
#' @keywords internal
.write_enrichment_notice <- function(comp_dir, comp_name, command, settings) {
  path <- file.path(comp_dir, sprintf("%s_%s_notice.json", comp_name, command))
  jsonlite::write_json(
    list(command = command, dndsR_version = as.character(utils::packageVersion("dndsR")),
         adjusted_for_length_and_dS = FALSE, calibration_notice = .enrichment_notice(command),
         evidence = "docs/validation/enrichment-confounding.md", settings = settings),
    path, auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  invisible(path)
}

#' Reject a threshold that a comparison with dN/dS could not use
#'
#' R compares a number with a string as text, so `dNdS > "q0.95"` is FALSE for every gene, and the
#' run looks like a real result with no positive genes.
#'
#' @param pos_threshold The value given.
#' @param command Command name, for the message.
#' @return `pos_threshold`, invisibly.
#' @keywords internal
.require_numeric_threshold <- function(pos_threshold, command) {
  if (!is.numeric(pos_threshold) || length(pos_threshold) != 1L || !is.finite(pos_threshold)) {
    stop("[", command, "] pos_threshold must be one number. Quantiles such as 'q0.95' are available in ",
         "go_enrichment and covariate_enrichment.", call. = FALSE)
  }
  invisible(pos_threshold)
}
