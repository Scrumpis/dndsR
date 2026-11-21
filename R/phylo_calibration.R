# Known calibration problems of the lineage-level tests, from simulation
# (docs/validation/phylo-lineage-calibration.md). They travel with the results: a `calibration`
# column on tested rows, a warning at the end of the run, and an entry in the provenance.

#' Calibration problems per method and output
#'
#' Each entry: `code` (written in the `calibration` column as `<output>:<code>`) and `note`
#' (the warning and provenance text). Outputs not listed have no known problem in the
#' simulations run; that is not evidence of calibration on real data.
#' @keywords internal
.lineage_calibration <- list(
  branch_omega = list(
    p_free_vs_m0 = list(code = "not_established_in_simulation",
      note = "p_free_vs_m0 is not established as calibrated: with taxa missing from genes (16-taxon simulation) it gave 8% false positives at 0.05 in one study and 4% in a replicate (6% of 600 genes, 95% CI 4.2-8.2%). It tests any variation among branches; per-branch omegas are descriptive estimates, not tests."),
    p_fg_omega_gt1 = list(code = "not_established_in_simulation",
      note = "p_fg_omega_gt1 is not established as calibrated at its boundary: with foreground omega exactly 1 it gave 5.5% and 6.8% false positives at 0.05 (16 and 8 taxa, 600 genes each; 6.2% of 1,200, 95% CI 4.9-7.7%). Below 1 it was conservative.")
  ),
  relax = list(
    p_value = list(code = "anticonservative_in_simulation",
      note = "RELAX's p_value is anticonservative: with no true difference it gave 8.3-9.0% false positives at 0.05 in three independent 16-taxon simulations (8.5% of 1,500 genes, 95% CI 7.2-10.1%). Neither dndsR, the simulator, nor a single starting point explained it; the likelihood-ratio statistics were larger than chi-square(1) (mean 1.28, expected 1). p_adj does not control the false-discovery rate at its nominal level.")
  )
)

#' The calibration entries that apply to a results table
#' @keywords internal
.calibration_for <- function(method, results) {
  cal <- .lineage_calibration[[method]]
  has <- function(o) o %in% names(results) && any(!is.na(results[[o]][results$status == "tested"]))
  cal <- cal[vapply(names(cal), has, logical(1))]
  if (length(cal)) cal else NULL
}

#' Add the `calibration` column to a results table, for the outputs it has
#' @keywords internal
.calibration_column <- function(results, method) {
  cal <- .calibration_for(method, results)
  if (is.null(cal)) return(results)
  tag <- paste(sprintf("%s:%s", names(cal), vapply(cal, `[[`, character(1), "code")), collapse = ";")
  results$calibration <- ifelse(results$status == "tested", tag, NA_character_)
  results
}

#' Known calibration problems, for the provenance
#' @keywords internal
.calibration_notes <- function(method, results) {
  cal <- .calibration_for(method, results)
  if (is.null(cal)) return(NULL)
  c(lapply(cal, `[[`, "note"), list(source = "docs/validation/phylo-lineage-calibration.md (simulation only)"))
}

#' Warn about known calibration problems
#' @keywords internal
.calibration_warn <- function(method, results) {
  cal <- .calibration_for(method, results)
  if (is.null(cal)) return(invisible(NULL))
  notes <- vapply(cal, `[[`, character(1), "note")
  warning(sprintf("[%s] Calibration: %s See docs/validation/phylo-lineage-calibration.md.", method, paste(notes, collapse = " ")), call. = FALSE)
}
