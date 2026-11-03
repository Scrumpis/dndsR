#' @keywords internal
cli_covariate_enrichment <- structure(
  function(...) do.call(covariate_enrichment, list(...)),
  target = "covariate_enrichment"
)
