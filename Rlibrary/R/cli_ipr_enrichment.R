#' @keywords internal
cli_ipr_enrichment <- structure(
  function(...) do.call(ipr_enrichment, list(...)),
  target = "ipr_enrichment"
)
