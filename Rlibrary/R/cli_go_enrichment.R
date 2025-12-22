#' @keywords internal
cli_go_enrichment <- structure(
  function(...) do.call(go_enrichment, list(...)),
  target = "go_enrichment"
)
