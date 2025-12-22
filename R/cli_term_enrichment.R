#' @keywords internal
cli_term_enrichment <- structure(
  function(...) do.call(term_enrichment, list(...)),
  target = "term_enrichment"
)
