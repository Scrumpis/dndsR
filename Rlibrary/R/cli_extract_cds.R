#' @keywords internal
cli_extract_cds <- structure(
  function(...) do.call(extract_cds, list(...)),
  target = "extract_cds"
)
