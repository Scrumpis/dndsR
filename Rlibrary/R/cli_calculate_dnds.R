#' @keywords internal
cli_calculate_dnds <- structure(
  function(...) do.call(calculate_dnds, list(...)),
  target = "calculate_dnds"
)
