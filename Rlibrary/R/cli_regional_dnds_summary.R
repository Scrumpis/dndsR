#' @keywords internal
cli_regional_dnds_summary <- structure(
  function(...) do.call(regional_dnds_summary, list(...)),
  target = "regional_dnds_summary"
)
