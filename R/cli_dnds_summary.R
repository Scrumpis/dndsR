#' @keywords internal
cli_dnds_summary <- structure(
  function(...) do.call(dnds_summary, list(...)),
  target = "dnds_summary"
)
