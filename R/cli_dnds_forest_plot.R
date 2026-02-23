#' @keywords internal
cli_dnds_forest_plot <- structure(
  function(...) do.call(dnds_forest_plot, list(...)),
  target = "dnds_forest_plot"
)
