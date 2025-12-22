#' @keywords internal
cli_dnds_ideogram <- structure(
  function(...) do.call(dnds_ideogram, list(...)),
  target = "dnds_ideogram"
)
