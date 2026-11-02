#' @keywords internal
cli_dnds_contrast <- structure(
  function(...) do.call(dnds_contrast, list(...)),
  target = "dnds_contrast"
)
