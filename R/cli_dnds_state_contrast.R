#' @keywords internal
cli_dnds_state_contrast <- structure(
  function(...) do.call(dnds_state_contrast, list(...)),
  target = "dnds_state_contrast"
)
