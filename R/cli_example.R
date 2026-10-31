#' @keywords internal
cli_example <- structure(
  function(...) do.call(dndsr_example, list(...)),
  target = "dndsr_example"
)

#' @keywords internal
cli_verify_example <- structure(
  function(...) do.call(verify_example, list(...)),
  target = "verify_example"
)
