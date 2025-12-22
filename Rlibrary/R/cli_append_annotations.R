#' @keywords internal
cli_append_annotations <- structure(
  function(...) do.call(append_annotations, list(...)),
  target = "append_annotations"
)
