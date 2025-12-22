#' @keywords internal
cli_split_comparisons <- structure(
  function(...) do.call(split_comparisons, list(...)),
  target = "split_comparisons"
)
