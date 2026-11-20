#' @keywords internal
cli_phylo_summarize <- structure(
  function(...) do.call(phylo_summarize, list(...)),
  target = "phylo_summarize"
)
