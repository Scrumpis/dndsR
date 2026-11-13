#' @keywords internal
cli_phylo_prepare <- structure(
  function(...) do.call(phylo_prepare, list(...)),
  target = "phylo_prepare"
)
