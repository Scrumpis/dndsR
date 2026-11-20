#' @keywords internal
cli_phylo_absrel <- structure(
  function(...) do.call(phylo_absrel, list(...)),
  target = "phylo_absrel"
)
