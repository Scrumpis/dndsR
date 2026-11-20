#' @keywords internal
cli_phylo_relax <- structure(
  function(...) do.call(phylo_relax, list(...)),
  target = "phylo_relax"
)
