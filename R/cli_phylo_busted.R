#' @keywords internal
cli_phylo_busted <- structure(
  function(...) do.call(phylo_busted, list(...)),
  target = "phylo_busted"
)
