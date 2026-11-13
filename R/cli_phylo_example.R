#' @keywords internal
cli_phylo_example <- structure(
  function(...) do.call(phylo_example, list(...)),
  target = "phylo_example"
)

#' @keywords internal
cli_verify_phylo_example <- structure(
  function(...) do.call(verify_phylo_example, list(...)),
  target = "verify_phylo_example"
)
