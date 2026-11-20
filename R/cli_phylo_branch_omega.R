#' @keywords internal
cli_phylo_branch_omega <- structure(
  function(...) do.call(phylo_branch_omega, list(...)),
  target = "phylo_branch_omega"
)
