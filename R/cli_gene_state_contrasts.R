#' @keywords internal
cli_gene_state_contrasts <- structure(
  function(...) do.call(gene_state_contrasts, list(...)),
  target = "gene_state_contrasts"
)
