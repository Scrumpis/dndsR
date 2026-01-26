#' @keywords internal
cli_gene_state_contrast <- structure(
  function(...) do.call(gene_state_contrast, list(...)),
  target = "gene_state_contrast"
)
