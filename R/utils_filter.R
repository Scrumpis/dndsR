#' Filter a dN/dS annotation table
#'
#' Internal helper used across dndsR to standardize filtering of dN/dS
#' annotation tables. Rows are retained if:
#'   1) The specified dN/dS column is finite.
#'   2) The value is strictly less than `max_dnds`.
#'   3) (Optional) The row satisfies an additional logical expression
#'      supplied via `filter_expr`.
#'
#' @param d A data.frame containing a dN/dS column.
#' @param dnds_col Character. Name of the column containing dN/dS values.
#'   Default is `"dNdS"`.
#' @param filter_expr Optional character string representing a logical
#'   expression evaluated within the filtered data frame (e.g.
#'   `"gene_type == 'CDS'"`). Must return a logical vector of length
#'   `nrow(d)` after base filtering.
#' @param max_dnds Numeric. Rows with dN/dS values greater than or equal
#'   to this threshold are removed. Default is `10`.
#' @param ... Additional arguments (currently unused; included for
#'   forward compatibility).
#'
#' @return A filtered `data.frame`.
#'
#' @keywords internal

.filter_dnds <- function(d, dnds_col = "dNdS", filter_expr = NULL, max_dnds = 10, ...) {
  if (!dnds_col %in% names(d)) {
    stop("Missing dNdS column '", dnds_col, "' in input table.")
  }

  x <- suppressWarnings(as.numeric(d[[dnds_col]]))
  keep <- is.finite(x) & x < max_dnds
  out <- d[keep, , drop = FALSE]

  if (!is.null(filter_expr)) {
    keep2 <- tryCatch(eval(parse(text = filter_expr), envir = out),
                      error = function(e) stop("filter_expr failed: ", e$message))
    if (!is.logical(keep2) || length(keep2) != nrow(out)) {
      stop("filter_expr must evaluate to a logical vector of length nrow(data).")
    }
    out <- out[keep2, , drop = FALSE]
  }

  out
}
