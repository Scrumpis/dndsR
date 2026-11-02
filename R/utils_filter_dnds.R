#' Filter a dN/dS annotation table
#'
#' Internal helper used across dndsR to standardize filtering of dN/dS
#' annotation tables. Rows are retained if:
#'   1) The specified dN/dS column is finite.
#'   2) The value is strictly less than `max_dnds`.
#'   3) (Optional) The row satisfies an additional logical expression
#'      supplied via `filter_expr`.
#'
#' A `filter_expr` that cannot be evaluated, or gives the wrong type or length, is an
#' error: running unfiltered could give plausible but wrong results. Rows where it gives
#' `NA` are removed.
#'
#' @param d A data.frame containing a dN/dS column.
#' @param dnds_col Character. Name of the column containing dN/dS values.
#'   Default is `"dNdS"`.
#' @param filter_expr Optional character string representing a logical
#'   expression evaluated within the data frame (e.g. `"gene_type == 'CDS'"`).
#'   If it evaluates to length 1 it is recycled; if it evaluates to length
#'   `nrow(d)` it is used directly.
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

  # coerce to numeric (characters allowed) and compute base keep mask
  x <- suppressWarnings(as.numeric(d[[dnds_col]]))
  keep <- !is.na(x) & is.finite(x) & (x < max_dnds)
  out <- d[keep, , drop = FALSE]

  if (!is.null(filter_expr) && nzchar(filter_expr)) {
    ok <- tryCatch(
      eval(parse(text = filter_expr), envir = d, enclos = parent.frame()),
      error = function(e) stop("filter_expr could not be evaluated: ", conditionMessage(e), call. = FALSE)
    )
    if (!is.logical(ok) || !(length(ok) %in% c(1L, nrow(d)))) {
      stop("filter_expr must give TRUE or FALSE for each row, or a single value.", call. = FALSE)
    }
    if (length(ok) == 1L) ok <- rep(ok, nrow(d))
    # A row is kept only if the filter is known to be TRUE for it.
    out <- d[keep & !is.na(ok) & ok, , drop = FALSE]
  }

  out
}
