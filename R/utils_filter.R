#' Filter a dN/dS annotation table
#'
#' Internal helper used across dndsR to standardize filtering of dN/dS
#' annotation tables. Rows are retained if:
#'   1) The specified dN/dS column is finite.
#'   2) The value is strictly less than `max_dnds`.
#'   3) (Optional) The row satisfies an additional logical expression
#'      supplied via `filter_expr`.
#'
#' If `filter_expr` fails to parse/evaluate or returns a non-logical / wrong-length
#' result, it will be ignored with a warning (rather than stopping).
#'
#' @param d A data.frame containing a dN/dS column.
#' @param dnds_col Character. Name of the column containing dN/dS values.
#'   Default is `"dNdS"`.
#' @param filter_expr Optional character string representing a logical
#'   expression evaluated within the data frame (e.g. `"gene_type == 'CDS'"`).
#'   If it evaluates to length 1 it is recycled; if it evaluates to length
#'   `nrow(d)` it is used directly. Any other result or errors cause a
#'   warning and the expression is ignored.
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

  # forgiving evaluation of optional filter_expr:
  # - evaluate in the original data.frame environment (envir = d)
  # - if evaluation errors or returns non-logical / wrong length -> warn and ignore
  if (!is.null(filter_expr) && nzchar(filter_expr)) {
    ok <- try(
      eval(parse(text = filter_expr),
           envir = d,
           enclos = parent.frame()),
      silent = TRUE
    )

    if (inherits(ok, "try-error")) {
      warning("filter_expr failed to evaluate; ignoring filter_expr.")
    } else if (!is.logical(ok)) {
      warning("filter_expr did not evaluate to a logical vector; ignoring filter_expr.")
    } else if (length(ok) == 1L) {
      out_keep <- rep_len(isTRUE(ok), nrow(out))
      out <- out[out_keep, , drop = FALSE]
    } else if (length(ok) == nrow(d)) {
      # ok is defined relative to full 'd'; subset it by `keep` to align with `out`
      ok_sub <- ok[keep]
      if (length(ok_sub) != nrow(out)) {
        # defensive: if lengths don't line up, warn and ignore
        warning("filter_expr produced unexpected length after subsetting; ignoring filter_expr.")
      } else {
        out <- out[ok_sub, , drop = FALSE]
      }
    } else {
      warning("filter_expr returned logical of length ", length(ok),
              " but expected 1 or nrow(d)=", nrow(d), "; ignoring filter_expr.")
    }
  }

  out
}
