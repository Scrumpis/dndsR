#' Apply a function across jobs, and fail if any job fails
#'
#' `parallel::mclapply()` reports worker failures as values and warnings, which a caller
#' can easily miss. This stops with a summary, so a failed job can never look like success.
#'
#' @param x Values to iterate over.
#' @param fun Function applied to each value.
#' @param threads Maximum number of forked workers. Windows always runs sequentially.
#' @param what Label used in the error message.
#' @return A list of results, in the order of `x`.
#' @keywords internal
.parallel_lapply <- function(x, fun, threads = 1L, what = "dndsR") {
  threads <- min(threads, length(x))
  if (threads <= 1L || .Platform$OS.type == "windows") {
    return(lapply(x, fun))
  }
  lost <- FALSE
  res <- withCallingHandlers(
    parallel::mclapply(x, fun, mc.cores = threads, mc.preschedule = FALSE),
    warning = function(w) {
      if (grepl("did not deliver|encountered errors|resulted in an error", conditionMessage(w))) {
        lost <<- grepl("did not deliver", conditionMessage(w)) || lost
        invokeRestart("muffleWarning")
      }
    }
  )
  .check_workers(res, what, lost)
  res
}

#' Stop if any parallel result is an error
#'
#' @param results List of worker results.
#' @param what Label used in the error message.
#' @param lost Whether a worker is known to have died without a result, for example
#'   because it ran out of memory.
#' @return `results`, invisibly.
#' @keywords internal
.check_workers <- function(results, what = "dndsR", lost = FALSE) {
  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed) || lost) {
    details <- vapply(results[failed], function(e) trimws(conditionMessage(attr(e, "condition"))), character(1))
    if (lost) details <- c(details, "a worker died without returning a result (possibly out of memory)")
    stop(
      sprintf("[%s] %d of %d jobs failed. Check the logs, or rerun with --threads 1.\n", what, sum(failed) + lost, length(results)),
      paste0("* ", unique(details), collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(results)
}
