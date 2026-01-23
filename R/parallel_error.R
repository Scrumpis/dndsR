#' Parallel error helper for CLI workflows
#'
#' Appends a short troubleshooting checklist to error messages when a
#' multithreaded step fails, guiding users to re-run single-threaded
#' for more informative diagnostics.
#'
#' @keywords internal
#' @noRd

append_parallel_checklist <- function(e, threads, step_name = "this step") {
  base <- conditionMessage(e)

  if (is.null(threads) || is.na(threads) || threads <= 1) {
    stop(base, call. = FALSE)
  }

  checklist <- c(
    "",
    "---- Parallel troubleshooting checklist ----",
    sprintf("Step: %s | Threads requested: %s", step_name, threads),
    "",
    "1) Re-run with --threads 1 to get a more informative error (if it is a real bug, it will reproduce unfiltered by mclapply parallelization generic error messaging).",
    "2) If --threads 1 succeeds, your previous failure was likely resource-related (RAM, open files, I/O, or a worker crash).",
    "   - Try a smaller thread count than before (e.g., 4 threads on 16 GB RAM).",
    "3) If it still fails with --threads 1, please report the single-thread error message/log output.",
    "",
    "Notes:",
    "- 'connection lost' / 'error reading from connection' often means a parallel worker died (commonly out-of-memory (OOM)).",
    "- Reducing threads often solves this issue",
    "If still getting connection lost, check input file paths (UPDATE)"
  )

  stop(paste(c(base, checklist), collapse = "\n"), call. = FALSE)
}

with_parallel_help <- function(expr, threads, step_name = "this step") {
  tryCatch(
    force(expr),
    error = function(e) append_parallel_checklist(e, threads = threads, step_name = step_name)
  )
}
