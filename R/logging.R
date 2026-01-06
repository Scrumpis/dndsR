#' Internal helper: run expression with stdout/stderr sunk to a log file
#' and optionally capture warnings with a summary appended to the log.
#'
#' Controlled by global options:
#'   options(dndsR.warnings = "off"|"summary"|"all")
#'   options(dndsR.warnings_max = NULL)  # NULL = no cap (print all)
#'
#' @keywords internal
.dndsr_with_log <- function(log_file,
                            tag = "dndsR",
                            header = NULL,
                            expr) {
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)

  con <- file(log_file, open = "wt")

  out_n <- sink.number(type = "output")
  msg_n <- sink.number(type = "message")

  warn_mode <- getOption("dndsR.warnings", "off")   # off|summary|all
  warn_mode <- tolower(as.character(warn_mode))[1L]
  if (is.na(warn_mode) || !nzchar(warn_mode)) warn_mode <- "off"

  # Default is "off" (no cap) -> print all warnings in summary
  warn_max_opt <- getOption("dndsR.warnings_max", NULL)
  warn_max <- suppressWarnings(as.integer(warn_max_opt))
  if (is.null(warn_max_opt) || is.na(warn_max) || warn_max < 1L) warn_max <- NA_integer_

  capture_warnings <- warn_mode %in% c("summary", "all")
  warn_buf <- character(0)

  on.exit({
    # Append warning summary at end of log (while sinks still active)
    if (capture_warnings) {
      cat("\n\n========== WARNINGS ==========\n")
      if (length(warn_buf)) {
        if (is.finite(warn_max) && length(warn_buf) > warn_max) {
          cat(sprintf("(showing first %d of %d warnings)\n\n", warn_max, length(warn_buf)))
          cat(paste0(warn_buf[seq_len(warn_max)], collapse = "\n"), "\n")
        } else {
          cat(paste0(warn_buf, collapse = "\n"), "\n")
        }
      } else {
        cat("(none)\n")
      }
      cat("======== END WARNINGS ========\n")
    }

    # unwind sinks back to original depth
    while (sink.number(type = "message") > msg_n) {
      try(sink(type = "message"), silent = TRUE)
    }
    while (sink.number(type = "output") > out_n) {
      try(sink(type = "output"), silent = TRUE)
    }
    try(close(con), silent = TRUE)
  }, add = TRUE)

  sink(con, type = "output")
  sink(con, type = "message")

  cat(sprintf("[%s] started=%s\n", tag, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  if (!is.null(header) && length(header)) {
    header <- as.character(header)
    cat(paste0(header, collapse = "\n"), "\n")
  }
  cat("\n")

  if (capture_warnings) {
    withCallingHandlers(
      { force(expr) },
      warning = function(w) {
        msg <- conditionMessage(w)
        call_txt <- conditionCall(w)
        if (!is.null(call_txt)) {
          msg <- sprintf("%s | call: %s", msg, deparse(call_txt)[1])
        }
        warn_buf <<- c(warn_buf, msg)

        # Optional live echo to console (stderr)
        if (identical(warn_mode, "all")) {
          cat(sprintf("[%s warning] %s\n", tag, conditionMessage(w)), file = stderr())
        }

        invokeRestart("muffleWarning")
      }
    )
  } else {
    # warnings behave normally (not intercepted)
    force(expr)
  }
}
