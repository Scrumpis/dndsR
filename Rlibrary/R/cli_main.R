#' @keywords internal
cli_list_commands <- function() {
  ns <- asNamespace("dndsR")  # safe after loadNamespace
  fns <- ls(ns, pattern = "^cli_", all.names = TRUE)
  sort(sub("^cli_", "", fns))
}

#' @keywords internal
cli_lookup_command <- function(cmd) {
  ns <- asNamespace("dndsR")
  fn_name <- paste0("cli_", cmd)
  if (!exists(fn_name, envir = ns, inherits = FALSE)) return(NULL)
  get(fn_name, envir = ns, inherits = FALSE)
}

#' @keywords internal
cli_print_usage <- function() {
  cmds <- cli_list_commands()
  cat("Usage: dndsr <command> [--key value] [--flag]\n\nCommands:\n")
  cat(paste0("  ", cmds), sep = "\n")
  cat("\nRun: dndsr <command> --help\n")
}

#' @keywords internal
cli_print_command_help <- function(cmd) {
  fn <- cli_lookup_command(cmd)
  if (is.null(fn)) {
    cat("Unknown command:", cmd, "\n\n")
    cli_print_usage()
    return(invisible(NULL))
  }

  topic <- attr(fn, "target", exact = TRUE)
  if (is.null(topic) || !nzchar(topic)) topic <- cmd

  h <- utils::help(topic, package = "dndsR")
  if (length(h) == 0) {
    cat("No documentation found for:", topic, "\n")
    cat("Tip: ensure the underlying function is documented with roxygen and exported,\n",
        "or map this CLI command to a documented topic.\n", sep = "")
    return(invisible(NULL))
  }

  utils::help(topic, package = "dndsR")  # prints help
  invisible(NULL)
}

#' dndsR command-line entrypoint
#'
#' Dispatches the `dndsr` command-line interface (CLI). This function is intended
#' to be called by the `dndsr` executable installed with the package (e.g.,
#' `inst/exec/dndsr`) or by the `dndsR-launcher`, rather than being used
#' interactively.
#'
#' The CLI discovers available commands by looking for functions named
#' `cli_<command>` in the `dndsR` namespace.
#'
#' @param argv Character vector of command-line arguments. Defaults to
#'   `commandArgs(trailingOnly = TRUE)`.
#'
#' @return Invisibly returns `NULL` on success. On failure, terminates the
#'   process with a non-zero status via `quit(status = 1)`.
#'
#' @details
#' Command help is printed using installed documentation for the underlying
#' target function (e.g. ipr_enrichment) via utils::help()
#'
#' @export
cli_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  suppressPackageStartupMessages({
    library(cli)
  })

  # Ensure package is installed/loaded
  if (!requireNamespace("dndsR", quietly = TRUE)) {
    cli::cli_alert_danger("dndsR is not installed. Run: dndsR-launcher install")
    quit(status = 1)
  }
  loadNamespace("dndsR")  # make asNamespace valid

  # Optional: bootstrap PATH/env (but this should be a normal package function)
  if (exists("cli_bootstrap_path", envir = asNamespace("dndsR"), inherits = FALSE)) {
    get("cli_bootstrap_path", envir = asNamespace("dndsR"))()
  }

  if (length(argv) == 0 || argv[1] %in% c("-h", "--help", "help")) {
    cli_print_usage()
    return(invisible(NULL))
  }

  cmd  <- argv[1]
  args <- argv[-1]

  if (cmd == "doctor") {
    if (exists("cli_doctor", envir = asNamespace("dndsR"), inherits = FALSE)) {
      get("cli_doctor", envir = asNamespace("dndsR"))()
      return(invisible(NULL))
    }
  }

  fn <- cli_lookup_command(cmd)
  if (is.null(fn)) {
    cli::cli_alert_danger("Unknown command: {.val {cmd}}")
    cli_print_usage()
    quit(status = 1)
  }

  if (any(args %in% c("-h", "--help", "help"))) {
    cli_print_command_help(cmd)
    return(invisible(NULL))
  }

  parsed <- cli_parse_args(args)  # must exist in package namespace or global

  if (!is.null(parsed$threads)) options(dndsR.threads = as.integer(parsed$threads))
  if (isTRUE(parsed$verbose))   options(dndsR.verbose = TRUE)

  tryCatch(
    do.call(fn, parsed),
    error = function(e) {
      cli::cli_alert_danger(conditionMessage(e))
      quit(status = 1)
    }
  )
}
