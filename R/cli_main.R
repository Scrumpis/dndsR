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
    cat(
      "Tip: ensure the underlying function is documented with roxygen and exported,\n",
      "or map this CLI command to a documented topic.\n",
      sep = ""
    )
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
#' @return Invisibly returns `NULL`. Errors are signaled with `stop()`.
#'
#' @details
#' Command help is printed using installed documentation for the underlying
#' target function (e.g. `ipr_enrichment`) via `utils::help()`.
#'
#' @export
cli_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  # Do NOT attach packages inside package code
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' is required for the dndsR CLI. Please install it.", call. = FALSE)
  }

  # This function is running inside dndsR; ensure namespace is loaded so asNamespace() is valid.
  # (In practice, calling dndsR::cli_main() already loads it, but this is harmless.)
  loadNamespace("dndsR")

  ns <- asNamespace("dndsR")

  # Optional bootstrap hook (kept internal)
  if (exists("cli_bootstrap_path", envir = ns, inherits = FALSE)) {
    get("cli_bootstrap_path", envir = ns, inherits = FALSE)()
  }

  if (length(argv) > 0 && identical(argv[1], "--")) {
    argv <- argv[-1]
  }
  
  if (length(argv) == 0 || argv[1] %in% c("-h", "--help", "help")) {
    cli_print_usage()
    return(invisible(NULL))
  }

  cmd  <- argv[1]
  args <- argv[-1]

  # Special-case "doctor" if you want it independent of the usual dispatch
  if (identical(cmd, "doctor")) {
    if (exists("cli_doctor", envir = ns, inherits = FALSE)) {
      get("cli_doctor", envir = ns, inherits = FALSE)()
      return(invisible(NULL))
    }
  }

  fn <- cli_lookup_command(cmd)
  if (is.null(fn)) {
    cli::cli_alert_danger("Unknown command: {.val {cmd}}")
    cli_print_usage()
    stop(sprintf("Unknown command: %s", cmd), call. = FALSE)
  }

  if (any(args %in% c("-h", "--help", "help"))) {
    cli_print_command_help(cmd)
    return(invisible(NULL))
  }

  # Make argument parsing explicit and discoverable
  if (!exists("cli_parse_args", envir = ns, inherits = FALSE)) {
    stop("Internal error: cli_parse_args() not found in dndsR namespace.", call. = FALSE)
  }
  if (!exists("parse_dnds_opts", envir = ns, inherits = FALSE)) {
  stop("Internal error: parse_dnds_opts() not found in dndsR namespace.", call. = FALSE)
  }
  parsed <- get("parse_dnds_opts", envir = ns, inherits = FALSE)(args = args)

  # Identify underlying target function (NOT the cli_* wrapper)
  target <- attr(fn, "target", exact = TRUE)
  if (is.null(target) || !nzchar(target)) target <- cmd
  target_fn <- get(target, envir = ns, inherits = TRUE)
  fmls <- names(formals(target_fn))
  has_dots <- "..." %in% fmls
  
  # Special-case: map CLI --threads to calculate_dnds(comp_cores)
  if (identical(target, "calculate_dnds")) {
    if (!is.null(parsed$threads) && is.null(parsed$comp_cores)) {
      parsed$comp_cores <- parsed$threads
    }
  }
  
  # Global options (do not forward as args)
  if (!is.null(parsed$threads)) options(dndsR.threads = as.integer(parsed$threads))
  if (isTRUE(parsed$verbose))   options(dndsR.verbose = TRUE)
  
  # Remove globals so they never reach target functions
  parsed$threads <- NULL
  parsed$verbose <- NULL
  
  # If target does not have ..., only pass arguments it explicitly accepts
  if (!has_dots) {
    parsed <- parsed[intersect(names(parsed), fmls)]
  }

  tryCatch(
    do.call(fn, parsed),
    error = function(e) {
      cli::cli_alert_danger(conditionMessage(e))
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  invisible(NULL)
}
