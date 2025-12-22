cli_list_commands <- function() {
  ns <- asNamespace("dndsR")
  fns <- ls(ns, pattern = "^cli_")
  sort(sub("^cli_", "", fns))
}

cli_lookup_command <- function(cmd) {
  ns <- asNamespace("dndsR")
  fn_name <- paste0("cli_", cmd)
  if (!exists(fn_name, envir = ns, inherits = FALSE)) return(NULL)
  get(fn_name, envir = ns, inherits = FALSE)
}

cli_print_usage <- function() {
  cmds <- cli_list_commands()
  cat("Usage: dndsr <command> [--key value] [--flag]\n\nCommands:\n")
  cat(paste0("  ", cmds), sep = "\n")
  cat("\nRun: dndsr <command> --help\n")
}

cli_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  suppressPackageStartupMessages({
    library(cli)
  })

  cli_bootstrap_path()

  if (length(argv) == 0 || argv[1] %in% c("-h","--help","help")) {
    cli_print_usage()
    return(invisible(NULL))
  }

  cmd <- argv[1]
  args <- argv[-1]

  if (cmd == "doctor") {
    cli_doctor()
    return(invisible(NULL))
  }

  fn <- cli_lookup_command(cmd)
  if (is.null(fn)) {
    cli::cli_alert_danger("Unknown command: {.val {cmd}}")
    cli_print_usage()
    quit(status = 1)
  }

  if (any(args %in% c("-h","--help","help"))) {
    cli_print_command_help(cmd)
    return(invisible(NULL))
  }

  parsed <- cli_parse_args(args)

  # global bridges (optional)
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
