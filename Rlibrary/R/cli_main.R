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

  # ---- load CLI helper utilities (works in cloned repo) ----
.local({
  # Find repo root (DESCRIPTION)
  find_root <- function(start = getwd()) {
    cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (k in 0:50) {
      if (file.exists(file.path(cur, "DESCRIPTION"))) return(cur)
      parent <- normalizePath(file.path(cur, ".."), winslash = "/", mustWork = FALSE)
      if (identical(parent, cur)) break
      cur <- parent
    }
    NULL
  }

  root <- find_root()
  if (!is.null(root)) {
    helpers <- file.path(root, "Rlibrary", "R", "cli_helpers")
    if (dir.exists(helpers)) {
      source(file.path(helpers, "path_utils.R"))
      source(file.path(helpers, "help_from_roxygen.R"))
    }
  }
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

.print_subcommand_help <- function(cmd, subcommands) {
  spec <- subcommands[[cmd]]
  fn_name <- spec$target %||% cmd

  src <- .fn_source_path(fn_name)
  if (is.null(src)) {
    stop("Can't locate source for ", fn_name, ".R.\n",
         "Expected: <repo>/Rlibrary/R/", fn_name, ".R\n",
         "Are you running from a cloned dndsR repo?")
  }
  .cli_help_from_source(fn_name, src, pkg_ns = asNamespace("dndsR"))
}

subcommands <- list(
  ipr_enrichment = list(
    opts = c(base_opts, opts_ipr_enrichment),
    fun  = run_ipr_enrichment,
    help = NULL,                 # no duplicated prose
    target = "ipr_enrichment"     # core function name
  ),
  go_enrichment = list(
    opts = c(base_opts, opts_go_enrichment),
    fun  = run_go_enrichment,
    help = NULL,
    target = "go_enrichment"
  ),
  # ... etc ...
  doctor = list(
    opts = c(base_opts),
    fun  = run_doctor,
    help = "Report backend (container vs host), shim state, and tool versions",
    target = NULL                # no roxygen-based help; keep CLI help
  )
)

# Subcommand help: dndsr <cmd> --help
if (length(args) && any(args %in% c("-h", "--help"))) {
  if (cmd %in% names(subcommands)) {
    # If this command has a target function, print roxygen-derived help; else fall back to optparse help
    if (!is.null(subcommands[[cmd]]$target)) {
      .print_subcommand_help(cmd, subcommands)
      quit(status = 0)
    }
  }
}
