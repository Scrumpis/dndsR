#' @keywords internal
.dnds_cli_internal <- c(
  "bootstrap_path",
  "cast",
  "list_commands",
  "lookup_command",
  "main",
  "parse_args",
  "print_command_help",
  "print_usage"
)

#' @keywords internal
cli_list_commands <- function() {
  ns <- asNamespace("dndsR")
  fns <- ls(ns, pattern = "^cli_", all.names = TRUE)
  cmds <- sub("^cli_", "", fns)
  cmds <- setdiff(cmds, .dnds_cli_internal)
  sort(cmds)
}

#' @keywords internal
cli_lookup_command <- function(cmd) {
  ns <- asNamespace("dndsR")
  fn_name <- paste0("cli_", cmd)
  if (!exists(fn_name, envir = ns, inherits = FALSE)) return(NULL)
  get(fn_name, envir = ns, inherits = FALSE)
}

# ---- Global CLI flags footer -------------------------------------------------

#' @keywords internal
.dndsr_global_flags <- list(
  # name = c("syntax", "description")
  threads      = c("-t, --threads <int>", "Threads to use (also sets options(dndsR.threads))."),
  verbose      = c("-v, --verbose",       "Verbose logging (sets options(dndsR.verbose = TRUE))."),
  warnings     = c("--warnings <mode>",   "Warnings: off|summary|all (sets options(dndsR.warnings))."),
  warnings_max = c("--warnings-max <int>","Max warnings to print (<=0 or non-int = no cap).")
)

#' @keywords internal
cli_print_global_flags <- function() {
  cat("\nGlobal options (available for all commands):\n")
  # Nice aligned two-column print
  rows <- lapply(.dndsr_global_flags, function(x) {
    stopifnot(length(x) == 2L)
    x
  })
  rows <- do.call(rbind, rows)

  width <- max(nchar(rows[, 1, drop = TRUE]), 1L)
  for (i in seq_len(nrow(rows))) {
    cat(sprintf("  %-*s  %s\n", width, rows[i, 1], rows[i, 2]))
  }
}

#-----------------------

#' @keywords internal
cli_print_usage <- function() {
  cmds <- cli_list_commands()
  cat("Usage: dndsr <command> [--key value] [--flag]\n\nCommands:\n")
  cat(paste0("  ", cmds), sep = "\n")
  cat("\nRun: dndsr <command> --help\n")

  # Always show the global footer on top-level help too
  cli_print_global_flags()
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

    # Still print global footer even when docs are missing
    cli_print_global_flags()
    return(invisible(NULL))
  }

  txt <- utils::capture.output(utils::help(topic, package = "dndsR"))
  if (length(txt)) cat(paste0(txt, collapse = "\n"), "\n")

  # Always append footer after command help
  cli_print_global_flags()
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
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' is required for the dndsR CLI. Please install it.", call. = FALSE)
  }

  loadNamespace("dndsR")
  ns <- asNamespace("dndsR")

  if (exists("cli_bootstrap_path", envir = ns, inherits = FALSE)) {
    get("cli_bootstrap_path", envir = ns, inherits = FALSE)()
  }

  if (length(argv) > 0 && identical(argv[1], "--")) {
    argv <- argv[-1]
  }

  if (length(argv) > 0 && identical(argv[1], "__commands")) {
    cat(paste(cli_list_commands(), collapse = "\n"), "\n")
    return(invisible(NULL))
  }

  if (length(argv) == 0 || argv[1] %in% c("-h", "--help", "help")) {
    cli_print_usage()
    return(invisible(NULL))
  }

  cmd  <- argv[1]
  args <- argv[-1]

  if (length(args) == 0) {
    cli_print_command_help(cmd)
    return(invisible(NULL))
  }

  if (any(args %in% c("-h", "--help", "--usage", "help"))) {
    cli_print_command_help(cmd)
    return(invisible(NULL))
  }

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

  # Global options
  if (!is.null(parsed$threads)) options(dndsR.threads = as.integer(parsed$threads))
  if (isTRUE(parsed$verbose))   options(dndsR.verbose = TRUE)

  if (!is.null(parsed$warnings)) {
    w <- tolower(as.character(parsed$warnings))[1L]
    if (!(w %in% c("off", "summary", "all"))) {
      stop("Invalid --warnings. Use: off|summary|all", call. = FALSE)
    }
    options(dndsR.warnings = w)
  }

  if (!is.null(parsed$warnings_max)) {
    wm <- suppressWarnings(as.integer(parsed$warnings_max))
    if (is.na(wm) || wm < 1L) {
      options(dndsR.warnings_max = NULL)  # NULL = no cap (print all)
    } else {
      options(dndsR.warnings_max = wm)
    }
  }

  # Forward --threads to targets that accept it (or have ...)
  if (!is.null(parsed$threads)) {
    if (!("threads" %in% fmls) && !has_dots) {
      parsed$threads <- NULL
    }
  }

  # Only forward verbose if supported
  if (!("verbose" %in% fmls) && !has_dots) parsed$verbose <- NULL

  # Never forward warnings controls to target functions
  parsed$warnings <- NULL
  parsed$warnings_max <- NULL

  # If target does not have ..., only pass args it explicitly accepts
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
