# HyPhy engine layer. Method-agnostic: finding the binary, running one analysis in its own
# folder with an explicit seed and thread count, and reading its JSON. Method-specific
# arguments and parsing live in the registry `.hyphy_methods`.

# HyPhy release this code was validated against. Others run, with a notice.
.hyphy_validated_version <- "2.5.103"

#' Find the HyPhy executable
#'
#' Order: the `hyphy` argument, `options(dndsR.hyphy)`, the `DNDSR_HYPHY` environment
#' variable, then `hyphy` on `PATH`.
#' @keywords internal
.hyphy_path <- function(hyphy = NULL) {
  candidates <- c(hyphy, getOption("dndsR.hyphy"), Sys.getenv("DNDSR_HYPHY", ""), Sys.which("hyphy"))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  stop("HyPhy was not found. Install it (conda install -c conda-forge -c bioconda hyphy=", .hyphy_validated_version,
       "), put 'hyphy' on PATH, or give its path with --hyphy or DNDSR_HYPHY. Published dndsR images up to v0.4.0 do not include it.",
       call. = FALSE)
}

#' HyPhy version, such as "2.5.103"
#' @keywords internal
.hyphy_version <- function(bin) {
  out <- tryCatch(suppressWarnings(system2(bin, "--version", stdout = TRUE, stderr = TRUE)), error = function(e) character())
  v <- regmatches(out, regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", out))
  if (!length(v)) stop("Could not read the HyPhy version from '", bin, " --version'.", call. = FALSE)
  v[1]
}

#' Run one HyPhy analysis in its own folder
#'
#' HyPhy writes `errors.log` into its working folder, so each run gets one. The random seed
#' and thread count are always set: in testing on HyPhy 2.5.103, a fixed seed and thread
#' count gave bit-identical results, and anything else did not.
#'
#' @param bin HyPhy executable.
#' @param method HyPhy analysis name, such as "busted".
#' @param args Character vector of keyword arguments, already quoted where needed.
#' @param workdir Folder for the run.
#' @param cpu Threads for HyPhy.
#' @param seed Random seed.
#' @return A list with `status` (exit code), `log` (path), `error` (text, or `NA`), and
#'   `command` (as run).
#' @keywords internal
.hyphy_run <- function(bin, method, args, workdir, cpu, seed) {
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(workdir, "hyphy.log")
  err_file <- file.path(workdir, "errors.log")
  if (file.exists(err_file)) unlink(err_file)
  full <- c(sprintf("CPU=%d", as.integer(cpu)), sprintf("ENV=RANDOM_SEED=%d", as.integer(seed)), method, args)
  old <- setwd(workdir)
  on.exit(setwd(old), add = TRUE)
  status <- suppressWarnings(system2(bin, full, stdout = "hyphy.log", stderr = "hyphy.log"))
  error <- NA_character_
  if (file.exists("errors.log")) {
    e <- trimws(readLines("errors.log", warn = FALSE))
    e <- sub("^Error:", "", e[nzchar(e)])
    error <- substr(paste(utils::head(e, 2), collapse = " "), 1, 300)
  }
  list(status = status, log = log, error = error, command = paste(c(shQuote(bin), full), collapse = " "))
}

#' Read HyPhy JSON and check the keys a method needs
#' @keywords internal
.hyphy_read_json <- function(path, required) {
  if (!file.exists(path)) stop("HyPhy wrote no JSON")
  j <- tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(j)) stop("HyPhy JSON could not be parsed")
  missing <- setdiff(required, names(j))
  if (length(missing)) stop("HyPhy JSON lacks: ", paste(missing, collapse = ", "))
  j
}

#' Rate classes from a HyPhy "Rate Distributions" entry
#' @keywords internal
.hyphy_rate_classes <- function(x, prefix) {
  if (is.null(x)) return(list())
  x <- x[order(as.integer(names(x)))]
  out <- list()
  for (k in seq_along(x)) {
    out[[sprintf("%s_omega_%d", prefix, k)]] <- as.numeric(x[[k]]$omega)
    out[[sprintf("%s_prop_%d", prefix, k)]] <- as.numeric(x[[k]]$proportion)
  }
  out
}

#' Registry of supported HyPhy methods
#'
#' Each entry gives the HyPhy analysis name, its settings with defaults equal to HyPhy's,
#' a function building keyword arguments, the JSON keys required, a parser returning a
#' one-row list of results, and the unit of testing ("gene" or "branch").
#' @keywords internal
.hyphy_methods <- list(
  busted = list(
    hyphy_name = "busted",
    unit = "gene",
    defaults = list(srv = "Yes", multiple_hits = "None", rates = 3L, syn_rates = 3L, error_sink = "No",
                    grid_size = 250L, starting_points = 1L, genetic_code = "Universal"),
    choices = list(srv = c("Yes", "No", "Branch-site"), multiple_hits = c("None", "Double", "Double+Triple"),
                   error_sink = c("No", "Yes"), genetic_code = "Universal"),
    args = function(s, alignment, tree, output) {
      a <- c("--alignment", shQuote(alignment), "--tree", shQuote(tree), "--branches", "FG",
             "--code", s$genetic_code, "--srv", s$srv, "--rates", s$rates)
      if (s$srv != "No") a <- c(a, "--syn-rates", s$syn_rates)
      c(a, "--multiple-hits", shQuote(s$multiple_hits), "--error-sink", s$error_sink,
        "--grid-size", s$grid_size, "--starting-points", s$starting_points, "--output", shQuote(output))
    },
    required = c("test results", "fits", "tested", "branch attributes", "input"),
    parse = function(j) {
      tr <- j[["test results"]]
      unc <- j$fits[["Unconstrained model"]]
      con <- j$fits[["Constrained model"]]
      if (is.null(tr$LRT) || is.null(tr[["p-value"]]) || is.null(unc)) stop("HyPhy JSON lacks the BUSTED test")
      ba <- j[["branch attributes"]][["0"]]
      lengths <- vapply(ba, function(b) if (is.null(b$unconstrained)) NA_real_ else as.numeric(b$unconstrained), numeric(1))
      c(list(
        lrt = as.numeric(tr$LRT),
        p_value = as.numeric(tr[["p-value"]]),
        constrained_fitted = !is.null(con),
        lnl_unconstrained = as.numeric(unc[["Log Likelihood"]]),
        lnl_constrained = if (is.null(con)) NA_real_ else as.numeric(con[["Log Likelihood"]]),
        tree_length = sum(lengths, na.rm = TRUE)
      ), .hyphy_rate_classes(unc[["Rate Distributions"]]$Test, "fg"),
      .hyphy_rate_classes(unc[["Rate Distributions"]]$Background, "bg"))
    }
  )
)

#' Branches HyPhy reports as tested
#' @keywords internal
.hyphy_tested_branches <- function(j) {
  t <- j$tested[["0"]]
  sort(names(t)[vapply(t, function(x) identical(tolower(x), "test"), logical(1))], method = "radix")
}

#' Resolve method settings against the registry
#' @keywords internal
.hyphy_settings <- function(method, given) {
  m <- .hyphy_methods[[method]]
  s <- m$defaults
  for (nm in names(given)) if (!is.null(given[[nm]])) s[[nm]] <- given[[nm]]
  for (nm in names(m$choices)) {
    v <- as.character(s[[nm]])
    hit <- m$choices[[nm]][tolower(m$choices[[nm]]) == tolower(v)]
    if (length(hit) != 1L) {
      stop(sprintf("Invalid %s '%s'. Use one of: %s.", nm, v, paste(m$choices[[nm]], collapse = ", ")), call. = FALSE)
    }
    s[[nm]] <- hit
  }
  for (nm in c("rates", "syn_rates", "grid_size", "starting_points")) {
    if (!is.null(s[[nm]])) {
      v <- suppressWarnings(as.integer(s[[nm]]))
      if (is.na(v) || v < 1L) stop("Invalid ", nm, ": ", s[[nm]], call. = FALSE)
      s[[nm]] <- v
    }
  }
  s
}
