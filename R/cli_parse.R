#' Shared CLI options for dndsR commands
#' @keywords internal
.dnds_cli_aliases <- function() {
  # map of tokens -> canonical key (snake_case)
  # shorthand definitions:
  c(
    "-c" = "comparison_file",
    "-C" = "comparison_file",

    "-o" = "output_dir",
    "-O" = "output_dir",

    "-t" = "threads",
    "-T" = "threads",

    "-w" = "warnings",
    "-W" = "warnings"
  )
}

#' Normalize long option tokens to kebab-case.
#' Converts --foo_bar and --foo_bar=1 into --foo-bar / --foo-bar=1.
#' Leaves values untouched.
#' @keywords internal
.cli_normalize_longopts <- function(args) {
  vapply(args, function(a) {
    if (grepl("^--", a)) gsub("_", "-", a, fixed = TRUE) else a
  }, character(1))
}

#' Defaults for shared options
#' @keywords internal
.dnds_cli_defaults <- function() {
  list(
    output_dir = "."
  )
}

#' @keywords internal
cli_cast <- function(x) {
  xl <- tolower(trimws(x))
  if (xl %in% c("true","false")) return(xl == "true")
  if (xl %in% c("na")) return(NA)
  if (xl %in% c("null","none")) return(NULL)

  # comma -> character vector (nice for --sides query,subject)
  if (grepl(",", x, fixed = TRUE)) return(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))

  if (grepl("^\\d+$", xl)) return(as.integer(xl))
  if (grepl("^\\d*\\.?\\d+(e[+-]?\\d+)?$", xl, ignore.case = TRUE)) return(as.numeric(xl))
  x
}

#' @keywords internal
cli_parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (!grepl("^--", a)) stop("Unexpected token: ", a)

    key <- sub("^--", "", a)
    key <- gsub("-", "_", key)  # accept kebab-case and snake_case

    # flag?
    if (i == length(args) || grepl("^--", args[[i+1]])) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- cli_cast(args[[i+1]])
      i <- i + 2
    }
  }
  out
}

#' Safe lookup for an alias key (returns NULL if missing)
#' @keywords internal
.cli_alias_get <- function(aliases, key) {
  v <- aliases[key]
  if (length(v) && !is.na(v[[1]]) && nzchar(v[[1]])) v[[1]] else NULL
}

#' Expand -C/-t style args into canonical --long_name form
#' Supports: -C file, -C=file, --comparison-file file, --comparison_file file
#' Unknown --long options are passed through unchanged (including their values).
#' @keywords internal
.cli_expand_aliases <- function(args, aliases) {
  out <- character()
  i <- 1

  is_opt <- function(x) grepl("^-", x)

  while (i <= length(args)) {
    a <- args[[i]]

    # handle -X=value or --long=value
    if (grepl("=", a, fixed = TRUE)) {
      parts <- strsplit(a, "=", fixed = TRUE)[[1]]
      k <- parts[[1]]
      v <- paste(parts[-1], collapse = "=") # keep any extra '='

      canon <- .cli_alias_get(aliases, k)
      if (!is.null(canon)) {
        out <- c(out, paste0("--", canon), v)
      } else if (grepl("^--", k)) {
        # keep as-is; will be normalized later
        out <- c(out, k, v)
      } else {
        stop("Unknown option: ", k)
      }
      i <- i + 1
      next
    }

    # handle -C value and --long value and flags (known aliases)
    canon <- .cli_alias_get(aliases, a)
    if (!is.null(canon)) {
      if (i == length(args) || is_opt(args[[i + 1]])) {
        out <- c(out, paste0("--", canon))
        i <- i + 1
      } else {
        out <- c(out, paste0("--", canon), args[[i + 1]])
        i <- i + 2
      }
      next
    }

    # allow any unknown --long option through unchanged (and consume value if present)
    if (grepl("^--", a)) {
      if (i < length(args) && !is_opt(args[[i + 1]])) {
        out <- c(out, a, args[[i + 1]])
        i <- i + 2
      } else {
        out <- c(out, a)
        i <- i + 1
      }
      next
    }

    stop("Unexpected token: ", a)
  }

  out
}

#' Parse CLI args with shared dndsR options + per-command additions
#'
#' extra_aliases: named character vector like c("-p"="pos_threshold")
#' defaults: named list of default values (merged after parsing)
#' @keywords internal
parse_dnds_opts <- function(extra_aliases = NULL,
                            defaults = NULL,
                            args = commandArgs(trailingOnly = TRUE)) {

  aliases <- .dnds_cli_aliases()
  if (!is.null(extra_aliases)) {
    aliases <- c(aliases, extra_aliases)
  }

  # 1) Expand -C/-t/-o etc (unknown --long pass through)
  expanded <- .cli_expand_aliases(args, aliases)

  # 2) Normalize --foo_bar -> --foo-bar (ONLY for longopts)
  expanded <- .cli_normalize_longopts(expanded)

  # 3) Parse into snake_case keys
  opt <- cli_parse_args(expanded)

  # 4) Defaults (command defaults override shared defaults if provided)
  dflt <- .dnds_cli_defaults()
  if (!is.null(defaults)) dflt <- c(dflt, defaults)

  for (nm in names(dflt)) {
    if (is.null(opt[[nm]])) opt[[nm]] <- dflt[[nm]]
  }

  opt
}
