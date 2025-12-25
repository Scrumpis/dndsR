#' Shared CLI options for dndsR commands
#' @keywords internal
.dnds_cli_aliases <- function() {
  # map of tokens -> canonical key (snake_case)
  # include both short and long forms you want to accept
  c(
    "-c" = "comparison_file",
    "-C" = "comparison_file",
    "--comparison_file" = "comparison_file",
    "--comparison-file" = "comparison_file",

    "-o" = "output_dir",
    "-O" = "output_dir",
    "--output_dir" = "output_dir",
    "--output-dir" = "output_dir",

    "-t" = "threads",
    "-T" = "threads",
    "--threads" = "threads"
  )
}

#' Defaults for shared options
#' @keywords internal
.dnds_cli_defaults <- function() {
  list(
    output_dir = ".",
    threads = 4
  )
}

#' @keywords internal
cli_cast <- function(x) {
  xl <- tolower(trimws(x))
  if (xl %in% c("true","false")) return(xl == "true")
  if (xl %in% c("na")) return(NA)
  if (xl %in% c("null","none")) return(NULL)

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

#' Expand -C/-t style args into canonical --long_name form
#' Supports: -C file, -C=file, --comparison-file file, --comparison_file file
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
      canon <- aliases[[k]]
      if (!is.null(canon)) {
        out <- c(out, paste0("--", canon), v)
      } else if (grepl("^--", k)) {
        # keep as-is (normalize later in cli_parse_args)
        out <- c(out, k, v)
      } else {
        stop("Unknown option: ", k)
      }
      i <- i + 1
      next
    }

    # handle -C value and --long value and flags
    canon <- aliases[[a]]
    if (!is.null(canon)) {
      # is it a flag? if next token missing or another option => TRUE
      if (i == length(args) || is_opt(args[[i + 1]])) {
        out <- c(out, paste0("--", canon))
        i <- i + 1
      } else {
        out <- c(out, paste0("--", canon), args[[i + 1]])
        i <- i + 2
      }
      next
    }

    # allow canonical --something through unchanged
    if (grepl("^--", a)) {
      out <- c(out, a)
      i <- i + 1
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
parse_dnds_opts <- function(extra_aliases = NULL, defaults = NULL, args = commandArgs(trailingOnly = TRUE)) {
  aliases <- .dnds_cli_aliases()
  if (!is.null(extra_aliases)) {
    aliases <- c(aliases, extra_aliases)
  }

  # expand -C/-t etc into canonical --keys so cli_parse_args can do its thing
  expanded <- .cli_expand_aliases(args, aliases)

  opt <- cli_parse_args(expanded)

  # apply defaults (command defaults override shared defaults if provided)
  dflt <- .dnds_cli_defaults()
  if (!is.null(defaults)) dflt <- c(dflt, defaults)

  for (nm in names(dflt)) {
    if (is.null(opt[[nm]])) opt[[nm]] <- dflt[[nm]]
  }

  opt
}
