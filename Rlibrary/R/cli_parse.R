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
