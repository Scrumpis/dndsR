# Rlibrary/R/help_roxygen_cli.R
# CLI helper to get roxygen help info

.read_roxygen_block <- function(path, fn_name) {
  lines <- readLines(path, warn = FALSE)

  rx_fn <- paste0("^\\s*", fn_name, "\\s*<-\\s*function\\s*\\(")
  i_fn <- grep(rx_fn, lines)
  if (!length(i_fn)) return(NULL)
  i_fn <- i_fn[1]

  # Walk upward collecting "#'" lines
  i <- i_fn - 1
  block <- character(0)
  while (i >= 1 && grepl("^\\s*#'", lines[i])) {
    block <- c(lines[i], block)
    i <- i - 1
  }
  if (!length(block)) return(NULL)

  # Strip "#'" prefix
  sub("^\\s*#'\\s?", "", block)
}

.parse_roxygen <- function(block) {
  if (is.null(block) || !length(block)) return(NULL)

  is_tag <- grepl("^@", block)
  nonempty <- nzchar(trimws(block))

  # Title: first non-empty non-tag line
  title_idx <- which(nonempty & !is_tag)
  title <- if (length(title_idx)) block[title_idx[1]] else ""

  # Description: lines until first tag, excluding title line
  first_tag <- which(is_tag)
  end_desc <- if (length(first_tag)) first_tag[1] - 1 else length(block)
  desc_lines <- block[seq_len(end_desc)]
  if (length(title_idx) && title_idx[1] <= length(desc_lines)) {
    desc_lines <- desc_lines[-title_idx[1]]
  }
  desc <- paste(trimws(desc_lines[nzchar(trimws(desc_lines))]), collapse = "\n")

  # Params: supports multi-line continuation, and "@param a,b ..." splits to both
  params <- list()
  i <- 1
  while (i <= length(block)) {
    ln <- block[i]
    if (grepl("^@param\\s+", ln)) {
      rest <- sub("^@param\\s+", "", ln)
      nm <- sub("\\s+.*$", "", rest)
      txt <- trimws(sub(paste0("^", nm, "\\s*"), "", rest))

      # Continuation lines until next tag
      j <- i + 1
      cont <- character(0)
      while (j <= length(block) && !grepl("^@", block[j])) {
        if (nzchar(trimws(block[j]))) cont <- c(cont, trimws(block[j]))
        j <- j + 1
      }
      if (length(cont)) txt <- paste(c(txt, cont), collapse = " ")

      # Split comma-delimited param names, e.g. "x_axis_min,x_axis_max"
      nms <- trimws(strsplit(nm, ",", fixed = TRUE)[[1]])
      for (one in nms) params[[one]] <- txt

      i <- j
      next
    }
    i <- i + 1
  }

  list(title = title, description = desc, params = params)
}

# Print CLI help for a function using:
#  - roxygen from source file
#  - defaults + arg list from formals()
.cli_help_from_source <- function(fn_name, source_file, pkg_ns = asNamespace("dndsR")) {
  fn <- get(fn_name, envir = pkg_ns)

  block <- .read_roxygen_block(source_file, fn_name)
  rx <- .parse_roxygen(block)

  fm <- formals(fn)
  arg_names <- names(fm)

  cat("\n", fn_name, "\n", sep = "")
  if (!is.null(rx) && nzchar(rx$title)) cat(rx$title, "\n")
  if (!is.null(rx) && nzchar(rx$description)) cat("\n", rx$description, "\n", sep = "")

  cat("\nUsage:\n  dndsr ", fn_name, " [options]\n\n", sep = "")

  cat("Options:\n")
  for (nm in arg_names) {
    flag <- paste0("--", gsub("_", "-", nm))
    def <- fm[[nm]]

    def_str <- paste(deparse(def), collapse = " ")
    def_str <- gsub("\\s+", " ", def_str)

    help <- if (!is.null(rx)) rx$params[[nm]] else NULL
    if (is.null(help)) help <- ""

    if (nzchar(help)) {
      cat(sprintf("  %-26s %s", flag, help))
    } else {
      cat(sprintf("  %-26s", flag))
    }

    # Add default
    cat(sprintf(" [default: %s]\n", def_str))
  }

  cat("\nNotes:\n  - Flags accept kebab-case and snake_case (e.g. --output-dir == --output_dir)\n\n")
  invisible(TRUE)
}
