# Rlibrary/R/cli_helpers/path_utils.R

# Find dndsR repo root by walking up to a DESCRIPTION file.
# Works when users clone the repo and run CLI from anywhere inside it.
.find_pkg_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)

  # If invoked via Rscript --file=..., anchor search near the script path too
  ca <- commandArgs(trailingOnly = FALSE)
  fi <- grep("^--file=", ca, value = TRUE)
  if (length(fi)) {
    script <- sub("^--file=", "", fi[1])
    if (nzchar(script)) {
      cur <- normalizePath(dirname(script), winslash = "/", mustWork = FALSE)
    }
  }

  for (k in 0:50) {
    cand <- file.path(cur, "DESCRIPTION")
    if (file.exists(cand)) return(cur)

    parent <- normalizePath(file.path(cur, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, cur)) break
    cur <- parent
  }

  NULL
}

# Resolve the dev-checkout source path for a given function name, e.g. ipr_enrichment -> Rlibrary/R/ipr_enrichment.R
.fn_source_path <- function(fn_name, pkg = "dndsR") {
  root <- .find_pkg_root()
  if (!is.null(root)) {
    p <- file.path(root, "Rlibrary", "R", paste0(fn_name, ".R"))
    if (file.exists(p)) return(p)
  }

  # Optional fallback: installed package copy (only works if you also ship R sources in extdata)
  p2 <- system.file("R", paste0(fn_name, ".R"), package = pkg)
  if (nzchar(p2) && file.exists(p2)) return(p2)

  NULL
}
