# Rlibrary/R/cli_helpers/path_utils.R
# CLI helper to resolve R script paths

# Try to get the path of the running Rscript file (when invoked as: Rscript /path/to/cli_main.R ...)
.cli_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  fi <- grep("^--file=", ca, value = TRUE)
  if (!length(fi)) return(NULL)
  p <- sub("^--file=", "", fi[1])
  if (!nzchar(p)) return(NULL)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

# Walk upward from a directory looking for DESCRIPTION.
.walk_up_to_description <- function(start_dir, max_up = 60L) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (k in seq_len(max_up)) {
    if (file.exists(file.path(cur, "DESCRIPTION"))) return(cur)
    parent <- normalizePath(file.path(cur, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, cur)) break
    cur <- parent
  }
  NULL
}

# Resolve repo root robustly:
# Priority:
#  1) DNDSR_HOME or DNDSR_REPO (explicit)
#  2) CLI script location (Rscript --file=...)
#  3) Working directory (fallback)
.find_pkg_root <- function(start = getwd()) {
  # 1) Explicit env var (best)
  env_root <- Sys.getenv("DNDSR_HOME", Sys.getenv("DNDSR_REPO", ""))
  if (nzchar(env_root)) {
    env_root <- normalizePath(env_root, winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(env_root, "DESCRIPTION"))) return(env_root)

    # allow pointing at subdir inside repo
    hit <- .walk_up_to_description(env_root)
    if (!is.null(hit)) return(hit)
  }

  # 2) Anchor on the CLI script path if available
  sp <- .cli_script_path()
  if (!is.null(sp)) {
    hit <- .walk_up_to_description(dirname(sp))
    if (!is.null(hit)) return(hit)
  }

  # 3) Fallback: walk from start (often getwd())
  hit <- .walk_up_to_description(start)
  if (!is.null(hit)) return(hit)

  NULL
}

# Resolve the dev-checkout source file path for a given function name
.fn_source_path <- function(fn_name, pkg = "dndsR") {
  root <- .find_pkg_root()
  if (!is.null(root)) {
    p <- file.path(root, "Rlibrary", "R", paste0(fn_name, ".R"))
    if (file.exists(p)) return(p)
  }

  # Optional fallback: installed package doc/help (not source)
  # This path_utils is meant for dev-checkout; if users install the package,
  # you'd typically use help() / ?fn instead of parsing source.
  NULL
}
