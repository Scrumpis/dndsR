#' @keywords internal
cli_bootstrap_path <- function() {
  dndsr_auto <- tolower(Sys.getenv("AUTO_SHIMS", "1"))
  run_prefix <- Sys.getenv("RUN_PREFIX", "")
  shims_dir  <- Sys.getenv("SHIMS_DIR", path.expand("~/.dndsr/shims"))

  if (dndsr_auto %in% c("1","true","yes") &&
      nzchar(run_prefix) &&
      dir.exists(shims_dir)) {
    Sys.setenv(PATH = paste(shims_dir, Sys.getenv("PATH"), sep = .Platform$path.sep))
  }

  invisible(TRUE)
}
