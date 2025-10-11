# Prepend shims/ to PATH so system() finds diamond/orthofinder/etc.
shim_dir <- normalizePath(file.path(getwd(), "shims"), winslash = "/", mustWork = FALSE)
if (dir.exists(shim_dir)) {
  Sys.setenv(PATH = paste(shim_dir, Sys.getenv("PATH"), sep = .Platform$path.sep))
}

# Share the same user library as the launcher
rlib <- path.expand("~/.dndsr/rlib")
if (!dir.exists(rlib)) dir.create(rlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(rlib, .libPaths())))
