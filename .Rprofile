# 1) Prepend shims/ to PATH
shim_dir <- normalizePath(file.path(getwd(), "shims"), winslash = "/", mustWork = FALSE)
if (dir.exists(shim_dir)) {
  Sys.setenv(PATH = paste(shim_dir, Sys.getenv("PATH"), sep = .Platform$path.sep))
}

# 2) Use the same R user library as the CLI
rlib <- path.expand("~/.dndsr/rlib")
if (!dir.exists(rlib)) dir.create(rlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(rlib, .libPaths())))
