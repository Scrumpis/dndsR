# R Installation for non-container usage in R
Current way to use in Rstudio, as container image does not contain rstudio and don't have way to use with contained dependencies without using rstudio server from container  
Will update this.
## Put your writable user library first
```
userlib <- file.path("~","Library","R","arm64","4.3","library")
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))
```

## Use Bioconductor’s repos and turn off interactive update prompts
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = userlib, dependencies = FALSE)
options(repos = BiocManager::repositories())
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
options(ask = FALSE)
```

## Install just what you need, into *userlib*, without pulling updates
```
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools", lib = userlib, dependencies = FALSE)
if (!requireNamespace("usethis", quietly = TRUE))
  install.packages("usethis", lib = userlib, dependencies = FALSE)

BiocManager::install(
  c("Biostrings","GenomicFeatures","GenomeInfoDb","Rsamtools"),
  ask = FALSE, update = FALSE
)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IHW")
```
