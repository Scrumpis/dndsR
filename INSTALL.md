# Installing dndsR

| Route | Best for | Status |
|---|---|---|
| Singularity or Apptainer | Linux and HPC clusters; no root needed | Singularity tested |
| Docker | Linux, macOS, and Windows with Docker | Not tested |
| Native R | Interactive R use | Package tests pass with the image's R packages; a from-scratch native install has not been tested |

The published image holds dependencies. `tools/dndsr-install` installs dndsR from your clone into a home-directory R library, so keep the clone in place.

## Singularity or Apptainer

See the [getting started tutorial](docs/tutorials/getting-started.md). Replace `singularity` with `apptainer` if needed.

The installer writes to `~/.config/dndsr`, `~/.dndsr/rlib`, and `~/.local/bin`, and adds shell completions. For zsh users it also edits `~/.zshrc`.

The current folder is visible inside the container. For input files elsewhere, bind their location at install time and use absolute paths:

```
tools/dndsr-install install --engine singularity --sif dndsr.sif --bind /path/to/data
```

To check an installation:

```
tools/dndsr-install doctor
dndsr doctor
```

## Docker

Not yet tested. The image is `linux/amd64`, so Apple-silicon Macs run it under emulation, which can be slow on large genomes.

```
docker pull scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831
tools/dndsr-install install --engine docker --image scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831 --docker-platform linux/amd64
```

## Native R

You need R 4.5 or later (older R versions cannot install the current Bioconductor packages this needs), a C++ compiler, and DIAMOND 2.1.8 or later on your `PATH`.

DIAMOND on Linux:

```
mkdir -p "$HOME/.local/bin"
curl -fsSL https://github.com/bbuchfink/diamond/releases/download/v2.1.14/diamond-linux64.tar.gz | tar -xz -C "$HOME/.local/bin" diamond
export PATH="$HOME/.local/bin:$PATH"
```

In R. This compiles many packages and can take an hour or more:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c(
  "AnnotationDbi", "BiocGenerics", "Biostrings", "GenomeInfoDb", "GenomicFeatures",
  "Rsamtools", "S4Vectors", "txdbmaker", "pwalign", "topGO", "GO.db", "IHW", "qvalue"
), update = FALSE, ask = FALSE)
install.packages(c("cli", "data.table", "ggplot2", "jsonlite", "scales", "showtext", "svglite", "sysfonts", "xml2", "RIdeogram", "testthat", "withr"))
remotes::install_github("drostlab/orthologr@06c3098adfba1e32135d904f72ea81a76902c664", upgrade = "never")
remotes::install_github("Scrumpis/dndsR", ref = "dev", upgrade = "never")
```

Then follow [Using dndsR from R](docs/tutorials/using-r.md), starting with `verify_example()`.
