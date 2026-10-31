#!/bin/bash
# RUN steps of tools/Dockerfile, in order, for building where Docker is unavailable.
# Regenerate with tools/dev/validation/test-image-build notes in docs/validation/installation.md.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive TZ=UTC LC_ALL=C.UTF-8 LANG=C.UTF-8 BSPM_DISABLE=1 MAMBA_ROOT_PREFIX=/opt/conda
export METABLASTR_REF=02dad6daf9f4a0675317177ff434ae9bbcec535f RDIAMOND_REF=d8b4006d5f2793ecf8b97a7560fa005b7a3420ee ORTHOLOGR_REF=06c3098adfba1e32135d904f72ea81a76902c664 MAMBA_VER=2.5.0 DIAMOND_VER=2.1.14
export PATH=/usr/local/bin:/opt/conda/envs/tools/bin:$PATH
set -eux; \
  apt-get update -qq; \
  apt-get install -y --no-install-recommends \
    ca-certificates gnupg2 curl wget git make g++ patch \
    pkg-config \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    libpcre2-dev libreadline-dev \
    libpng-dev libjpeg-dev libtiff5-dev libwebp-dev \
    libcairo2-dev libxt-dev libfontconfig1-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev bzip2 \
    librsvg2-dev \
  && rm -rf /var/lib/apt/lists/*
R -q -e "options(warn=2, Ncpus=parallel::detectCores()); install.packages(c('remotes','BiocManager'))"
R -q -e "options(warn=2, Ncpus=parallel::detectCores()); \
  install.packages(c( \
    'cli','data.table','lme4','xml2','ggplot2','dplyr','dtplyr','R.utils','testthat','withr', \
    'RIdeogram','doParallel','foreach','ape','Rdpack','benchmarkme','svglite', \
    'sysfonts','showtext','showtextdb' \
  )); \
  missing <- c('RIdeogram','svglite','sysfonts','showtext')[ \
    !vapply(c('RIdeogram','svglite','sysfonts','showtext'), requireNamespace, logical(1), quietly=TRUE) ]; \
  if (length(missing)) stop('Missing CRAN packages after install: ', paste(missing, collapse=', '))"
R -q -e "options(warn=2, Ncpus=parallel::detectCores()); \
  BiocManager::install(c( \
    'AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures', \
    'GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO', \
    'IHW','qvalue','rtracklayer','pwalign','txdbmaker' \
  ), update=FALSE, ask=FALSE); \
  missing <- c('AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures','GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO','IHW','qvalue','rtracklayer','pwalign','txdbmaker')[ \
    !vapply(c('AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures','GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO','IHW','qvalue','rtracklayer','pwalign','txdbmaker'), requireNamespace, logical(1), quietly=TRUE) ]; \
  if (length(missing)) stop('Missing Bioc packages after install: ', paste(missing, collapse=', '))"
R -q -e "options(warn=2); remotes::install_github(c('drostlab/metablastr@${METABLASTR_REF}','drostlab/rdiamond@${RDIAMOND_REF}','drostlab/orthologr@${ORTHOLOGR_REF}'), upgrade='never', dependencies=FALSE); \
    missing <- c('metablastr','rdiamond','orthologr')[!vapply(c('metablastr','rdiamond','orthologr'), requireNamespace, logical(1), quietly=TRUE)]; \
    if (length(missing)) stop('Missing GitHub packages after install: ', paste(missing, collapse=', '))"
Rscript -e "stopifnot( \
  requireNamespace('RIdeogram', quietly=TRUE), \
  requireNamespace('topGO', quietly=TRUE), \
  requireNamespace('GO.db', quietly=TRUE), \
  requireNamespace('GenomicFeatures', quietly=TRUE) \
)"
set -eux; \
  for f in /etc/R/Rprofile.site /usr/lib/R/etc/Rprofile.site; do \
    if [ -f "$f" ]; then \
      sed -i \
        -e '/bspm::enable/d' \
        -e '/options(bspm\./d' \
        "$f"; \
      echo 'options(bspm.enabled = FALSE)' >> "$f"; \
    fi; \
  done
set -eux; \
  curl -fsSL "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VER}" -o /tmp/micromamba.tar.bz2; \
  tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/bin --strip-components=1 bin/micromamba; \
  rm -f /tmp/micromamba.tar.bz2
micromamba create -y -n tools -c conda-forge -c bioconda \
      python=3.10.19 \
      orthofinder=2.5.4 \
      mafft=7.526 mcl=22.282 fasttree=2.2.0 \
  && micromamba clean -a -y
set -eux; \
  curl -fsSL -o /tmp/diamond.tgz \
    "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VER}/diamond-linux64.tar.gz"; \
  tar -xvzf /tmp/diamond.tgz -C /usr/local/bin diamond; \
  chmod 0755 /usr/local/bin/diamond; \
  rm -f /tmp/diamond.tgz; \
  /usr/local/bin/diamond version
R CMD INSTALL --clean /opt/dndsR \
  && install -m 0755 /opt/dndsR/inst/scripts/dndsr /usr/local/bin/dndsr \
  && Rscript -e "stopifnot(requireNamespace('dndsR', quietly=TRUE))" \
  && dndsr __commands
