# ------------------------------------------------------------
# Base: R 4.3.3 (Debian Bookworm)
# ------------------------------------------------------------
FROM rocker/r-ver:4.3.3 AS rbase

ENV DEBIAN_FRONTEND=noninteractive TZ=UTC \
    LC_ALL=C.UTF-8 LANG=C.UTF-8

# System libs for CRAN/BioC builds
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates gnupg curl wget git make g++ \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    libpcre2-dev libreadline-dev \
    libpng-dev libjpeg-dev libtiff5-dev \
    libcairo2-dev libxt-dev libfontconfig1-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev \
 && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# Micromamba for OrthoFinder/DIAMOND toolchain
# ------------------------------------------------------------
ENV MAMBA_ROOT_PREFIX=/opt/conda
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest \
 | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba
RUN micromamba create -y -n tools -c conda-forge -c bioconda \
      python=3.10 \
      orthofinder=2.5.4 \
      diamond=2.1.8 \
      mafft \
      mcl \
      fasttree \
  && micromamba clean -a -y
ENV PATH=$MAMBA_ROOT_PREFIX/envs/tools/bin:$PATH

# ------------------------------------------------------------
# R dependencies (CRAN + Bioconductor)
# ------------------------------------------------------------
RUN R -q -e "options(repos=c(CRAN='https://cloud.r-project.org')); install.packages(c('remotes','BiocManager'))"

# Pin Bioconductor release compatible with R 4.3
RUN R -q -e "BiocManager::install(version = '3.18'); \
  BiocManager::install(c( \
    'AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures', \
    'GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO', \
    'IHW','qvalue','rtracklayer','pwalign' \
  ), update=FALSE, ask=FALSE)"

# CRAN packages used broadly
RUN R -q -e "install.packages(c('cli','data.table','optparse','xml2','ggplot2','dplyr'))"

# Optional but nice for docs/tests
RUN R -q -e "install.packages(c('testthat','knitr','rmarkdown'))"

# RIdeogram (optional; guarded in code)
RUN R -q -e "install.packages('RIdeogram')"

# orthologr + its GitHub deps
RUN R -q -e "install.packages(c('doParallel','foreach','ape','Rdpack','benchmarkme','devtools'))"
RUN R -q -e "remotes::install_github(c('drostlab/metablastr','drostlab/rdiamond','drostlab/orthologr'), upgrade='never')"

# ------------------------------------------------------------
# Build & install dndsR from source with cached deps
# ------------------------------------------------------------
WORKDIR /usr/local/src/dndsR
COPY DESCRIPTION ./
# COPY NAMESPACE ./   # uncomment if present for better caching
RUN R -q -e "remotes::install_deps('.', dependencies=TRUE, upgrade='never')"

COPY . .
RUN R CMD build . && R CMD INSTALL *.tar.gz

# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------
RUN mkdir -p /usr/local/bin
COPY cli/dndsR.R /usr/local/bin/dndsR
RUN chmod +x /usr/local/bin/dndsR

WORKDIR /work
ENTRYPOINT ["dndsR"]
CMD ["--help"]
