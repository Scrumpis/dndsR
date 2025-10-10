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
# Micromamba (multi-arch safe)
# ------------------------------------------------------------
ENV MAMBA_ROOT_PREFIX=/opt/conda
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends bzip2 && rm -rf /var/lib/apt/lists/*

ARG TARGETPLATFORM
ARG MAMBA_VER=latest
RUN set -eux; \
  case "${TARGETPLATFORM}" in \
    linux/amd64)  MAMBA_ARCH=linux-64 ;; \
    linux/arm64)  MAMBA_ARCH=linux-aarch64 ;; \
    *) echo "Unsupported platform: ${TARGETPLATFORM}" >&2; exit 1 ;; \
  esac; \
  curl -fsSL "https://micro.mamba.pm/api/micromamba/${MAMBA_ARCH}/${MAMBA_VER}" -o /tmp/micromamba.tar.bz2; \
  tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/bin --strip-components=1 bin/micromamba; \
  rm -f /tmp/micromamba.tar.bz2

# ------------------------------------------------------------
# Tools env (NO diamond here to avoid old 0.9 pulls)
# ------------------------------------------------------------
RUN micromamba create -y -n tools -c conda-forge -c bioconda \
      python=3.10 \
      orthofinder=2.5.4 \
      mafft \
      mcl \
      fasttree \
  && micromamba clean -a -y

# ------------------------------------------------------------
# DIAMOND v2.x per-arch:
#  - amd64: official GitHub binary (flat tarball -> single 'diamond' file)
#  - arm64: Bioconda (>=2.1.12) since GitHub doesn’t publish arm64 tarballs
# ------------------------------------------------------------
ARG DIAMOND_VER=2.1.14
RUN set -eux; \
  case "${TARGETPLATFORM}" in \
    linux/amd64) \
      curl -fsSL -o /tmp/diamond.tgz \
        "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VER}/diamond-linux64.tar.gz"; \
      # tarball contains a single file named 'diamond' at top level
      tar -xvzf /tmp/diamond.tgz -C /usr/local/bin diamond; \
      chmod 0755 /usr/local/bin/diamond; \
      rm -f /tmp/diamond.tgz; \
      /usr/local/bin/diamond version \
    ;; \
    linux/arm64) \
      micromamba install -y -n tools -c bioconda -c conda-forge "diamond>=2.1.12"; \
      micromamba clean -a -y; \
      micromamba run -n tools diamond version \
    ;; \
    *) echo "Unsupported platform for DIAMOND: ${TARGETPLATFORM}" >&2; exit 1 ;; \
  esac

ENV PATH=/usr/local/bin:$MAMBA_ROOT_PREFIX/envs/tools/bin:$PATH

# ------------------------------------------------------------
# R dependencies (CRAN + Bioconductor)
# ------------------------------------------------------------
RUN R -q -e "options(repos=c(CRAN='https://cloud.r-project.org')); install.packages(c('remotes','BiocManager'))"

RUN R -q -e "BiocManager::install(version = '3.18'); \
  BiocManager::install(c( \
    'AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures', \
    'GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO', \
    'IHW','qvalue','rtracklayer','pwalign' \
  ), update=FALSE, ask=FALSE)"

RUN R -q -e "install.packages(c('cli','data.table','optparse','xml2','ggplot2','dplyr'))"
RUN R -q -e "install.packages(c('testthat','knitr','rmarkdown'))"
RUN R -q -e "install.packages('RIdeogram')"

# orthologr + deps
RUN R -q -e "install.packages(c('doParallel','foreach','ape','Rdpack','benchmarkme','devtools'))"
RUN R -q -e "remotes::install_github(c('drostlab/metablastr','drostlab/rdiamond','drostlab/orthologr'), upgrade='never')"

# ------------------------------------------------------------
# Build & install dndsR
# ------------------------------------------------------------
WORKDIR /usr/local/src/dndsR
COPY DESCRIPTION ./
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
