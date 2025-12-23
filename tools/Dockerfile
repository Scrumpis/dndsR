# syntax=docker/dockerfile:1.6
FROM rocker/r2u:jammy

ENV DEBIAN_FRONTEND=noninteractive TZ=UTC \
    LC_ALL=C.UTF-8 LANG=C.UTF-8
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# ============================================================
# Base system libs (graphics + dev headers)
# ============================================================
RUN set -eux; \
  apt-get update -qq; \
  apt-get install -y --no-install-recommends \
    ca-certificates gnupg2 curl wget git make g++ \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    libpcre2-dev libreadline-dev \
    libpng-dev libjpeg-dev libtiff5-dev \
    libcairo2-dev libxt-dev libfontconfig1-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev bzip2 \
  && rm -rf /var/lib/apt/lists/*

# ============================================================
# Disable bspm auto-enable (breaks under Singularity / no apt access)
# ============================================================
RUN set -eux; \
  for f in /etc/R/Rprofile.site /usr/lib/R/etc/Rprofile.site; do \
    if [ -f "$f" ]; then \
      sed -i \
        -e '/bspm::enable/d' \
        -e '/options(bspm\./d' \
        "$f"; \
      echo 'options(bspm.enabled = FALSE)' >> "$f"; \
    fi; \
  done
ENV BSPM_DISABLE=1

# ============================================================
# Micromamba environment (python + OrthoFinder toolchain)
# ============================================================
ENV MAMBA_ROOT_PREFIX=/opt/conda
ARG MAMBA_VER=latest
RUN set -eux; \
  curl -fsSL "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VER}" -o /tmp/micromamba.tar.bz2; \
  tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/bin --strip-components=1 bin/micromamba; \
  rm -f /tmp/micromamba.tar.bz2

RUN micromamba create -y -n tools -c conda-forge -c bioconda \
      python=3.10 \
      orthofinder=2.5.4 \
      mafft mcl fasttree \
  && micromamba clean -a -y

# ============================================================
# DIAMOND binary
# ============================================================
ARG DIAMOND_VER=2.1.14
RUN set -eux; \
  curl -fsSL -o /tmp/diamond.tgz \
    "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VER}/diamond-linux64.tar.gz"; \
  tar -xvzf /tmp/diamond.tgz -C /usr/local/bin diamond; \
  chmod 0755 /usr/local/bin/diamond; \
  rm -f /tmp/diamond.tgz; \
  /usr/local/bin/diamond version

ENV PATH=/usr/local/bin:$MAMBA_ROOT_PREFIX/envs/tools/bin:$PATH

# ============================================================
# R/Bioconductor dependencies
# ============================================================
RUN R -q -e "install.packages(c('remotes','BiocManager'), Ncpus=2)"

RUN R -q -e "BiocManager::install(c( \
    'AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures', \
    'GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO', \
    'IHW','qvalue','rtracklayer','pwalign','txdbmaker' \
  ), update=FALSE, ask=FALSE, Ncpus=2)"

RUN R -q -e "install.packages(c( \
    'cli','data.table','optparse','xml2','ggplot2','dplyr','testthat','knitr','rmarkdown', \
    'RIdeogram','doParallel','foreach','ape','Rdpack','benchmarkme','devtools','svglite', \
    'sysfonts','showtext','showtextdb' \
  ), Ncpus=2)"

RUN R -q -e "remotes::install_github(c( \
    'drostlab/metablastr','drostlab/rdiamond','drostlab/orthologr' \
  ), upgrade='never')"

WORKDIR /work
CMD ["bash"]

LABEL \
  org.opencontainers.image.source="https://github.com/scrumpis/dndsr" \
  org.opencontainers.image.description="dndsR runtime (R 4.5 + Bioconductor + OrthoFinder toolchain)"
