# syntax=docker/dockerfile:1.6
# ------------------------------------------------------------
# Base: R 4.3.3 (Debian Bookworm)
# ------------------------------------------------------------
FROM rocker/r-ver:4.3.3

ENV DEBIAN_FRONTEND=noninteractive TZ=UTC \
    LC_ALL=C.UTF-8 LANG=C.UTF-8

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Base syslibs (incl. bzip2 for micromamba tarball)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates gnupg curl wget git make g++ \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    libpcre2-dev libreadline-dev \
    libpng-dev libjpeg-dev libtiff5-dev \
    libcairo2-dev libxt-dev libfontconfig1-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev bzip2 && \
    rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# Micromamba (linux-64 only) + tools env
# ------------------------------------------------------------
ENV MAMBA_ROOT_PREFIX=/opt/conda
ARG MAMBA_VER=latest

RUN set -eux; \
  curl -fsSL "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VER}" -o /tmp/micromamba.tar.bz2; \
  tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/bin --strip-components=1 bin/micromamba; \
  rm -f /tmp/micromamba.tar.bz2

RUN micromamba create -y -n tools -c conda-forge -c bioconda \
      python=3.10 \
      orthofinder=2.5.4 \
      mafft \
      mcl \
      fasttree \
  && micromamba clean -a -y

# DIAMOND v2.x (official linux64 tarball)
ARG DIAMOND_VER=2.1.14
RUN set -eux; \
  curl -fsSL -o /tmp/diamond.tgz \
    "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VER}/diamond-linux64.tar.gz"; \
  tar -xvzf /tmp/diamond.tgz -C /usr/local/bin diamond; \
  chmod 0755 /usr/local/bin/diamond; \
  rm -f /tmp/diamond.tgz; \
  /usr/local/bin/diamond version

ENV PATH=/usr/local/bin:$MAMBA_ROOT_PREFIX/envs/tools/bin:$PATH

# ------------------------------------------------------------
# R dependencies (CRAN + Bioconductor + orthologr deps)
# (No dndsR source is copied into the image)
# ------------------------------------------------------------
RUN R -q -e "options(repos=c(CRAN='https://cloud.r-project.org')); install.packages(c('remotes','BiocManager'), Ncpus=2)"

RUN R -q -e "BiocManager::install(version = '3.18'); \
  BiocManager::install(c( \
    'AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures', \
    'GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO', \
    'IHW','qvalue','rtracklayer','pwalign' \
  ), update=FALSE, ask=FALSE)"

RUN R -q -e "install.packages(c('cli','data.table','optparse','xml2','ggplot2','dplyr','testthat','knitr','rmarkdown','RIdeogram','doParallel','foreach','ape','Rdpack','benchmarkme','devtools'), Ncpus=2)"

RUN R -q -e "remotes::install_github(c('drostlab/metablastr','drostlab/rdiamond','drostlab/orthologr'), upgrade='never')"

# ------------------------------------------------------------
# Runtime bootstrapper: install dndsR from GitHub or local mount, then exec CLI
# - DNDSR_REF   : tag/branch/SHA (default: main)
# - DNDSR_LOCAL : path to local package dir to install (e.g., /work)
# ------------------------------------------------------------
RUN cat >/usr/local/bin/bootstrap_dndsr.sh <<'EOF' && chmod +x /usr/local/bin/bootstrap_dndsr.sh
#!/usr/bin/env bash
set -euo pipefail

REF="${DNDSR_REF:-main}"
LOCAL="${DNDSR_LOCAL:-}"

# Prefer local source if provided
if [[ -n "${LOCAL}" && -d "${LOCAL}" ]]; then
  echo "[dndsr] Installing local package from ${LOCAL} ..."
  R -q -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes'); remotes::install_local('${LOCAL}', upgrade='never')"
else
  echo "[dndsr] Installing from GitHub scrumpis/dndsr@${REF} ..."
  R -q -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes'); remotes::install_github('scrumpis/dndsr', ref='${REF}', upgrade='never')"
fi

# Try to find a CLI script shipped with the package (e.g., inst/cli/dndsR.R)
CLI_PATH=$(R -q -e "cat(system.file('cli','dndsR.R', package='dndsR'))" 2>/dev/null || true)

if [[ -n "${CLI_PATH}" && -f "${CLI_PATH}" ]]; then
  echo "[dndsr] Launching CLI: ${CLI_PATH}"
  exec Rscript "${CLI_PATH}" "$@"
else
  echo "[dndsr] CLI script not found in package; trying to call dndsR::main()"
  exec Rscript -e 'if (!"dndsR" %in% rownames(installed.packages())) stop("dndsR not installed"); if (!exists("main", asNamespace("dndsR"))) stop("No dndsR::main() found"); dndsR::main()' "$@"
fi
EOF

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/bootstrap_dndsr.sh"]
CMD ["--help"]

# Optional metadata
LABEL org.opencontainers.image.source="https://github.com/scrumpis/dndsr" \
      org.opencontainers.image.description="dndsR runtime (deps only) that installs the package at run time from GitHub or a local mount"
