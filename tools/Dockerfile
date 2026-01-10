# syntax=docker/dockerfile:1.6
FROM rocker/r2u:jammy

ENV DEBIAN_FRONTEND=noninteractive TZ=UTC \
    LC_ALL=C.UTF-8 LANG=C.UTF-8
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# ============================================================
# Base system libs (graphics + dev headers)
#   - Add pkg-config + librsvg2-dev + libwebp-dev to fix rsvg/ragg/RIdeogram
# ============================================================
RUN set -eux; \
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
# Fail-fast + faster compilation (when it must compile)
# ============================================================
ENV MAKEFLAGS="-j$(nproc)"

# ============================================================
# R dependencies
# ============================================================

# remotes + BiocManager
RUN R -q -e "options(warn=2, Ncpus=parallel::detectCores()); install.packages(c('remotes','BiocManager'))"

# ---- CRAN ----
RUN R -q -e "options(warn=2, Ncpus=parallel::detectCores()); \
  install.packages(c( \
    'cli','data.table','optparse','xml2','ggplot2','dplyr','testthat','knitr','rmarkdown', \
    'RIdeogram','doParallel','foreach','ape','Rdpack','benchmarkme','svglite', \
    'sysfonts','showtext','showtextdb' \
  )); \
  missing <- c('RIdeogram','svglite','sysfonts','showtext')[ \
    !vapply(c('RIdeogram','svglite','sysfonts','showtext'), requireNamespace, logical(1), quietly=TRUE) ]; \
  if (length(missing)) stop('Missing CRAN packages after install: ', paste(missing, collapse=', '))"

# ---- Bioconductor ----
RUN R -q -e "options(warn=2, Ncpus=parallel::detectCores()); \
  BiocManager::install(c( \
    'AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures', \
    'GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO', \
    'IHW','qvalue','rtracklayer','pwalign','txdbmaker' \
  ), update=FALSE, ask=FALSE); \
  missing <- c('AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures','GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO','IHW','qvalue','rtracklayer','pwalign','txdbmaker')[ \
    !vapply(c('AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures','GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO','IHW','qvalue','rtracklayer','pwalign','txdbmaker'), requireNamespace, logical(1), quietly=TRUE) ]; \
  if (length(missing)) stop('Missing Bioc packages after install: ', paste(missing, collapse=', '))"

# ---- GitHub runtime deps ----
# Install metablastr + rdiamond normally, then build orthologr from source with our patch applied.
COPY inst/patches/orthologr-pwalign.patch /tmp/orthologr-pwalign.patch
ARG ORTHOLOGR_REF=HEAD

RUN set -eux; \
  R -q -e "options(warn=2); remotes::install_github(c('drostlab/metablastr','drostlab/rdiamond'), upgrade='never')"; \
  \
  rm -rf /tmp/orthologr; \
  git clone https://github.com/drostlab/orthologr /tmp/orthologr; \
  \
  # allow future pinning a known-good orthologr commit/tag if build breaks
  # build with: docker build --build-arg ORTHOLOGR_REF=<sha-or-tag> ...
  if [ -n "${ORTHOLOGR_REF:-}" ] && [ "${ORTHOLOGR_REF}" != "HEAD" ]; then \
    git -C /tmp/orthologr checkout -q "${ORTHOLOGR_REF}"; \
  fi; \
  \
  # fail early if the patch no longer applies cleanly
  patch --dry-run --forward -p1 -d /tmp/orthologr < /tmp/orthologr-pwalign.patch; \
  patch --forward -p1 -d /tmp/orthologr < /tmp/orthologr-pwalign.patch; \
  R CMD INSTALL /tmp/orthologr; \
  \
  R -q -e "options(warn=2); \
    missing <- c('metablastr','rdiamond','orthologr')[!vapply(c('metablastr','rdiamond','orthologr'), requireNamespace, logical(1), quietly=TRUE)]; \
    if (length(missing)) stop('Missing GitHub packages after install: ', paste(missing, collapse=', '))"

# Ensure patch was applied
RUN R -q -e " \
  stopifnot(requireNamespace('orthologr', quietly=TRUE)); \
  fn <- get('pairwise_aln', asNamespace('orthologr')); \
  txt <- paste(deparse(body(fn)), collapse='\n'); \
  stopifnot(grepl('pwalign::pairwiseAlignment', txt)); \
  stopifnot(grepl('pwalign::pattern', txt)); \
  stopifnot(grepl('pwalign::subject', txt)); \
  stopifnot(grepl('pwalign::writePairwiseAlignments', txt)); \
  cat('orthologr patch OK\n'); \
"

# Final sanity check
RUN Rscript -e "stopifnot( \
  requireNamespace('RIdeogram', quietly=TRUE), \
  requireNamespace('topGO', quietly=TRUE), \
  requireNamespace('GO.db', quietly=TRUE), \
  requireNamespace('GenomicFeatures', quietly=TRUE) \
)"

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

WORKDIR /work
CMD ["bash"]

LABEL \
  org.opencontainers.image.source="https://github.com/scrumpis/dndsr" \
  org.opencontainers.image.description="dndsR runtime dependencies (R/Bioc + OrthoFinder toolchain)"
