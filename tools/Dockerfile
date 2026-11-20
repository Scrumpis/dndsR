# syntax=docker/dockerfile:1.6
# Pinned by digest. To update, change the digest here and rerun the validation in CI.
FROM rocker/r2u@sha256:b51f03901dd738476588240f0c37f359062c755b509213d4273b994e4923faea

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
# R dependencies
# ============================================================

# remotes + BiocManager
RUN R -q -e "options(warn=2, Ncpus=parallel::detectCores()); install.packages(c('remotes','BiocManager'))"

# ---- CRAN ----
RUN R -q -e "options(warn=2, Ncpus=parallel::detectCores()); \
  install.packages(c( \
    'cli','data.table','lme4','xml2','ggplot2','dplyr','dtplyr','R.utils','testthat','withr', \
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

# ---- GitHub runtime deps, pinned by commit ----
# orthologr has included the pwalign fix upstream since 2026, so no patch is needed.
ARG ORTHOLOGR_REF=06c3098adfba1e32135d904f72ea81a76902c664
ARG METABLASTR_REF=02dad6daf9f4a0675317177ff434ae9bbcec535f
ARG RDIAMOND_REF=d8b4006d5f2793ecf8b97a7560fa005b7a3420ee
# Source archives of those commits from codeload.github.com, checked by SHA-256, then their
# CRAN/Bioconductor dependencies and R CMD INSTALL in dependency order. This uses neither the
# GitHub API (rate-limited without a token) nor any credential.
ARG METABLASTR_SHA256=61a4c3eac5ec6218585fba068d336a140681a2e75b799e132d02949c6b1f2303
ARG RDIAMOND_SHA256=1172107720a0935163a0429515dbc76d455c26bda1de684c4c39ff61fd8a763f
ARG ORTHOLOGR_SHA256=037a53d2a0da203f459c333f5215d0eb4cde699ccd88585a744fa523065fe6fa
RUN set -eux; \
  d=$(mktemp -d); cd "$d"; \
  for spec in "metablastr ${METABLASTR_REF} ${METABLASTR_SHA256}" "rdiamond ${RDIAMOND_REF} ${RDIAMOND_SHA256}" "orthologr ${ORTHOLOGR_REF} ${ORTHOLOGR_SHA256}"; do \
    set -- $spec; \
    curl -fsSL --retry 5 --retry-all-errors -o "$1.tar.gz" "https://codeload.github.com/drostlab/$1/tar.gz/$2"; \
    echo "$3  $1.tar.gz" | sha256sum -c -; \
    tar -xzf "$1.tar.gz"; \
    R -q -e "options(warn = 2, repos = BiocManager::repositories()); \
      deps <- setdiff(remotes::local_package_deps('$1-$2', dependencies = NA), c(rownames(installed.packages()), 'R', 'metablastr', 'rdiamond')); \
      if (length(deps)) install.packages(deps)"; \
    R CMD INSTALL "$1-$2"; \
  done; \
  cd /; rm -rf "$d"; \
  R -q -e "missing <- c('metablastr','rdiamond','orthologr')[!vapply(c('metablastr','rdiamond','orthologr'), requireNamespace, logical(1), quietly=TRUE)]; \
    if (length(missing)) stop('Missing GitHub packages after install: ', paste(missing, collapse=', '))"

# Final sanity check
RUN Rscript -e "stopifnot( \
  requireNamespace('RIdeogram', quietly=TRUE), \
  requireNamespace('topGO', quietly=TRUE), \
  requireNamespace('GO.db', quietly=TRUE), \
  requireNamespace('GenomicFeatures', quietly=TRUE) \
)"

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
# Micromamba environment (python, OrthoFinder toolchain, HyPhy)
# ============================================================
ENV MAMBA_ROOT_PREFIX=/opt/conda
ARG MAMBA_VER=2.5.0
RUN set -eux; \
  curl -fsSL "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VER}" -o /tmp/micromamba.tar.bz2; \
  tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/bin --strip-components=1 bin/micromamba; \
  rm -f /tmp/micromamba.tar.bz2

# HyPhy is the engine for phylo_busted, phylo_absrel, and phylo_relax; PAML's codeml for
# phylo_branch_omega. Both are pinned to the releases dndsR was validated against:
# docs/validation/phylo-busted-wrapper.md, docs/validation/phylo-lineage.md. KaKs_Calculator
# computes the dN/dS methods other than Comeron: docs/validation/estimators.md.
RUN micromamba create -y -n tools -c conda-forge -c bioconda \
      python=3.10.19 \
      orthofinder=2.5.4 \
      mafft=7.526 mcl=22.282 fasttree=2.2.0 \
      hyphy=2.5.103 \
      paml=4.10.10 \
      kakscalculator2=2.0.1 \
  && micromamba clean -a -y \
  && $MAMBA_ROOT_PREFIX/envs/tools/bin/hyphy --version \
  && test -x $MAMBA_ROOT_PREFIX/envs/tools/bin/codeml \
  && $MAMBA_ROOT_PREFIX/envs/tools/bin/KaKs_Calculator -h | grep -q "Version: 2.0"

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
# dndsR itself, so the image runs without any other installation
# ============================================================
COPY . /opt/dndsR
RUN R CMD INSTALL --clean /opt/dndsR \
  && install -m 0755 /opt/dndsR/inst/scripts/dndsr /usr/local/bin/dndsr \
  && Rscript -e "stopifnot(requireNamespace('dndsR', quietly=TRUE))" \
  && dndsr __commands \
  && Rscript -e "ip <- installed.packages()[, c('Package', 'Version')]; write.table(ip[order(ip[, 'Package']), ], '/opt/dndsr-package-versions.tsv', sep = '\t', quote = FALSE, row.names = FALSE); cat(R.version.string, file = '/opt/dndsr-r-version.txt')"

WORKDIR /work
CMD ["dndsr", "--help"]

LABEL \
  org.opencontainers.image.source="https://github.com/scrumpis/dndsr" \
  org.opencontainers.image.description="dndsR and its runtime dependencies"
