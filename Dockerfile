# ------------------------------------------------------------
# Base: R 4.3.3 (Debian Bookworm) + system build deps
# ------------------------------------------------------------
FROM rocker/r-ver:4.3.3 AS rbase

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC

# System libraries needed by CRAN/BioC packages (xml2, curl, SSL, compilers, etc.)
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
# Add micromamba (tiny conda) to manage Orthofinder/DIAMOND, etc.
# ------------------------------------------------------------
ENV MAMBA_ROOT_PREFIX=/opt/conda
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin bin/micromamba --strip-components=1

# Create a clean env with precise versions (Orthofinder 2.5.4, DIAMOND 2.1.8)
# and common dependencies Orthofinder likes to find on PATH.
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
# R package installation (CRAN + Bioconductor)
# ------------------------------------------------------------
# Pre-seed repos for fast installs
RUN R -q -e "options(repos=c(CRAN='https://cloud.r-project.org')); install.packages(c('remotes','BiocManager'))"

# Install CRAN packages (Imports + runtime)
# Note: 'stats' and 'utils' are base; 'RIdeogram' spelling matches CRAN.
RUN R -q -e "install.packages(c( \
  'cli','data.table','optparse','readr','xml2','ggplot2','dplyr' \
))"

# Install Bioconductor packages (Imports + optional)
# GO.db is large; pin BiocManager to current Bioc release for R-4.3
RUN R -q -e "BiocManager::install(c( \
  'AnnotationDbi','BiocGenerics','Biostrings','GenomicFeatures', \
  'GenomeInfoDb','GO.db','IRanges','Rsamtools','S4Vectors','topGO', \
  'IHW','qvalue' \
), update=FALSE, ask=FALSE)"

# RIdeogram (CRAN; separate to surface errors clearly)
RUN R -q -e "install.packages('RIdeogram')"

# Suggested packages (nice to have; harmless on HPC)
RUN R -q -e "install.packages(c('testthat','knitr','rmarkdown'))"

# CRAN helpers
RUN R -q -e "options(repos=c(CRAN='https://cloud.r-project.org')); install.packages(c('remotes','BiocManager'))"

# Bioconductor dependencies used by orthologr chain
RUN R -q -e "BiocManager::install(c( \
  'Biostrings','GenomicRanges','GenomicFeatures','Rsamtools','rtracklayer','pwalign' \
), update=FALSE, ask=FALSE)"

# CRAN deps used by orthologr
RUN R -q -e "install.packages(c('doParallel','foreach','ape','Rdpack','benchmarkme','devtools'))"

# GitHub deps for orthologr
RUN R -q -e "remotes::install_github(c('drostlab/metablastr','drostlab/rdiamond','drostlab/orthologr'), upgrade='never')"

# ------------------------------------------------------------
# Build and install your package from source
# (Copy *only* DESCRIPTION/NAMESPACE first for better layer caching)
# ------------------------------------------------------------
WORKDIR /usr/local/src/dndsR
COPY DESCRIPTION ./
# if you have NAMESPACE at the root, copy it too for caching:
# COPY NAMESPACE ./

# Install (and cache) just the dependencies computed from DESCRIPTION
RUN R -q -e "remotes::install_deps('.', dependencies=TRUE, upgrade='never')"

# Now bring the rest of the source and install the package
COPY . .
RUN R CMD build . && \
    R CMD INSTALL *.tar.gz

# ------------------------------------------------------------
# Provide a friendly CLI; assume you keep the script at tools/dndsR.R
# If your script lives elsewhere, change the COPY path accordingly.
# ------------------------------------------------------------
# Prefer installing a tiny wrapper that loads package and calls the script.
# If your repo already has tools/dndsR.R as the exact CLI we finalized, copy it in:
RUN mkdir -p /usr/local/bin
COPY cli/dndsR.R /usr/local/bin/dndsR
RUN chmod +x /usr/local/bin/dndsR

# ------------------------------------------------------------
# Defaults
# ------------------------------------------------------------
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8
WORKDIR /work
ENTRYPOINT ["dndsR"]
CMD ["--help"]
