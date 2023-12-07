FROM rocker/r-ver:4.3
LABEL maintainer=" Brent T. Schlegel [bts76@pitt.edu]" \
    description="Environment and dependencies for BCR Nextflow pipeline"

# Install ps, for Nextflow. https://www.nextflow.io/docs/latest/tracing.html
RUN apt-get update && \
    apt-get install -y procps \
    libglpk40 \
    pandoc \
    libcurl4-openssl-dev \
    r-cran-curl

# Install required R packages
ARG R_DEPS="c( \
    'BiocManager', \
    'dplyr', \
    'ggplot2', \
    'viridisLite', \
    'colorRamp2', \
    'Platypus', \
    'optparse', \
    'stringr', \
    'tidyr', \
    'devtools', \
    'Seurat', \
    'utils', \
    'scales', \
    'ade4', \
    'phylobase' \
)"
ARG R_BIOC_DEPS="c( \
    'Biostrings', \
    'msa', \
    'GenomicRanges', \
    'GenomicAlignments',\
    'BrepPhylo', \
    'org.Mm.eg.db', \
    'org.Hs.eg.db', \
    'edgeR', \
    'fgsea', \
    'ggtree' \
    )"

RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE, dependencies=TRUE)"
RUN Rscript -e "BiocManager::install(${R_BIOC_DEPS}, dependencies=TRUE)"
RUN Rscript -e "install.packages('NMF', clean=TRUE)"
RUN Rscript -e "devtools::install_github('fmichonneau/phylobase')"
RUN Rscript -e "devtools::install_github('Fraternalilab/BrepPhylo',dependencies = TRUE)"
RUN Rscript -e "devtools::install_github('jokergoo/ComplexHeatmap')"
CMD ["R"]
