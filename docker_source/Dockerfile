# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_23

# Install R packages
RUN R -e 'BiocManager::install(c("tidyverse", "ggExtra", "cowplot", "S4Vectors", "readxl", "ggcorrplot", "forestploter", "ggmanh", "ggrepel", "ggbreak", "ggpubr", "poolr"))'

# Install Bioc packages
RUN R -e 'BiocManager::install(c("Biobase", "SummarizedExperiment", "sva", "isva", "SmartSVA", "limma", "UpSetR"))'

# Install Bioc packages (1.1)
RUN R -e 'BiocManager::install(c("enrichR"))'

# Install Bioc packages (1.2)
RUN R -e 'BiocManager::install(c("GEOquery", "affy", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))'

# Install Bioc packages (1.3)
RUN R -e 'BiocManager::install(c("meta", "eulerr"))'

# Install Bioc packages (1.4)
RUN R -e 'BiocManager::install(c("recount3"))'


# Install compositions packages (1.5)
RUN R -e 'BiocManager::install(c("compositions"))'

