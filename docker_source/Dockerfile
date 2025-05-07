# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install R packages
RUN R -e 'BiocManager::install(c("tidyverse", "ggExtra", "cowplot", "S4Vectors", "readxl", "ggcorrplot", "forestploter", "ggmanh", "ggrepel", "ggbreak", "mediation"))'

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

# Install Bioc packages (1.5)
RUN R -e 'BiocManager::install(c("poolr"))'

# Install MR (1.6)
RUN R -e 'remotes::install_github("MRCIEU/TwoSampleMR")'
RUN R -e 'remotes::install_github("MRCIEU/MRInstruments")'

# Install ggpubr (1.7)
RUN R -e 'BiocManager::install(c("ggpubr"))'

# Update enrichR version (1.8)
RUN R -e 'BiocManager::install(c("enrichR"))'
