#'##############################################################################
#' Template for running the dimorphic gene analysis in GenerationR
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(limma)
library(SummarizedExperiment)

## Load data ####
load("results/preprocessFiles/gexp_genR_SE_filt_sva_sex.Rdata") ## Output of preprocess_GenerationR.R

genR_SE ## This should be the name of the SE object in the Rdata file.
#' This block of code removes genes in the Y chromosome
genR_SE_noY <- genR_SE[as.character(seqnames(genR_SE)) %in% paste0("chr", c(1:22, "X")), ]


## Normalize count data
mod_genR <- model.matrix(formula(paste("~ sex + age + blood_cell_proportions +", , 
                                       paste(paste0("SV", 1:37), collapse = "+"))),
                         data = colData(genR_SE_noY)) 
genR_voom <- voom(assay(genR_SE_noY), mod_genR)
save(genR_voom, file = "paper_sex/results/genR_voom.Rdata")
#' Modify mod_genR to include the number of SVs found in the data. Cell type proportions
#' should also be added to mod_genR model.

## Run analysis
### Main model
lm.genR_main <- lmFit(genR_voom, mod_genR) %>%
  eBayes()
tab.genR_main <- topTable(lm.genR_main, coef = 2, n = Inf, confint = TRUE)
tab.genR_main$Symbol <- rowData(genR_SE_noY)[rownames(tab.genR_main), "gene_name"]
#' GTEx genes were already mapped to gene Symbol. This line of code retrieves the
#' information and assign it to a column in the results table.
tab.genR_main$chromosome <- as.character(seqnames(genR_SE_noY[rownames(tab.genR_main)]))


### Scaled model
genR_voom$E <- t(scale(t(genR_voom$E)))
#' To ensure comparability of effect sizes among the cohorts, the data will be
#' centered and scaled. Be careful because this step will overwrite the 
#' voom-transformed RNAseq data.
lm.genR_scaled <- lmFit(genR_voom, mod_genR) %>%
  eBayes()
tab.genR_scaled <- topTable(lm.genR_scaled, coef = 2, n = Inf, confint = TRUE)
tab.genR_scaled$Symbol <- rowData(genR_SE_noY)[rownames(tab.genR_scaled), "gene_name"]
tab.genR_scaled$chromosome <- as.character(seqnames(genR_SE_noY[rownames(tab.genR_scaled)]))
save(tab.genR_scaled, tab.genR_main, file = "paper_sex/results/genR_analysis.Rdata")
#' This file should be shared with us for the manuscript preparation.