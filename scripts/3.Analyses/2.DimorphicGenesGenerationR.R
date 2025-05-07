#'##############################################################################
#' Template for running the dimorphic gene analysis in GenerationR
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(limma)
library(SummarizedExperiment)

## Load data ####
load("results/preprocessFiles/gexp_gtex_sva_sex.Rdata") ## Output of preprocess_GenerationR.R

## Remove bad genes
pseudo_reg <- GRanges(c("chrX:60001-2699520", "chrX:154931044-155260560")) 
## Coordinates of pseudoautosomal region in hg19. Modify it if using GRCh38 genome build
pseudo_genes_gtex <- subsetByOverlaps(rowRanges(rse_gtex_filt), pseudo_reg)
rse_gtex_noY <- rse_gtex_filt[!rownames(rse_gtex_filt) %in% names(pseudo_genes_gtex) & 
                                as.character(seqnames(rse_gtex_filt)) %in% paste0("chr", c(1:22, "X")), ]
#' This block of code removes genes in the pseudoautosomal region


## Normalize count data
mod_gtex <- model.matrix(formula(paste("~ sex + gtex.age +", 
                                       paste(paste0("SV", 1:37), collapse = "+"))),
                         data = colData(rse_gtex_noY)) 
gtex_voom <- voom(assay(rse_gtex_noY), mod_gtex)
save(gtex_voom, file = "paper_sex/results/gtex_voom.Rdata")

## Run analysis
gtex_voom$E <- t(scale(t(gtex_voom$E)))
#' To ensure comparability of effect sizes among the cohorts, the data will be
#' centered and scaled. Be careful because this step will overwrite the 
#' voom-transformed RNAseq data.
lm.gtex <- lmFit(gtex_voom, mod_gtex) %>%
  eBayes()
tab.gtex <- topTable(lm.gtex, coef = 2, n = Inf, confint = TRUE)
tab.gtex$Symbol <- rowData(rse_gtex_noY)[rownames(tab.gtex), "gene_name"]
#' GTEx genes were already mapped to gene Symbol. This line of code retrieves the
#' information and assign it to a column in the results table.
tab.gtex$chromosome <- as.character(seqnames(rse_gtex_noY[rownames(tab.gtex)]))
save(tab.gtex, file = "paper_sex/results/gtex_analysis.Rdata")
