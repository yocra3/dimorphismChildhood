#'##############################################################################
#' Template for preprocessing data for GenerationR
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(edgeR)
library(limma)
library(SummarizedExperiment)
library(sva)

## Load data ####
load("Gexp.Rdata")
#' RNAseq data should be stored in a SummarizedExperiment (SE). The SE should 
#' contain raw counts and a column in the colData called sex with two levels: male and female.
#' In case the data is not in a SE object, you can check the SummarizedExperiment 
#' package documentation: https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html.

## Remove low expressed genes
keep <- filterByExpr(assay(rse_gtex), group = rse_gtex$sex, min.count = 10, min.prop = 0.7)
rse_gtex_filt <- rse_gtex[keep, ]
#' This block of code removes low expressed genes. 

## Run SVA
mod1_gtex <- model.matrix(~sex + gtex.age, colData(rse_gtex_filt) )
mod0_gtex <- model.matrix(~ gtex.age, colData(rse_gtex_filt))
svseq <- svaseq(assay(rse_gtex_filt), mod1_gtex, mod0_gtex)$sv
colnames(svseq) <- paste0("SV", seq_len(ncol(svseq))) ## Rename SV columns
rownames(svseq) <- colnames(rse_gtex_filt) ## Set sample names as SV rownames
colData(rse_gtex_filt) <- cbind(colData(rse_gtex_filt), svseq) ## Add SVs to colData
#' This code compute surrogate variables (SV) for RNAseq data. SV are then added 
#' to the SE as additional colData variables. 
save(rse_gtex_filt, file = "results/preprocessFiles/gexp_gtex_sva_sex.Rdata")
