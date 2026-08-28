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

genR_SE ## Name your SE object as genR_SE.

## Remove low expressed genes
keep <- filterByExpr(assay(genR_SE), group = genR_SE$sex, min.count = 10, min.prop = 0.7)
genR_SE_filt <- genR_SE[keep, ]
#' This block of code removes low expressed genes. 

## Run SVA

## Subtitute blood_cell_proportions with the actual names of the columns in your colData that contain the blood cell proportions.
mod1_gtex <- model.matrix(~sex + age + blood_cell_proportions, colData(genR_SE_filt) ) 
mod0_gtex <- model.matrix(~ age + blood_cell_proportions, colData(genR_SE_filt))
svseq <- svaseq(assay(genR_SE_filt), mod1_gtex, mod0_gtex)$sv
colnames(svseq) <- paste0("SV", seq_len(ncol(svseq))) ## Rename SV columns
rownames(svseq) <- colnames(genR_SE_filt) ## Set sample names as SV rownames
colData(genR_SE_filt) <- cbind(colData(genR_SE_filt), svseq) ## Add SVs to colData
#' This code compute surrogate variables (SV) for RNAseq data. SV are then added 
#' to the SE as additional colData variables. 
save(genR_SE_filt, file = "results/preprocessFiles/gexp_genR_SE_filt_sva_sex.Rdata")
