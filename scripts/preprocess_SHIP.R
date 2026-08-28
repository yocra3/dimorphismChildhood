#'##############################################################################
#' Template for preprocessing data for GenerationR
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(SummarizedExperiment)
library(sva)

## Load data ####
load("Gexp.Rdata")
#' Gene expression data should be stored in a SummarizedExperiment (SE). The SE should 
#' contain log2 expression values (as normalized by the cohort) and collapsed to a gene SYMBOL.
#' The colData should contain the covariates of interest, including sex, age and cell type proportions.
#' In case the data is not in a SE object, you can check the SummarizedExperiment 
#' package documentation: https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html.
#' Genes with low call rate or bad quality should have been already removed.
ship_SE ## Name your SE object as ship_SE.

## Run SVA

mat <- assay(ship_SE)
pd <- colData(ship_SE)

## Subtitute blood_cell_proportions and technical_variables with the actual 
## names of the columns in your colData that contain the blood cell proportions.
mod1_ship <- model.matrix(~sex + age + blood_cell_proportions + technical_variables, pd ) 
mod0_ship <- model.matrix(~ age + blood_cell_proportions + technical_variables, pd)

n.sv <- num.sv(mat, mod1_ship, method="leek") + 1


sv.obj <- sva(mat, mod1_ship, mod0 = mod0_ship, n.sv = n.sv)

## Add SVs to SE
ship_sv <- sv.obj$sv
colnames(ship_sv) <- paste0("SV", seq_len(ncol(ship_sv)))
rownames(ship_sv) <- colnames(ship_SE)

colData(ship_SE) <- cbind(colData(ship_SE), ship_sv)

# Save
save(ship_SE, file = "results/preprocessFiles/SHIP_SE_SV.RData")
