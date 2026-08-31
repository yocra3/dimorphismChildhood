#'##############################################################################
#' Template for running the dimorphic gene analysis in GenerationR
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(limma)
library(SummarizedExperiment)

## Load data ####
load("results/preprocessFiles/SHIP_SE_SV.RData") ## Output of preprocess_GenerationR.R

ship_SE ## This should be the name of the SE object in the Rdata file.
#' This block of code removes genes in the Y chromosome
ship_SE_noY <- ship_SE[as.character(seqnames(ship_SE)) %in% paste0("chr", c(1:22, "X")), ]


## Run analysis
### Main model
## Substitute blood_cell_proportions and technical_variables with the actual names of the columns in your colData 
## Set the number of SVs according to the computation in preprocess_SHIP.R
mod_ship <- model.matrix(formula(paste("~ sex + age + blood_cell_proportions + technical_variables", 
                                              paste(paste0("SV", 1:3), collapse = "+"))),
                            colData(ship_SE_noY ))

lm.main <- lmFit(assay(ship_SE_noY), mod_ship) %>%
  eBayes()
tab.ship.main <- topTable(lm.main, coef = 2, n = Inf, confint = TRUE)
tab.ship.main$chr <- as.character(seqnames(ship_SE_noY[rownames(tab.ship.main), ]))
tab.ship.main$Symbol <-  rowData(ship_SE_noY)[rownames(tab.ship.main), ]$GeneSymbolDB2


### Scaled model
lm.scaled <- lmFit(t(scale(t(assay(ship_SE_noY)))), mod_ship) %>%
  eBayes()
tab.ship.scaled <- topTable(lm.scaled, coef = 2, n = Inf, confint = TRUE)
tab.ship.scaled$chr <- as.character(seqnames(ship_SE_noY[rownames(tab.ship.scaled), ]))
tab.ship.scaled$Symbol <-  rowData(ship_SE_noY)[rownames(tab.ship.scaled), ]$GeneSymbolDB2

## Stratified models
### Pre-menopause
## Filter your SE to have only pre-menopause samples as reported in the
## analysis plan. Remember to select men individuals inside women's age range.
ship_SE_pre

## This model is the same than in the main model but including oral contraceptive.
mod_ship_pre <- model.matrix(formula(paste("~ sex + age + blood_cell_proportions + technical_variables + oral_contraceptive", 
                                              paste(paste0("SV", 1:3), collapse = "+"))),
                            colData(ship_SE_pre ))
lm.pre_main <- lmFit(assay(ship_SE_pre), mod_ship_pre) %>%
  eBayes()
tab.ship.pre_main <- topTable(lm.pre_main, coef = 2, n = Inf, confint = TRUE)
tab.ship.pre_main$chr <- as.character(seqnames(ship_SE_pre[rownames(tab.ship.pre_main), ]))
tab.ship.pre_main$Symbol <-  rowData(ship_SE_pre)[rownames(tab.ship.pre_main), ]$GeneSymbolDB2

lm.pre_scaled <- lmFit(t(scale(t(assay(ship_SE_pre)))), mod_ship_pre) %>%
  eBayes()
tab.ship.pre_scaled <- topTable(lm.pre_scaled, coef = 2, n = Inf, confint = TRUE)
tab.ship.pre_scaled$chr <- as.character(seqnames(ship_SE_pre[rownames(tab.ship.pre_scaled), ]))
tab.ship.pre_scaled$Symbol <-  rowData(ship_SE_pre)[rownames(tab.ship.pre_scaled), ]$GeneSymbolDB2


### Post-menopause
## Filter your SE to have only post-menopause samples as reported in the
## analysis plan. Remember to select men individuals inside women's age range.
ship_SE_post

## This model is the same than in the main model 
mod_ship_post <- model.matrix(formula(paste("~ sex + age + blood_cell_proportions + technical_variables", 
                                              paste(paste0("SV", 1:3), collapse = "+"))),
                            colData(ship_SE_post ))
lm.post_main <- lmFit(assay(ship_SE_post), mod_ship_post) %>%
  eBayes()
tab.ship.post_main <- topTable(lm.post_main, coef = 2, n = Inf, confint = TRUE)
tab.ship.post_main$chr <- as.character(seqnames(ship_SE_post[rownames(tab.ship.post_main), ]))
tab.ship.post_main$Symbol <-  rowData(ship_SE_post)[rownames(tab.ship.post_main), ]$GeneSymbolDB2

lm.post_scaled <- lmFit(t(scale(t(assay(ship_SE_post)))), mod_ship_post) %>%
  eBayes()
tab.ship.post_scaled <- topTable(lm.post_scaled, coef = 2, n = Inf, confint = TRUE)
tab.ship.post_scaled$chr <- as.character(seqnames(ship_SE_post[rownames(tab.ship.post_scaled), ]))
tab.ship.post_scaled$Symbol <-  rowData(ship_SE_post)[rownames(tab.ship.post_scaled), ]$GeneSymbolDB2


save(tab.ship.scaled, tab.ship.main, tab.ship.pre_main, tab.ship.pre_scaled, tab.ship.post_main, tab.ship.post_scaled,
 file = "results/ship_analysis.Rdata")
#' This file should be shared with us for the manuscript preparation.