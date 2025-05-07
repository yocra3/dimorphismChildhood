#'##############################################################################
#' Compute descriptives of the cohorts
#' docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.3 R
#'##############################################################################

## Load data and packages ####
library(tidyverse)
library(SummarizedExperiment())

load("results/preprocessFiles/gexp_GSE43488_sex.Rdata")
load("results/preprocessFiles/gexp_gse36382_sex.Rdata")
load("results/preprocessFiles/gexp_gse33828_sex.Rdata")
load("results/preprocessFiles/gexp_gtex_sva_sex.Rdata")

## Sub Table 2
IQR <- function(vec) c(median = median(vec, na.rm = TRUE), 
      range = quantile(vec, probs = c(0.25, 0.75), na.rm = TRUE))

gtex_iqr <- IQR(sapply(strsplit(blood_gtex_filt$gtex.age, "-"), function(x) mean(as.numeric(x))))
gseiqrs <- lapply(list( gse43488_controls$age, gse33828_se$age, gse36382_se$age), IQR)
iqrs <- c(gseiqrs, list(gtex_iqr))
names(iqrs) <- c("GSE43488", "GSE33828", "GSE36382", "GTEx") 
iqrs$GSE36382 <- c(median = NA, iqrs$GSE36382)

iqr_text <- sapply(iqrs, function(x) 
    sprintf("%.2f (%.2f-%.2f)", x[1], x[2], x[3]))

propSex <- function(vec){

    sprintf("%d (%.1f%%)", sum(vec == "female"), mean(vec == "female")*100)
}

datasets <- list(gse43488_controls, gse33828_se, gse36382_se, blood_gtex_filt)
sum_tab <- tibble(N = sapply(datasets, ncol),
    age = iqr_text, sex = sapply(datasets, function(x) propSex(x$sex))) %>%
    t()
colnames(sum_tab) <- c("GSE43488", "GSE33828", "GSE36382", "GTEx")
write.table(sum_tab, file = "results/descriptives/descriptives_gexp_cohorts.txt", col.names = TRUE, 
            sep = "\t", quote = FALSE, row.names = FALSE)