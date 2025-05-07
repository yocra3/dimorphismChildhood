###############################################################################
#' Create SummarizedExperiment for HELIX gene expression data
#'  docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.8 R
###############################################################################

## Load data and libraries ####
library(Biobase, verbose = FALSE)
library(SummarizedExperiment, verbose = FALSE)
library(sva, verbose = FALSE)
library(isva, verbose = FALSE)
library(limma, verbose = FALSE)    
library(SmartSVA, verbose = FALSE) 

## Load SummarizedExperiment
load("./data/transcriptome_subcohort_notfitr_inclsex_v3.RData")
gexp <- transcriptome_subcohort_notfitr_inclsex


## Remove probes with call rate < 85 -> >90 loses XIST
gexp_filt <- gexp[fData(gexp)$CallRate > 85, ]

## Select European individuals and individuals with cell counts
gexp_filt <- gexp_filt[, !is.na(gexp_filt$NK_6) & gexp_filt$h_ethnicity_c == "Caucasian"]


# Making SE
se <- makeSummarizedExperimentFromExpressionSet(gexp_filt)

## Change ranges to be TSS
gr <- rowRanges(se)
gr$TC_Start <- start(gr)
gr$TC_end <- end(gr)
start(gr) <- gr$TSS_Affy
end(gr) <- gr$TSS_Affy
rowRanges(se) <- gr

save(se, file = "results/preprocessFiles/Expression_SE_raw.RData")

