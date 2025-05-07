###############################################################################
#' Select samples in common between phenotypes and gene expression in HELIX
#' docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.3 R
###############################################################################

library(dplyr)
library(SummarizedExperiment)
library(tidyr)

## Load data ####
pheno.raw <- read.csv("data/ap103_31ene.2023.csv")
rownames(pheno.raw) <- pheno.raw$HelixID

load("results/preprocessFiles/Expression_SE_raw.RData")
colnames(se) <- se$HelixID

## Filter SE to have only children <= 10 yo
se_age <- se[, se$age_sample_years <= 10 & se$cohort != "EDEN"]


## Select phenotypes after discussing with Mariona
sel_phenos <- c("hs_c_height", "hs_c_weight", "hs_c_bmi", "hs_bp_sys", "hs_bp_dia", 
                "hs_waist", "hs_midup_arm", "hs_w2h",  
                "hs_skf_sum2", "hs_asthma", "hs_eczema", "hs_algy_food", "hs_rhin_ly",
                "hs_ADHD_raw", "hs_correct_raven", "hs_Gen_Int",
                "hs_Cognit_raw", "hs_Gen_Ext", "hs_Hyper_raw", "h_edumc")

phenos.filt <- pheno.raw[intersect(se_age$HelixID, rownames(pheno.raw)), sel_phenos ]
se.filt <- se_age[,  rownames(phenos.filt)]
save(se.filt, file = "results/preprocessFiles/Expression_SE_raw_filtered.Rdata")

## Redefine categories
all_phenos <- cbind(phenos.filt, colData(se.filt))

## Convert height to cm
all_phenos$height <- all_phenos$hs_c_height*100

all_phenos$asthma <- ifelse(all_phenos$hs_asthma == "0", "Never", "Yes")
all_phenos$eczema <- ifelse(all_phenos$hs_eczema == "0", "Never", 
                            ifelse(all_phenos$hs_eczema == "1", "Yes", NA))
all_phenos$algy_food <- ifelse(all_phenos$hs_algy_food == "0", "Never", 
                            ifelse(all_phenos$hs_algy_food == "1", "Yes", NA))
all_phenos$rhin <- ifelse(all_phenos$hs_rhin_ly == "0", "Never", 
                            ifelse(all_phenos$hs_rhin_ly == "1", "Yes", NA))
all_phenos$mat_educ <- factor(c("Low", "Middle", "High")[all_phenos$h_edumc], 
  levels = c("Low", "Middle", "High"))
save(all_phenos, file = "results/preprocessFiles/filtered_phenotypes.Rdata")


