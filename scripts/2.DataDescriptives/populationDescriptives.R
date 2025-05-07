###############################################################################
#' Create table with population descriptives
#'  docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.4 R
###############################################################################

library(tidyverse)
library(MASS)

## Load data ####
load("results/preprocessFiles/filtered_phenotypes.Rdata")

## Select phenotypes 
sel_phenos <- c("height", "hs_c_weight", "hs_c_bmi", "hs_bp_sys", "hs_bp_dia", 
                "hs_waist", "hs_midup_arm", "hs_w2h",  
                "hs_skf_sum2", "hs_asthma", "hs_eczema", "hs_algy_food", "hs_rhin_ly",
                "hs_ADHD_raw", "hs_correct_raven", "hs_Gen_Int",
                "hs_Cognit_raw", "hs_Gen_Ext", "hs_Hyper_raw", "h_edumc")

phenos <- sel_phenos[!sel_phenos %in% c("HelixID", "cohort", "smoking", "h_edumc", "mat_educ") ]
phenos_cat <- c("asthma", "eczema", "algy_food", "rhin")
phenos_count <- c("hs_correct_raven", "hs_Gen_Int", "hs_Cognit_raw",
                  "hs_Gen_Ext", "hs_Hyper_raw", "hs_ADHD_raw")
phenos_cont <- phenos[!phenos %in% c(phenos_cat, phenos_count, "hs_asthma",   
                                     "hs_eczema", "hs_algy_food", "hs_rhin_ly")]
names(phenos_cat) <- phenos_cat

#'##############################################################################
#' Compute phenotype descriptives
#'##############################################################################
all_phenos$cohort <- droplevels(all_phenos$cohort)

chisq.test(table(all_phenos$cohort, all_phenos$mat_educ))

summary(lm(age_sample_years ~ cohort, all_phenos))

## Define functions
getSum <- function(vec, type = "continuous"){
  if (type == "continuous"){
    c(median = median(vec, na.rm = TRUE), 
      range = quantile(vec, probs = c(0.25, 0.75), na.rm = TRUE))
  } else if (type == "categorical"){
    t <- table(vec, useNA = "ifany")
    data.frame(names = names(t),  tab = as.vector(t), 
               props = as.vector(prop.table(t)))
  }
}
getSumCohort <- function(tab, varname, type = "continuous"){
  tab$var <- tab[, varname]
  df <- tab %>%
    data.frame() %>%
    group_by(cohort)
  if (type == "continuous"){
    df %>% 
      summarize(median = median(var, na.rm = TRUE), 
                range.25 = quantile(var, probs = 0.25, na.rm = TRUE),
                range.75 = quantile(var, probs = 0.75, na.rm = TRUE)) %>%
      mutate(val = sprintf("%.2f (%.2f-%.2f)", median, range.25, range.75),
             cohort = as.character(cohort)) %>%
      dplyr::select(cohort, val) %>%
      rbind(c("pval", format(anova(lm(var ~ cohort, df))[1, 5], scientific = TRUE, digits = 3)))
  
  } else if (type == "categorical"){
    df %>%
      group_by(cohort, var) %>%
      summarize(N = n()) %>%
      group_by(cohort) %>%
      mutate(val = paste0(N, " (", round(N/sum(N)*100, 1), "%)")) %>%
      dplyr::select(cohort, var, val) %>%
      spread(var, val) 
    }
}
getCohortPvals <- function(tab, varname){
  tab$var <- tab[, varname]
  tab$cohort <- droplevels(tab$cohort)
 chisq.test(table(tab$var, tab$cohort))$p.value
}

getCohortPvalsCont <- function(tab, varname, model){
  tab$var <- tab[, varname]
  anova(lm(model, tab))
}


getCohortPvalsCount <- function(tab, varname, model){
  tab$var <- tab[, varname]
  anova(glm.nb(var ~ cohort, tab), 
  glm.nb(model, tab), 
  test = "Chisq")
}


getCohortPvalsCat <- function(tab, varname, model){
  tab$var <- factor(tab[, varname])
  anova(glm(var ~ cohort, tab, family = "binomial"), 
  glm(model, tab, family = "binomial"), 
  test = "Chisq")
}


pheno_vec_name <- c( height = "Height",
                     hs_c_weight = "Weight",
                     hs_waist = "Waist circumf.",
                     hs_midup_arm =	"Mid-upper arm circumf.",
                     hs_w2h	= "Waist-to-Height ratio",
                     hs_c_bmi = "BMI",                
                     hs_skf_sum2 = "Skinfold Thickness",
                   
                     hs_bp_sys = "Systolic BP",
                     hs_bp_dia = "Diastolic BP",
   
                     asthma = "Asthma",
                     eczema = "Eczema",
                     algy_food	= "Food Allergy",
                     rhin	= "Rhinitis",
                     
                     hs_correct_raven =	"Score Raven Test",
                     hs_Gen_Int	= "Internalizing scales",
                     hs_Gen_Ext	= "Externalizing scales",
                     hs_Cognit_raw = "Cognitive Problems",
                     hs_Hyper_raw = "Hyperactivity score",
                     hs_ADHD_raw  = "ADHD symptomatology"     
)


## Whole cohort descriptives
### Categorical variables
catVecs <- lapply(phenos_cat, function(x) all_phenos[, x])
catSums <- lapply(catVecs, getSum, type = "categorical")
names(catSums) <- phenos_cat
catSums <- lapply(phenos_cat, function(x){
  catSums[[x]]$names <- paste(x, catSums[[x]]$names)
  catSums[[x]]
})

### Create table
catTab <- Reduce(rbind, catSums)
#### Convert proportions to percentages
catTab$text <- paste0(sprintf("%d (%.2f", catTab$tab, catTab$props*100), "%)")

### Continuous variables
contVecs <- lapply(c(phenos_cont, phenos_count), function(x) all_phenos[, x])
contSums <- lapply(contVecs, getSum, type = "continuous")

#### Create table
contTab <- data.frame(Reduce(rbind, contSums))
#### Convert proportions to percentages
contVals <- sprintf("%.2f (%.2f-%.2f)", contTab$median, contTab$range.25., contTab$range.75.)
contTab <- data.frame(names = c(phenos_cont, phenos_count),
                      text = contVals)

allTab <- rbind(catTab[, c("names", "text")],
                contTab)

## Stratified by cohort ####
cohortCat <- lapply(phenos_cat, getSumCohort, tab = all_phenos, type = "categorical") 
names(cohortCat) <- phenos_cat

cohortCatTab <- lapply(phenos_cat, function(x){
  
  tab <- cohortCat[[x]][, -1]
  tab[is.na(tab)] <- "0 (0%)"
  res <- t(tab) %>%
    data.frame()
  colnames(res) <- cohortCat[[x]]$cohort
  rownames(res)[rownames(res) == "<NA>"] <- "NA"
  res$names = paste(x, rownames(res))
  res
}) %>%
  Reduce(f = rbind)
cohortCatTab$pval <- NA
cohortCatTab$age <- NA
cohortCatTab$mat_educ <- NA


cohortCont <- lapply(c(phenos_cont, phenos_count), getSumCohort, tab = all_phenos) %>%
  Reduce(f = function(x, y) left_join(x, y, by = "cohort"))
colnames(cohortCont)[-1] <- c(phenos_cont, phenos_count)

cohortContTab <- data.frame(t(cohortCont[, -1]))
colnames(cohortContTab) <- cohortCont$cohort
cohortContTab$names <- rownames(cohortContTab)



lapply(phenos_cont, getCohortPvalsCont, tab = all_phenos, model = formula(var ~ cohort + age_sample_years) )
pvals_cont <- lapply(phenos_cont, getCohortPvalsCont, tab = all_phenos, model = formula(var ~ cohort + age_sample_years + mat_educ) )
names(pvals_cont) <- phenos_cont

extra_cols_count <- tibble(names = phenos_count, 
  age = sapply(phenos_count, function(var) 
    getCohortPvalsCount(tab = subset(all_phenos, !is.na(age_sample_years)), varname = var,  model = formula(var ~ cohort + age_sample_years ) )$`Pr(Chi)`[2]),
    mat_educ = sapply(phenos_count, function(var) 
    getCohortPvalsCount(tab = subset(all_phenos, !is.na(mat_educ)), varname = var,  model = formula(var ~ cohort + mat_educ) )$`Pr(Chi)`[2])
) 

extra_cols_cont <- Reduce(rbind, 
  lapply(phenos_cont, function(x) tibble(names = x, age = pvals_cont[[x]]$`Pr(>F)`[2], mat_educ = pvals_cont[[x]]$`Pr(>F)`[3]))
)

extra_cols <- rbind(extra_cols_cont, extra_cols_count)
cohortContTab2 <- left_join(cohortContTab, extra_cols, by = "names")

allTab$varName <- pheno_vec_name[sapply(strsplit(allTab$names, " "), `[`, 1)]

cohortComb <- rbind(cohortContTab2, cohortCatTab) %>%
  left_join(x = allTab[, c("varName", "names", "text")], by = "names")  

# Replace NA by missing
## Table 1
cohortComb$names <- gsub("NA", "Missing", cohortComb$names)
colnames(cohortComb)[1:3] <- c("Phenotype", "", "Total")
write.table(cohortComb, file = "results/descriptives/phenotypes_descriptives_cohort.txt", col.names = TRUE, 
            sep = "\t", quote = FALSE, row.names = FALSE)



sapply(phenos_cat, getCohortPvals, tab = all_phenos)




lapply(phenos_cat, getCohortPvalsCat, tab = all_phenos, model = formula(var ~ cohort ) )
sapply(phenos_cat, function(var) getCohortPvalsCat(varname = var, tab = subset(all_phenos, !is.na(mat_educ)), model = formula(var ~ cohort +  mat_educ))$`Pr(>Chi)`[2] )
sapply(phenos_cat, function(var) getCohortPvalsCat(varname = var, tab = subset(all_phenos, !is.na(age_sample_years)), model = formula(var ~ cohort +  age_sample_years))$`Pr(>Chi)`[2] )
lapply(phenos_cat, getCohortPvalsCat, tab = subset(all_phenos, !is.na(mat_educ)), model = formula(var ~ cohort +  mat_educ + age_sample_years) )

#'##############################################################################
#' Compute covariates descriptives
#'##############################################################################


## Categorical variables ####
### Describe Cohort ("cohort") and maternal education
cat_cov_Vecs <- lapply(c("cohort", "mat_educ"), function(x) all_phenos[, x])
cat_cov_Sums <- lapply(cat_cov_Vecs, getSum, type = "categorical")

## Create table
cat_cov_Tab <- Reduce(rbind, cat_cov_Sums)
## Convert proportions to percentages
cat_cov_Tab$text <- paste0(sprintf("%d (%.2f", cat_cov_Tab$tab, cat_cov_Tab$props*100), "%)")

## Continuous variables ####
### Describe Age ("age_sample_years") and 
### Cell types ("NK", "Bcell", "CD4T", "CD8T", "Mono", "Neu")
cont_cov_Vars <- c("age_sample_years", "NK_6", "Bcell_6", "CD4T_6", "CD8T_6", "Mono_6", "Gran_6")
cont_cov_Vecs <- lapply(cont_cov_Vars, function(x) all_phenos[, x])
cont_cov_Sums <- lapply(cont_cov_Vecs, getSum, type = "continuous")

## Create table
cont_cov_Tab <- data.frame(Reduce(rbind, cont_cov_Sums))
## Convert proportions to percentages
cont_cov_Vals <- sprintf("%.2f (%.2f-%.2f)", cont_cov_Tab$median, cont_cov_Tab$range.25., cont_cov_Tab$range.75.)
cont_cov_Tab <- data.frame(names = cont_cov_Vars, text = cont_cov_Vals)

all_cov_Tab <- rbind(cat_cov_Tab[, c("names", "text")],
                cont_cov_Tab)


## Statistics stratified by Sex
getSumSex <- function(tab, varname, type = "continuous"){
  tab$var <- tab[, varname]
  df <- tab %>%
    data.frame() %>%
    group_by(e3_sex)
  if (type == "continuous"){
    df %>% 
      summarize(median = median(var, na.rm = TRUE), 
                range.25 = quantile(var, probs = 0.25, na.rm = TRUE),
                range.75 = quantile(var, probs = 0.75, na.rm = TRUE)) %>%
      mutate(val = sprintf("%.2f (%.2f-%.2f)", median, range.25, range.75)) %>%
      select(e3_sex, val)
    
  } else if (type == "categorical"){
    df %>%
      group_by(e3_sex, var) %>%
      summarize(N = n()) %>%
      group_by(e3_sex) %>%
      mutate(val = paste0(N, " (", round(N/sum(N)*100, 1), "%)")) %>%
      select(e3_sex, var, val) %>%
      spread(var, val)
  }
}
sexCat <- lapply(c("cohort", "mat_educ"), getSumSex, type = "categorical", tab = all_phenos) %>%
   Reduce(f = function(x, y) left_join(x, y, by = "e3_sex"), x = .)
sexCont <- lapply(cont_cov_Vars, getSumSex, tab = all_phenos) %>%
  Reduce(f = function(x, y) left_join(x, y, by = "e3_sex"), x = .) %>%
  left_join(sexCat, ., by = "e3_sex")
colnames(sexCont)[-c(1:10)] <- cont_cov_Vars

allvec <- all_cov_Tab$text
names(allvec) <- all_cov_Tab$names
sexComb <- rbind(allvec[names(sexCont)[-1]], sexCont[, -1]) %>%
  t %>%
  rbind(c(nrow(all_phenos), sum(all_phenos$e3_sex == "female"),  sum(all_phenos$e3_sex == "male")))
rownames(sexComb)[nrow(sexComb)] <- "Total"
colnames(sexComb) <- c("All", "Girls", "Boys")

## Sup Table 1
write.table(sexComb, file = "results/descriptives/PopDescrip_cohort.txt", col.names = TRUE, 
            sep = "\t", quote = FALSE, row.names = TRUE)

chisq.test(table(all_phenos$cohort, all_phenos$e3_sex))
chisq.test(table(all_phenos$mat_educ, all_phenos$e3_sex))


contPvals <- function(varname, tab){
  tab$var <- tab[, varname]
  summary(lm(var ~ e3_sex, tab))$coefficients[2, 4]
}
sapply(cont_cov_Vars, contPvals, tab = all_phenos)
