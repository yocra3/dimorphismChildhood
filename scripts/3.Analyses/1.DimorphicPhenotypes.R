#'##############################################################################
#'##############################################################################
#' Definition of dimorphic phenotypes in HELIX
#'  docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.8 R
#'##############################################################################
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(ggExtra)
library(cowplot)
library(S4Vectors)
library(readxl)
library(ggcorrplot)
#library(forestploter)
library(MASS)
library(meta)

## Load data ####
load("results/preprocessFiles/filtered_phenotypes.Rdata")

## Process phenotype data ####
## Select phenotypes included in the manuscript
sel_phenos <- c("height", "hs_c_weight", "hs_c_bmi", "hs_bp_sys", "hs_bp_dia", 
                "hs_waist", "hs_midup_arm", "hs_w2h",  
                "hs_skf_sum2", "hs_asthma", "hs_eczema", "hs_algy_food", "hs_rhin_ly",
                "hs_ADHD_raw", "hs_correct_raven", "hs_Gen_Int",
                "hs_Cognit_raw", "hs_Gen_Ext", "hs_Hyper_raw", "h_edumc")


phenos <- sel_phenos[!sel_phenos %in% c("HelixID", "cohort", "smoking", "h_edumc", "mat_educ") ]

## Divide phenotypes depending on their statistical properties 
phenos_cat <- c("asthma", "eczema", "algy_food", "rhin")
phenos_count <- c("hs_ADHD_raw", "hs_Gen_Int", "hs_Cognit_raw",
                  "hs_Gen_Ext", "hs_Hyper_raw")
phenos_cont <- phenos[!phenos %in% c(phenos_cat, phenos_count, "hs_asthma", "hs_correct_raven",  
                                     "hs_eczema", "hs_algy_food", "hs_rhin_ly")]

## Define names for manuscript
pheno_vec_name <- c( height = "Height",
                     hs_c_weight = "Weight",
                     hs_waist = "Waist circumf.",
                     hs_midup_arm =	"Mid-upper arm circumf.",
                     hs_w2h	= "Waist-to-height ratio",
                     hs_c_bmi = "BMI",                
                     hs_skf_sum2 = "Skinfold thickness",
                   
                     hs_bp_sys = "Systolic BP",
                     hs_bp_dia = "Diastolic BP",
   
                     asthma = "Asthma",
                     eczema = "Eczema",
                     algy_food	= "Food allergy",
                     rhin	= "Rhinitis",
                     
                     hs_correct_raven =	"Non-verbal intelligence",
                     hs_Gen_Int	= "Internalizing scale",  
                     hs_Gen_Ext	= "Externalizing scale",
                     hs_Cognit_raw = "Inattention index",
                     hs_Hyper_raw = "Hyperactivity index",
                     hs_ADHD_raw  = "ADHD index"     
)

#'##############################################################################
#' Define dimorphic phenotypes
#'##############################################################################

## Check dispersion in count data
a <- lapply(c(phenos_count, "hs_correct_raven") , function(phe){
 
  model <- glm.nb(formula(paste(phe, " ~ e3_sex + cohort + age_sample_years  + mat_educ")), 
              all_phenos)
  residual_deviance <- model$deviance
  df <- model$df.residual
  dispersion_ratio <- residual_deviance / df
  cat("Phenotype:", phe, "\n")
  cat("Residual Deviance:", residual_deviance, "\n")
  cat("Degrees of Freedom:", df, "\n")
  cat("Dispersion Ratio:", dispersion_ratio, "\n")
  observed_variance <- var(all_phenos[[phe]], na.rm = TRUE)
  observed_mean <- mean(all_phenos[[phe]], na.rm = TRUE)
  vmr <- observed_variance / observed_mean

  cat("Variance-to-Mean Ratio:", vmr, "\n")

})

## Run association test depending on variable types
dim_phenos_cont <- lapply(phenos_cont, function(phe){
  message(phe)
  lm(formula(paste(phe, "~ e3_sex + cohort + age_sample_years + mat_educ")), all_phenos)
})

dim_phenos_cat <- lapply(phenos_cat, function(phe){
  message(phe)

  glm(formula(paste("factor(", phe, ", ordered = TRUE) ~ e3_sex + cohort + age_sample_years  + mat_educ")), 
              all_phenos, family = "binomial")
  
})
dim_phenos_count <- lapply(phenos_count, function(phe){
  message(phe)
  
  glm.nb(formula(paste(phe, " ~ e3_sex + cohort + age_sample_years  + mat_educ")), 
              all_phenos)
  
})

raven_test <- glm(hs_correct_raven ~ e3_sex + cohort + age_sample_years  + mat_educ,
              all_phenos, family = "poisson")

phenos_sex_assoc <- c(dim_phenos_cont, dim_phenos_cat, dim_phenos_count, list(raven_test))
names(phenos_sex_assoc) <- c(phenos_cont, phenos_cat, phenos_count, "hs_correct_raven")

### Phenotypes summary ####
### Group phenotypes by categories 
df_pheno_annot <- data.frame(var = names(pheno_vec_name), Phenotype = pheno_vec_name, 
                             Category = rep(c("Anthropometric phenotypes", "Cardiovascular Traits", "Atopic Traits", "Cognition and Neuro-behaviour"), c(7, 2, 4, 6))) %>%
  mutate(var_type = ifelse(var %in% phenos_cont, "Continous", 
                           ifelse(var %in% phenos_cat, "Categorical", "Count")))
df_pheno_annot$Category <- factor(df_pheno_annot$Category, levels = c("Anthropometric phenotypes", "Cardiovascular Traits", "Atopic Traits", "Cognition and Neuro-behaviour"))
rownames(df_pheno_annot) <- df_pheno_annot$var

df_pheno_assocs <- tibble(var = names(phenos_sex_assoc), 
                          coef = sapply(phenos_sex_assoc,  function(x) x$coefficients[2]),
                          coef_low = sapply(phenos_sex_assoc,  function(x) confint(x)[2, 1]),
                          coef_high = sapply(phenos_sex_assoc,  function(x) confint(x)[2, 2]),
                          pval = sapply(phenos_sex_assoc,  function(x) summary(x)$coefficients[2, 4])) %>%
  left_join(df_pheno_annot, by = "var") %>%
  mutate(FDR = p.adjust(pval, "fdr"),
         Significance = ifelse(FDR < 0.05, "Significant", "Non-significant"),
         Significance = factor(Significance, levels = c("Significant", "Non-significant"))) 

phenos.dim <- subset(df_pheno_assocs, FDR < 0.05)$var
save(df_pheno_assocs, file = "results/phenotypes/dimorphic_phenos.Rdata")

### Plot of all associations ####

### Figure 2A
plot_coef_anthro <- df_pheno_assocs %>%
  filter(Category == "Anthropometric phenotypes") %>%
  arrange(Phenotype) %>%
  mutate(Phenotype = factor(Phenotype, levels = c("Mid-upper arm circumf.", "Waist-to-height ratio", "Waist circumf.", "Skinfold thickness", "BMI", "Weight", "Height")),
         Direction = ifelse(Significance == "Significant", 
                            ifelse(coef > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
  ggplot(aes(x = Phenotype, y = coef, color = Direction)) +
  geom_point(size = 3) +
  scale_color_manual(breaks = c("Higher Boys", "Higher Girls", "No differences"),
                     values = c("blue", "maroon", "black")) +
  geom_errorbar(aes(x = Phenotype, ymin = coef_low, ymax = coef_high)) +
  theme_bw(20) +
    coord_flip() +
  ylab("Coefficient") +
  geom_hline(yintercept = 0) +
  ggtitle("Anthropometric phenotypes")  +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.text.y.right =  element_text(angle = 90, hjust = 0.5),
        axis.ticks.y.right = element_blank())


plot_coef_blood <- df_pheno_assocs %>%
  filter(Category == "Cardiovascular Traits") %>%
  arrange(Phenotype) %>%
  mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)),
         Direction = ifelse(Significance == "Significant", 
                            ifelse(coef > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
  ggplot(aes(x = Phenotype, y = coef, color = Direction)) +
  geom_point(size = 3) +
  scale_color_manual(breaks = c("Higher Boys", "Higher Girls", "No differences"),
                     values = c("blue", "maroon", "black")) +
  geom_errorbar(aes(x = Phenotype, ymin = coef_low, ymax = coef_high)) +
  theme_bw(20) +
    coord_flip() +
  ylab("Coefficient") +
  geom_hline(yintercept = 0) +
  ggtitle("Cardiovascular Traits")  +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.text.y.right =  element_text(angle = 90, hjust = 0.5),
        axis.ticks.y.right = element_blank())



plot_coef_immune <- df_pheno_assocs %>%
  filter(Category == "Atopic Traits") %>%
  arrange(Category, Phenotype) %>%
  arrange(Phenotype) %>%
  mutate(coef = exp(coef),
                coef_low = exp(coef_low),
              coef_high = exp(coef_high),
          Phenotype = factor(Phenotype, levels = unique(Phenotype)),
         Direction = ifelse(Significance == "Significant", 
                            ifelse(coef > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
  ggplot(aes(x = Phenotype, y = coef, color = Direction)) +
  geom_point(size = 3) +
  scale_color_manual(breaks = c("Higher Boys", "Higher Girls", "No differences"),
                     values = c("blue", "maroon", "black")) +
  geom_errorbar(aes(x = Phenotype, ymin = coef_low, ymax = coef_high)) +
  theme_bw(20) +
  geom_hline(yintercept = 1) +
  ylab("Odds Ratio") +
    coord_flip() +
  scale_y_log10() +
  ggtitle("Atopic Traits")  +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.text.y.right =  element_text(angle = 90, hjust = 0.5),
        axis.ticks.y.right = element_blank())

plot_coef_neuro <- df_pheno_assocs %>%
  filter(var_type == "Count") %>%
  arrange(Category, Phenotype) %>%
  mutate(coef = exp(coef),
          coef_low = exp(coef_low),
          coef_high = exp(coef_high),
          Phenotype = factor(Phenotype, 
          levels = c("Inattention index", "Hyperactivity index", "ADHD index", 
            "Internalizing scale", "Externalizing scale", 
              "Non-verbal intelligence")),
         Direction = ifelse(Significance == "Significant", 
                            ifelse(coef > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
  ggplot(aes(x = Phenotype, y = coef, color = Direction)) +
  geom_point(size = 3) +
  scale_color_manual(breaks = c("Higher Boys", "Higher Girls", "No differences"),
                     values = c("blue", "maroon", "black")) +
  geom_errorbar(aes(x = Phenotype, ymin = coef_low, ymax = coef_high)) +
  theme_bw(20)+
  ylab("Risk Ratio") +
  ggtitle("Cognition and Neurobehaviour")  +
  coord_flip() +
  geom_hline(yintercept = 1) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.text.y.right =  element_text(angle = 90, hjust = 0.5),
        axis.ticks.y.right = element_blank()) 


legend <- get_legend(df_pheno_assocs %>%
  mutate(Direction = ifelse(Significance == "Significant", 
                            ifelse(coef > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
  ggplot(aes(x = Phenotype, y = coef, color = Direction)) +
                       theme_bw(20) +
                       geom_point() +
  scale_color_manual(name = "Effect",
                     breaks = c("Higher Boys", "Higher Girls", "No differences"),
                     values = c("blue", "maroon", "black"))
)
  
plot_coef_all <- plot_grid(
        plot_grid(
          plot_coef_anthro, 
          plot_coef_blood, 
          ncol = 1, rel_heights = c(2.5, 1)
        ),
        plot_grid(
        plot_coef_neuro,
        plot_coef_immune,
        ncol = 1, rel_heights = c(1.5, 1)
        )
,
  legend,
  ncol = 3, rel_widths = c(3, 3, 1)
)
png("figures/dimorphic_phenotypes_coef.png", width = 5000, height = 2700, res = 300)
plot_coef_all
dev.off()


group_by(df_pheno_assocs, Category) %>%
  summarize(Total = n(),
            Significative = sum(FDR < 0.05))

## Compute meta-analysis
cohorts <- unique(all_phenos$cohort)
names(cohorts) <- cohorts
cohorts_name <- recode(cohorts, SAB = "INMA")

dim_phenos_cont_cohort <- lapply(phenos_cont, function(phe){
  message(phe)
  lapply(cohorts, function(coh){
    lm(formula(paste(phe, "~ e3_sex + age_sample_years + mat_educ")), 
       all_phenos, subset = cohort == coh)
  })
})
names(dim_phenos_cont_cohort) <- phenos_cont

dim_phenos_cont_meta <- lapply(dim_phenos_cont_cohort, function(models){
  te <- sapply(models, function(m) summary(m)$coefficients[2, 1])
  sete <- sapply(models, function(m) summary(m)$coefficients[2, 2])
  n <- sapply(models, nobs)
  metagen(
    TE = te,
    seTE = sete,
    n.e = n,
    studlab = as.character(cohorts_name[cohorts]),
    sm = "MD",
    method.tau = "REML"
   )
})

phenos_count_cohort <- lapply(phenos_count, function(phe){
  message(phe)
  lapply(cohorts, function(coh){
      glm.nb(formula(paste(phe, " ~ e3_sex + age_sample_years  + mat_educ")), 
              all_phenos, subset = cohort == coh)
  })
})
names(phenos_count_cohort) <- phenos_count

phenos_count_meta <- lapply(phenos_count_cohort, function(models){
  te <- sapply(models, function(m) summary(m)$coefficients[2, 1])
  sete <- sapply(models, function(m) summary(m)$coefficients[2, 2])
  n <- sapply(models, nobs)
  metagen(
    TE = te,
    seTE = sete,
    n.e = n,
    studlab = as.character(cohorts_name[cohorts]),
    sm = "RR",
    method.tau = "REML"
  )
})

phenos_cat_cohort <- lapply(phenos_cat, function(phe){
  message(phe)
  lapply(cohorts, function(coh){
    glm(formula(paste("factor(", phe, ", ordered = TRUE) ~ e3_sex + age_sample_years + mat_educ")),
        all_phenos, subset = cohort == coh, family = "binomial")
  })
})
names(phenos_cat_cohort) <- phenos_cat

phenos_cat_meta <- lapply(names(phenos_cat_cohort), function(phe){
  models <- phenos_cat_cohort[[phe]]
  te <- sapply(models, function(m) summary(m)$coefficients[2, 1])
  sete <- sapply(models, function(m) summary(m)$coefficients[2, 2])
  n <- sapply(models, nobs)
  if (phe == "algy_food" && "BIB" %in% names(te)) {
    te["BIB"] <- NA
    sete["BIB"] <- NA
    n["BIB"] <- NA
  }
  
  metagen(
    TE = te,
    seTE = sete,
    n.e = n,
    studlab = as.character(cohorts_name[cohorts]),
    sm = "OR",
    method.tau = "REML"
  )
})
names(phenos_cat_meta) <- names(phenos_cat_cohort)

raven_cohort <- lapply(cohorts, function(coh){
  glm(hs_correct_raven ~ e3_sex + age_sample_years + mat_educ,
      all_phenos, subset = cohort == coh, family = "poisson")
})
names(raven_cohort) <- cohorts

raven_meta <- metagen(
  TE = sapply(raven_cohort, function(m) summary(m)$coefficients[2, 1]),
  seTE = sapply(raven_cohort, function(m) summary(m)$coefficients[2, 2]),
  n.e = sapply(raven_cohort, nobs),
  studlab = as.character(cohorts_name[cohorts]),
  sm = "RR",
  method.tau = "REML"
)


plot_forest <- function(meta, var){
  title <- pheno_vec_name[var]
  grid::grid.grabExpr(
    forest(meta,
           smlab = title,
           prediction = TRUE,
           leftcols = c("studlab", "n.e"),
           leftlabs = c("Cohort", "N"),
           rightcols = c("effect", "ci"),
           test.overall.common = TRUE,
           test.overall.random = TRUE,
           label.left = "Higher Girls",
           label.right = "Higher Boys",
           addrows.below.overall = 3)
  )
}

### Figure 2B

meta_objs <- c(dim_phenos_cont_meta, phenos_count_meta, phenos_cat_meta, list(hs_correct_raven= raven_meta))

main_phenos <- c("height", "hs_c_weight", "hs_skf_sum2", "hs_bp_sys", 
                 "hs_ADHD_raw", "hs_Cognit_raw", "hs_Gen_Ext", "hs_Hyper_raw")


main_forest <- lapply(main_phenos, function(phe){
  plot_forest(meta_objs[[phe]], phe)
})


forest_panel <- plot_grid(plotlist = main_forest, nrow = 3)
png("figures/dimorphic_phenos_forest.png", width = 6500, height = 3500, res = 300)
forest_panel
dev.off()

sup_phenos <- c("hs_c_bmi", "hs_waist", "hs_midup_arm", "hs_w2h", "hs_bp_dia",
              "asthma", "eczema", "algy_food", "rhin",
              "hs_correct_raven", "hs_Gen_Int")

sup_forest <- lapply(sup_phenos, function(phe){
  plot_forest(meta_objs[[phe]], phe)
})

png("figures/sup_phenos_forest.png", height = 5000, width = 6500, res = 300)
plot_grid(plotlist = sup_forest, nrow = 4)
dev.off()



## Figure 2
png("figures/dimorphic_phenos_panel.png", width = 6000, height = 6000, res = 300)
plot_grid(plot_coef_all, forest_panel, rel_heights = c(1, 1.7), nrow = 2, labels = "AUTO")
dev.off()



## Compute correlation between phenotypes ####
phenos_cor_mat <- all_phenos[, df_pheno_assocs$var]

## Convert categorical phenos to numeric
phenos_cor_mat[, subset(df_pheno_assocs, Category == "Atopic Traits")$var] <- phenos_cor_mat[, subset(df_pheno_assocs, Category == "Atopic Traits")$var] == "Yes"
sel_phenos_cor <- cor(phenos_cor_mat, use = "complete")
colnames(sel_phenos_cor) <- rownames(sel_phenos_cor) <- pheno_vec_name[colnames(sel_phenos_cor)]

pheno_cor <- ggcorrplot(sel_phenos_cor, method = "circle", hc.order = TRUE) 
pheno_l <- ggplot_build(pheno_cor)
cat_order <- pheno_l$layout$panel_params[[1]]$y.sec$limits

text_cols <- RColorBrewer::brewer.pal(4, "Set1")
text_cols[2] <- "darkgrey"
df_pheno_annot$color <- text_cols[as.numeric(df_pheno_annot$Category)]

col_df <- df_pheno_annot
rownames(col_df) <- df_pheno_annot$Phenotype

## Sup Figure 1
pheno_cor <- pheno_cor +
  ggplot2::theme(axis.text.x = ggplot2::element_text(color = col_df[cat_order, "color"]), 
    axis.text.y =  ggplot2::element_text(color = col_df[cat_order, "color"]))

png("figures/selected_phenos_cor.png", width = 3000, height = 3000, res = 300)
pheno_cor
dev.off()


