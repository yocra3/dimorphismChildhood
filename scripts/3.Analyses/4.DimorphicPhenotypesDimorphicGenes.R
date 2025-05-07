#'##############################################################################
#' Phenotypes vs genes (Mediation analysis)
#' docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.8 R
#'##############################################################################

## Load data and packages ####
library(mediation)
library(limma)
library(readxl)
library(tidyverse)
library(poolr)
library(cowplot)
library(parallel)
library(MASS)
library(ggmanh)
library(ggrepel)
library(ggbreak)


## Load phenotypes
load("results/preprocessFiles/filtered_phenotypes.Rdata")
load("results/phenotypes/dimorphic_phenos.Rdata")
phenos.dim <- subset(df_pheno_assocs, FDR < 0.05)$var

all_phenos$cohort <- droplevels(all_phenos$cohort)


## Load Gene expression
load("./results/preprocessFiles/Expression_SE_raw_filtered_SV.RData")

se.filt$cohort <- droplevels(se.filt$cohort)

### Remove genes in chrY and genes 
se.noY <- se.filt[seqnames(se.filt) != "chrY", ]


## Load dimorphic genes
load("results/genes/HELIX_dimorphic_genes.Rdata")
load("results/genes/children_analysis.Rdata")
load("results/genes/adult_analysis.Rdata")

child_comb_sig <- subset(child_comb, adj.P.Val.com < 0.05)
dim.gexp <- child_comb_sig$TC

## Load genes annotation
xci_genes <- read_xlsx("data/XCI_genes.PMID29022598.xlsx", skip = 1) 
xci_genes <- xci_genes[, 1:7]
colnames(xci_genes) <- c("Symbol", "Ensembl", "Chr", "Start", "End", "Type", "XCI")

## Add SVs to all_phenos df
all_phenos_sv <- left_join(all_phenos, 
  select(colData(se.noY) %>% data.frame(), HelixID, starts_with("SV")), 
  by = "HelixID")
rownames(all_phenos_sv) <- all_phenos_sv$HelixID

## Dimorphic phenotypes vs dimorphic genes ####
pheno_assocs <- lapply(phenos.dim, function(phe){
  message(phe)
  mod <- model.matrix(formula(paste("~", phe, "+ cohort +  age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6 +",
  paste(paste0("SV", 1:19), collapse = "+"))), 
                      all_phenos_sv)
  lmfit <- lmFit(assay(se.noY[, rownames(mod)]), mod) %>% eBayes()
  topTable(lmfit, coef = 2, n = Inf)
})

pheno_assocs.dim <- lapply(pheno_assocs, function(x){
  
  sel <- x[dim.gexp,]
  sel$FDR <- p.adjust(sel$P.Value, method = "fdr")
  arrange(sel, P.Value)
})
names(pheno_assocs.dim) <- phenos.dim

### Manhattan ####
### Sup Fig 3
ntest_gene_dim <- meff(cor(t(assay(se.noY[dim.gexp, ]))), method = "galwey")


plotManhattan <- function(tab, title){
  tab$Position <- start(rowRanges(se.noY[rownames(tab, )]))
  tab$TC <- rownames(tab)
  tab$chr <- as.character(seqnames(se.noY[rownames(tab), ]))
  tab$chr <- factor(tab$chr, levels = paste0("chr", c(1:22, "X")))
  tab$chrNum <- factor(gsub("chr", "", tab$chr), levels = c(1:22, "X"))
  tab$Genes <- rowData(se.noY[tab$TC, ])$GeneSymbolDB2
  tab$Genes <- ifelse(tab$P.Value < 0.05/ntest_gene_dim, tab$Genes, "")
  tab$Genes[tab$Genes == "NA"] <- ""
  tab$pval <- ifelse(tab$P.Value < 1e-200, 1e-200, tab$P.Value)
  
  
  man_plot <- manhattan_plot(x = tab, pval.colname = "pval", 
                             signif = c(0.05/ntest_gene_dim),
                             chr.colname = "chrNum", pos.colname = "Position",
                             rescale = FALSE,
                             plot.title = title,
                             label.colname = "Genes") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0, 8)

}



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
                     hs_Cognit_raw = "Innatention index",
                     hs_Hyper_raw = "Hyperactivity index",
                     hs_ADHD_raw  = "ADHD index"     
)

man_order <- c("height", "hs_bp_sys", "hs_skf_sum2", "hs_Gen_Ext",
   "hs_ADHD_raw", "hs_Cognit_raw", "hs_Hyper_raw")

manhattans <- Map(plotManhattan, pheno_assocs.dim[man_order], pheno_vec_name[man_order])

png("figures/manhattan_phenotype.png", res = 300, height = 3000, width = 6000)
plot_grid(plotlist = manhattans, ncol = 3)
dev.off()

### Test mediation ####
pheno_gene.tab2 <- Reduce(rbind, lapply(phenos.dim, function(x) {
  sel <- subset(pheno_assocs.dim[[x]], P.Value < 0.05/ntest_gene_dim)
  if(nrow(sel) == 0){
    return(data.frame(TC = character(0), pheno = character(0), coef_pheno_tc = numeric(0)))
  }
  data.frame(TC = rownames(sel), pheno = x, coef_pheno_tc = sel$logFC)
}))
length(unique(pheno_gene.tab2$TC))
table(pheno_gene.tab2$pheno)

pheno_gene.tab_filt <- subset(pheno_gene.tab2, !pheno %in% c("hs_Hyper_raw", "hs_Cognit_raw"))
length(unique(pheno_gene.tab_filt$TC))
table(pheno_gene.tab_filt$pheno)


runMediations <- function(pheno, tc){
  pheno <- as.character(pheno)
  tc <- as.character(tc)
  df <- data.frame(pheno = all_phenos[, pheno],
                     tc = as.numeric(assay(se.noY[tc, ])))
  df <- cbind(df, colData(se.noY)) %>%
    filter(!is.na(pheno))

  mod.med <- lm(formula(paste("tc ~ e3_sex + cohort +  age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6 +",
  paste(paste0("SV", 1:19), collapse = "+"))), df)
    
  if (pheno %in% c("hs_ADHD_raw", "hs_Hyper_raw", "hs_Cognit_raw")){
     mod.out <-  glm.nb(pheno ~ tc + e3_sex + cohort + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, 
              df)
  } else {
    mod.out <- lm(pheno ~ tc + e3_sex + cohort + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, df)
  }
  med <- mediate(mod.med, mod.out, treat = "e3_sex", mediator = "tc", sims = 800)
  med
}

pheno_gene.tab2$coef_tc_sex <- tab.gexp[pheno_gene.tab2$TC, "logFC"]
pheno_gene.tab2 <- left_join(pheno_gene.tab2, select(df_pheno_assocs, var, coef), 
  by = c("pheno" = "var"))
pheno_gene.tab2 <- mutate(pheno_gene.tab2, 
                         coherent = ifelse(sign(coef_tc_sex) * sign(coef_pheno_tc) == sign(coef ), "Coherent", "Discordant")
)

pheno_gene.coherent2 <- subset(pheno_gene.tab2, coherent == "Coherent")
mediation.list2 <- mclapply(seq_len(nrow(pheno_gene.coherent2)), function(i) {
  runMediations(pheno = pheno_gene.coherent2[i, ]$pheno, tc = pheno_gene.coherent2[i, ]$TC)
}, mc.cores = 15)


mediation.df2 <- as_tibble(pheno_gene.coherent2) %>%
  mutate(med.prop = sapply(mediation.list2, function(x) x$n.avg),
         med.pvalue = sapply(mediation.list2, function(x) x$n.avg.p),
         med.fdr = p.adjust(med.pvalue, "fdr"))
mediation.df2$Symbol <- rowData(se.noY[mediation.df2$TC, ])$GeneSymbolDB2
save(mediation.list2, mediation.df2, file = "results/mediation/mediation_phenos_metachild.Rdata")


n_gene <- table(mediation.df2$TC)
subset(mediation.df2, TC %in% names(n_gene[n_gene > 1]))  %>%
  arrange(TC) %>% data.frame() %>% select(TC, pheno, Symbol)

## Sup Table 4
mediation.df2_out <- mutate(mediation.df2, 
  Phenotype = pheno_vec_name[pheno],
  coef_pheno_sex = coef, 
  med_prop = med.prop,
  med_pvalue = med.pvalue,
  med_fdr = med.fdr) %>%
  select(TC, Symbol, Phenotype, starts_with("coef"), starts_with("med"))
write.table(mediation.df2_out, file = "results/mediation/mediation_coding.txt", 
  quote = FALSE, row.names = FALSE, sep = "\t")


mediation_sig2 <- subset(mediation.df2, med.fdr < 0.05)

write.table(mediation_sig2, file = "results/mediation/mediation_coding_sig.txt", 
  quote = FALSE, row.names = FALSE)

subset(child_comb, Symbol %in% mediation_sig2$Symbol) %>% data.frame() %>% 
  select(Symbol, starts_with(c("logFC", "P.Value")))

subset(tab.adult, Symbol %in% mediation_sig2$Symbol) %>% data.frame() %>% 
  select(Symbol, starts_with(c("logFC", "P.Value")))

comb_pheno <- t(assay(se.noY[mediation_sig2$TC,])) %>%  
  as_tibble() %>%
  mutate(HelixID = colnames(se.noY)) %>%
  inner_join(all_phenos, by = "HelixID")

raw_mod <- summary(lm(hs_skf_sum2 ~ e3_sex + cohort + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, comb_pheno))
med_mod <- summary(lm(hs_skf_sum2 ~ TC0X000980.hg.1 + TC04001283.hg.1 + TC06004128.hg.1  + TC0X001077.hg.1 + e3_sex + cohort + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, comb_pheno))
1 - med_mod$coefficients[6, 1]/ raw_mod$coefficients[2, 1]

### Gene descriptives ####
mediator_genes <- unique(mediation_sig2$TC)

mediator_genes_df <- data.frame(TC = mediator_genes) %>%
  mutate(chr = as.character(seqnames(se.noY[TC, ])),
         Coding = ifelse(TC %in% codingTCs, "Coding", "non-Coding"), 
         Symbol = rowData(se.noY)[TC, ]$GeneSymbolDB2, 
         Symbol = ifelse(Symbol == "NA" | Symbol == "", TC, Symbol)) %>% 
  left_join(xci_genes[, c("Symbol", "Type", "XCI")],
            by = "Symbol")


## N genes
length(mediator_genes)

## N genes autosomic
table(mediator_genes_df$chr == "chrX")


### XCI
table(mediator_genes_df$XCI)


## Coding
table(mediator_genes_df$Coding)

### Plot mediations ####
## Sup Figures 4-7 
plot_mediations <- function(pheno,pheno_name,  TC, gene_name){

  comb_pheno <- tibble(tc = as.numeric(assay(se.noY[TC,])),
                      HelixID = colnames(se.noY)) %>%
    inner_join(all_phenos, by = "HelixID") %>%
    mutate(Sex = ifelse(e3_sex == "male", "Boys", "Girls"))
  comb_pheno$Pheno <- comb_pheno[[pheno]]

  plot_sex_tc <- ggplot(comb_pheno, aes(x = Sex, y = tc)) +
    geom_boxplot() +
    xlab("Sex") +
    geom_label(
      aes(x = Inf, y = Inf),  
      label = sprintf("P = %.1e", tab.gexp[TC, ]$P.Value), 
      hjust = 1, vjust = 1
    ) +
    ylab(paste(gene_name, "Expression")) +  
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))


  plot_sex_pheno <- ggplot(comb_pheno, aes(x = Sex, y = Pheno)) +
    geom_boxplot() +
    xlab("Sex") +
    ylab(pheno_name) +
    geom_label(
      aes(x = Inf, y = Inf),  
      label = sprintf("P = %.1e", (subset(df_pheno_assocs, var == pheno)$pval)), 
      hjust = 1, vjust = 1
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  plot_tc_pheno <- ggplot(comb_pheno, aes(x = tc, y = Pheno)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab(paste(gene_name, "Expression")) +
    ylab(pheno_name) +
    geom_label(
      aes(x = Inf, y = Inf),  
      label = sprintf("P = %.1e", pheno_assocs.dim_cod[[pheno]][TC, ]$P.Value), 
      hjust = 1, vjust = 1
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))


    title <- ggdraw() + 
    draw_label(paste("Sex -", gene_name, "-", pheno_name),
               fontface = 'bold',  hjust = 0.5  ) 
  plot_grid(title, 
    plot_grid(plot_sex_pheno, plot_sex_tc, plot_tc_pheno, nrow = 1, labels = LETTERS[1:3]),
  ncol = 1, rel_heights = c(1, 10))
}

mediation_plots <- lapply(seq_len(nrow(mediation_sig2)), function(i){

  pheno <- mediation_sig2$pheno[i]
  pheno_name <- pheno_vec_name[pheno]
  tc <- mediation_sig2$TC[i]
  symbol <- mediation_sig2$Symbol[i]

  plot_mediations(pheno, pheno_name, tc, symbol)
})


png("figures/mediation/efhc2_sk.png", height = 900, width = 3000, res = 300)
mediation_plots[[1]]
dev.off()

png("figures/mediation/cxcl5_sk.png", height = 900, width = 3000, res = 300)
mediation_plots[[2]]
dev.off()

png("figures/mediation/hla_sk.png", height = 900, width = 3000, res = 300)
mediation_plots[[3]]
dev.off()

png("figures/mediation/alas2_sk.png", height = 900, width = 3000, res = 300)
mediation_plots[[4]]
dev.off()
