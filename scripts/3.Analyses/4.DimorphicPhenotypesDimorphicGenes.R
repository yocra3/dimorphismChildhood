#'##############################################################################
#' Phenotypes vs genes (Mediation analysis)
#' docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.5 R
#'##############################################################################

## Load data and packages ####
#library(mediation)
library(limma)
library(SummarizedExperiment)
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
load("data/gexpAnnotation.Rdata")

se.filt$cohort <- droplevels(se.filt$cohort)

### Remove genes in chrY and genes 
se.noY <- se.filt[seqnames(se.filt) != "chrY", ]


## Load dimorphic genes
load("results/genes/HELIX_cohort_genes_summary.Rdata")
load("results/genes/array_analysis.Rdata")
load("results/genes/gtex_analysis.Rdata")
load("results/genes/GenerationR_analysis_2026.Rdata")

dim.gexp <- subset(gene_cohort_summary, dimorphic == "Dimorphic")$TC

## Load genes annotation
xci_genes <- read_xlsx("data/XCI_genes.PMID29022598.xlsx", skip = 1) 
xci_genes <- xci_genes[, 1:7]
colnames(xci_genes) <- c("Symbol", "Ensembl", "Chr", "Start", "End", "Type", "XCI")

## Add SVs to all_phenos df
all_phenos_sv <- left_join(all_phenos, 
  dplyr::select(colData(se.noY) %>% data.frame(), HelixID, starts_with("SV")), 
  by = "HelixID")
rownames(all_phenos_sv) <- all_phenos_sv$HelixID

## Dimorphic phenotypes vs dimorphic genes ####
pheno_assocs <- lapply(phenos.dim, function(phe){
  message(phe)
  mod <- model.matrix(formula(paste("~", phe, "+ cohort +  age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6 +",
  paste(paste0("SV", 1:3), collapse = "+"))), 
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
save(pheno_assocs.dim, file = "results/genes_phenos/dimorphic_phenos.Rdata")
### Manhattan ####
### Sup Fig 3
# ntest_gene_dim <- meff(cor(t(assay(se.noY[dim.gexp, ]))), method = "galwey")


plotManhattan <- function(tab, title, p.thres = 1e-4){
  p.thres <- max(tab$P.Value[tab$FDR < 0.05], na.rm = TRUE)
  tab$Position <- start(rowRanges(se.noY[rownames(tab, )]))
  tab$TC <- rownames(tab)
  tab$chr <- as.character(seqnames(se.noY[rownames(tab), ]))
  tab$chr <- factor(tab$chr, levels = paste0("chr", c(1:22, "X")))
  tab$chrNum <- factor(gsub("chr", "", tab$chr), levels = c(1:22, "X"))
  tab$Genes <- rowData(se.noY[tab$TC, ])$GeneSymbolDB2
  tab$Genes <- ifelse(tab$P.Value < p.thres, tab$Genes, "")
  tab$Genes[tab$Genes == "NA"] <- ""
  tab$pval <- ifelse(tab$P.Value < 1e-200, 1e-200, tab$P.Value)
  
  
  man_plot <- manhattan_plot(x = tab, pval.colname = "pval", 
                             signif = c(p.thres),
                             chr.colname = "chrNum", pos.colname = "Position",
                             rescale = FALSE,
                             plot.title = title,
                             label.colname = "Genes") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0, 8)
  man_plot$layers[[2]] <- NULL ## Remove significance line
  man_plot
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
                     hs_Cognit_raw = "Inattention index",
                     hs_Hyper_raw = "Hyperactivity index",
                     hs_ADHD_raw  = "ADHD index"     
)

man_order <- c("height", "hs_c_weight", "hs_bp_sys", "hs_skf_sum2", "hs_Gen_Ext",
   "hs_ADHD_raw", "hs_Cognit_raw", "hs_Hyper_raw")

manhattans <- Map(plotManhattan, pheno_assocs.dim[man_order], pheno_vec_name[man_order])

png("figures/genes_pheno/manhattan_phenotype.png", res = 300, height = 3000, width = 6000)
plot_grid(plotlist = manhattans, ncol = 3)
dev.off()

## Combine all phenotype associations into a single table
pheno_assocs.dim_tab <- Map(function(tab, phe){
  tab$TC <- rownames(tab)
  tab$Trait <- pheno_vec_name[phe]
  tab
}, pheno_assocs.dim, names(pheno_assocs.dim)) %>%
  Reduce(rbind, .)

pheno_assocs.dim_tab <-  left_join(pheno_assocs.dim_tab, dplyr::select(expAnnot, TC, GeneSymbol_Affy, Coding, seqname), by = "TC") %>%
  mutate(Symbol = GeneSymbol_Affy, Chromosome = seqname) 

write.table(pheno_assocs.dim_tab %>% 
 dplyr::select(Trait, TC, Symbol, Chromosome, logFC, P.Value, FDR ),
   file = "results/genes_phenos/dimorphic_phenos_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)


### Explore genes ####
pheno_gene.tab2 <- Reduce(rbind, lapply(phenos.dim, function(x) {
  sel <- subset(pheno_assocs.dim[[x]], FDR < 0.05)
  if(nrow(sel) == 0){
    return(data.frame(TC = character(0), pheno = character(0), pheno_tc_coef = numeric(0), pheno_tc_pval = numeric(0), pheno_tc_FDR = numeric(0)))
  }
  data.frame(TC = rownames(sel), pheno = x, pheno_tc_coef = sel$logFC,  pheno_tc_pval = sel$P.Value, pheno_tc_FDR = sel$FDR)
}))
length(unique(pheno_gene.tab2$TC))
table(pheno_gene.tab2$pheno)

pheno_gene.tab2$Symbol <-  rowData(se.noY)[pheno_gene.tab2$TC, ]$GeneSymbolDB2
pheno_gene.tab2$chromosome <-  left_join(pheno_gene.tab2, data.frame(TC = rownames(se.noY), chromosome = as.character(seqnames(se.noY))), by = "TC")$chromosome


pheno_gene.tab_filt <- subset(pheno_gene.tab2, !pheno %in% c("hs_Hyper_raw", "hs_Cognit_raw"))
length(unique(pheno_gene.tab_filt$TC))
table(pheno_gene.tab_filt$pheno)
table(pheno_gene.tab_filt$TC)
nrow(pheno_gene.tab_filt)
table(pheno_gene.tab_filt$chromosome)
table(filter(gene_cohort_summary, TC %in% unique(pheno_gene.tab_filt$TC))$chr)

## Comparison with adults ####
tab.adult.genepheno <- filter(tab.array, GeneSymbol %in% unique(pheno_gene.tab2$Symbol) )
tab.gtex.genepheno <- filter(tab.gtex, Symbol %in% unique(pheno_gene.tab2$Symbol) )



## Run dimorphism in tissues ####
load("results/preprocessFiles/gtex_tissues_sva_sex.Rdata")

runTissueGTEx <- function(rse){

  sv_cols <- grep("^SV", colnames(colData(rse)), value = TRUE)
  mod_gtex <- model.matrix(formula(paste("~ sex + AGE + gtex.smtsisch + gtex.smrin +", paste(sv_cols, collapse = "+"))),
                                data = colData(rse))
  rse_noY <- rse[as.character(seqnames(rse)) %in% paste0("chr", c(1:22, "X")), ]

  ## Run Voom and limma
  gtex_voom <- voom(assay(rse_noY), mod_gtex)
  lm.gtex <- lmFit(gtex_voom, mod_gtex) %>%
    eBayes()
  tab.gtex <- topTable(lm.gtex, coef = 2, n = Inf, confint = TRUE)
  tab.gtex$Symbol <- rowData(rse_noY)[rownames(tab.gtex), "gene_name"]
  tab.gtex$chromosome <- as.character(seqnames(rse_noY[rownames(tab.gtex)]))
  tab.gtex
}
tissues_tabs <- lapply(gtex_tissues_sva, runTissueGTEx)
save(tissues_tabs, file = "results/genes_phenos/gtex_tissues_analyses.Rdata")

dim_symbol <- unique(pheno_gene.tab2$Symbol)
tissue_symbol <- lapply(dim_symbol, function(sym){
  tab_list <- lapply(names(tissues_tabs), function(tab_name){
    sub_df <- subset(tissues_tabs[[tab_name]], Symbol == sym) %>% dplyr::select(logFC, P.Value)
    if (nrow(sub_df) == 0){
      sub_df <- data.frame(logFC = NA, P.Value = NA, Tissue = tab_name)
    } else{
      sub_df$Tissue <- tab_name
    }
    sub_df
  })
  Reduce(rbind, tab_list)
})
names(tissue_symbol) <- dim_symbol


tissue_names <- c(
  "Adipose - Subcutaneous"                    = "Adipose_SC",
  "Adipose - Visceral (Omentum)"              = "Adipose_Vis",
  "Brain - Hippocampus"                       = "Brain_Hippocampus",
  "Brain - Frontal Cortex (BA9)"              = "Brain_PFC_BA9",
  "Brain - Anterior cingulate cortex (BA24)"  = "Brain_ACC_BA24",
  "Brain - Caudate (basal ganglia)"           = "Brain_Caudate",
  "Brain - Putamen (basal ganglia)"           = "Brain_Putamen",
  "Brain - Nucleus accumbens (basal ganglia)" = "Brain_NAcc",
  "Brain - Cerebellum"                        = "Brain_Cerebellum",
  "Liver"                                     = "Liver",
  "Muscle - Skeletal"                         = "Muscle",
  "Skin - Sun Exposed (Lower leg)"            = "Skin_Leg",
  "Skin - Not Sun Exposed (Suprapubic)"       = "Skin_Suprapubic"
)

tissue_tab <- lapply(dim_symbol[!dim_symbol %in% c("", NA)], function(gene){
  df <- tissue_symbol[[gene]]
  df$Tissue <- tissue_names[df$Tissue]
  df$Symbol <- gene
  df_wide <- pivot_wider(df, names_from = Tissue, values_from = c(logFC, P.Value))
}) %>% 
Reduce(rbind, .) 

pheno_gene.tissue <- left_join(pheno_gene.tab2, tissue_tab, by = "Symbol") %>%
  left_join(gene_cohort_summary %>% dplyr::select(TC, logFC,P.Value, chr), by = "TC") %>%
  left_join(tab.array %>% dplyr::select(GeneSymbol, logFC, P.Value), by = c("Symbol" = "GeneSymbol"), 
  suffix = c("", "_adult")) %>%
  left_join(tab.genR_main %>% dplyr::select(Symbol, logFC, P.Value), by = "Symbol",
  suffix = c("", "_genR")) %>%
  left_join(tab.gtex %>% dplyr::select(Symbol, logFC, P.Value), by = "Symbol", 
  suffix = c("_helix", "_gtex")) 

pheno_tissues_map <- list(
  hs_c_weight = c("Adipose_SC", "Adipose_Vis", "Brain_Hippocampus", "Liver", "Muscle"),
  hs_skf_sum2 = c("Adipose_SC", "Skin_Suprapubic", "Skin_Leg", "Muscle"),
  hs_ADHD_raw = c("Brain_PFC_BA9", "Brain_ACC_BA24", "Brain_Caudate", "Brain_Putamen", "Brain_NAcc", "Brain_Cerebellum")
)
pheno_tissues_map$hs_Cognit_raw <- pheno_tissues_map$hs_ADHD_raw
pheno_tissues_map$hs_Hyper_raw <- pheno_tissues_map$hs_ADHD_raw

sel_phenos <- unique(pheno_gene.tissue$pheno)

phenos_tabs <- lapply(sel_phenos, function(phe){

  gene_df <- filter(pheno_gene.tissue, pheno == phe) %>%
    dplyr::select(TC, pheno, pheno_tc_coef, pheno_tc_pval, Symbol, chr, logFC_helix, P.Value_helix, logFC_genR, P.Value_genR,
    logFC_adult, P.Value_adult, logFC_gtex, P.Value_gtex, ends_with(pheno_tissues_map[[phe]]))
}) 
names(phenos_tabs) <- sel_phenos

for (phe in sel_phenos){
  out_tab <- phenos_tabs[[phe]] %>%
    filter(!is.na(Symbol) & Symbol != "NA") %>%
    mutate(pheno = pheno_vec_name[phe]) %>%
    dplyr::select(pheno, TC, Symbol, chr, everything()) %>%
    rename_with(~ gsub("gtex", "Whole_Blood", .x)) %>%
    rename_with(~ gsub("^pheno_tc", "Trait_TC", .x)) %>%
    rename(Trait = pheno, Chromosome = chr)
  write.table(out_tab, file = paste0("results/genes_phenos/", phe, "_gene_table.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE)
}

## Genes dimorphic in adults
lapply(phenos_tabs, function(tab){
  filter(tab, P.Value_adult < 0.05 & sign(logFC_helix) == sign(logFC_adult) )$Symbol
}) %>% unlist() %>% unique() %>% length()

## Genes dimorphic in specific tissues
lapply(phenos_tabs[1:2], function(tab){
  tissue_names <- c("Adipose_SC", "Adipose_Vis", "Liver", "Muscle")
  matches <- sapply(tissue_names, function(t){
    pcol <- paste0("P.Value_", t)
    lcol <- paste0("logFC_", t)
    if (!all(c(pcol, lcol) %in% colnames(tab))) return(rep(FALSE, nrow(tab)))
    tab[[pcol]] < 0.05 & sign(tab[[lcol]]) == sign(tab$logFC_helix)
  })
  tab$TC[apply(matches, 1, any, na.rm = TRUE)]
}) %>% unlist() %>% unique() %>% length()

lapply(phenos_tabs[2], function(tab){
  tissue_names <- c("Skin_Leg", "Skin_Suprapubic")
  matches <- sapply(tissue_names, function(t){
    pcol <- paste0("P.Value_", t)
    lcol <- paste0("logFC_", t)
    if (!all(c(pcol, lcol) %in% colnames(tab))) return(rep(FALSE, nrow(tab)))
    tab[[pcol]] < 0.05 & sign(tab[[lcol]]) == sign(tab$logFC_helix)
  })
  tab$TC[apply(matches, 1, any, na.rm = TRUE)]
}) %>% unlist() %>% unique() %>% length()

## Genes dimorphic in brain tissues
lapply(phenos_tabs[c("hs_ADHD_raw", "hs_Cognit_raw", "hs_Hyper_raw")], function(tab){
  tissue_names <- c("Brain_PFC_BA9", "Brain_ACC_BA24", "Brain_Caudate", "Brain_Putamen", "Brain_NAcc", "Brain_Cerebellum")
  matches <- sapply(tissue_names, function(t){
    pcol <- paste0("P.Value_", t)
    lcol <- paste0("logFC_", t)
    if (!all(c(pcol, lcol) %in% colnames(tab))) return(rep(FALSE, nrow(tab)))
    tab[[pcol]] < 0.05 & sign(tab[[lcol]]) == sign(tab$logFC_helix)
  })
  tab$TC[apply(matches, 1, any, na.rm = TRUE)]
}) %>% unlist() %>% unique() %>% length()

## chrX genes coherent with HELIX in every tissue (excluding adults, GenR and Whole Blood)
lapply(phenos_tabs, function(tab){
  tab_x <- filter(tab, chr == "chrX")
  logfc_cols <- grep("^logFC_", colnames(tab_x), value = TRUE)
  tissue_names <- setdiff(gsub("^logFC_", "", logfc_cols), c("helix", "genR", "adult", "gtex"))
  matches <- sapply(tissue_names, function(t){
    tab_x[[paste0("P.Value_", t)]] < 0.05 & sign(tab_x[[paste0("logFC_", t)]]) == sign(tab_x$logFC_helix)
  })
  if (is.null(dim(matches))) matches <- matrix(matches, ncol = length(tissue_names))
  tab_x$TC[apply(matches, 1, function(x) all(!is.na(x)) && all(x))]
}) %>% unlist() %>% unique() %>% length()

labs_tissues <- c(
  "helix" = "HELIX",
  "genR" = "Generation R",
  "adult" = "Adults",
  "gtex" = "Whole Blood",
  "Adipose_SC" = "Adipose Sub.",
  "Adipose_Vis" = "Adipose Visc.",
  "Liver" = "Liver",
  "Muscle" = "Muscle",
  "Skin_Leg" = "Skin (Leg)",
  "Skin_Suprapubic" = "Skin (Suprapubic)",
  "Brain_Hippocampus" = "Brain (Hippocampus)",
  "Brain_PFC_BA9" = "Brain (Frontal Cortex)",
  "Brain_ACC_BA24" = "Brain (Anterior cingulate cortex)",
  "Brain_Caudate" = "Brain (Caudate)",
  "Brain_Putamen" = "Brain (Putamen)",
  "Brain_NAcc" = "Brain (Nucleus accumbens)",
  "Brain_Cerebellum" = "Brain (Cerebellum)"
)

plot_gene_concordance <- function(tab, main){
  tab %>%
    dplyr::select(Symbol, starts_with("logFC"), starts_with("P.Value_"), chr) %>%
    left_join(dplyr::select(xci_genes, Symbol, XCI), by = "Symbol") %>%
    mutate(chromosome = ifelse(chr == "chrX", "chrX", "Autosomic"),
    Gene = ifelse(!is.na(XCI) & XCI == "escape", paste0(Symbol, "*"), Symbol)) %>%
    pivot_longer(cols = matches("^(logFC|P\\.Value)_"), names_to = c(".value", "Tissue"), 
    names_pattern = "(logFC|P\\.Value)_(.+)") %>%
    mutate(Tissue = gsub("logFC_", "", Tissue),
    Tissue = labs_tissues[Tissue],
    Tissue = factor(Tissue, levels = c("HELIX", "Generation R", "Adults", "GTEx", setdiff(unique(Tissue), c("HELIX", "Generation R", "Adults", "GTEx")))),
    significance = case_when(
        P.Value < 1e-3 ~ "***",
        P.Value < 0.01 ~ "**",
        P.Value < 0.05 ~ "*",
        TRUE ~ ""
      ),
    Direction = ifelse(logFC > 0, "Boys", "Girls"),) %>%
    ggplot(aes(x = Tissue, y = Gene, fill = logFC, color = Direction)) +
    geom_tile() +
    ggtitle(main) +
    geom_text(aes(label = significance), size = 3,key_glyph = "text") +
    scale_fill_gradient2(low = "maroon", mid = "white", high = "blue", midpoint = 0,
    limits = c(-1, 1),  oob = scales::squish) +
    scale_color_manual(values = c("Boys" = "blue", "Girls" = "maroon"),
      na.value = "grey",
      na.translate = FALSE,
      guide = guide_legend(
    override.aes = list(fill = "white", label = "*")
  )) +
    theme_bw() +
    facet_grid(chromosome ~., scales = "free", space = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
) 
}

## Weight
png("figures/genes_pheno/heatmap_weight_tissues.png", width = 1500, height = 1600, res = 300)
phenos_tabs$hs_c_weight %>%
  filter(!Symbol %in% c("NA", "PRKX;RNU6-146P")) %>%
  plot_gene_concordance(main = "Weight")
dev.off()


weight_genes <- phenos_tabs$hs_c_weight %>%
  filter(Symbol %in% c("ZC3H12D", "TRAK2", "TCP11L2", 
  "PSME4", "PGM2L1", "PCBP2", "LEPROTL1", "IFIT1B", "ESPN", "EBF1", "DDX18", "CXCL5", "CTSE",
  "AP2A1", "UBA1", "STS", "MED14", "EIF2S3", "EIF1AX", "EFHC2", "CA5BP1", "ALAS2" )) %>%
  plot_gene_concordance(main = "Weight") +
  theme(legend.position = "none")


png("figures/genes_pheno/heatmap_weight_tissues_selGenes.png", width = 1500, height = 1500, res = 300)
weight_genes
dev.off()

## Skinfold Thickness
png("figures/genes_pheno/heatmap_skf_tissues.png", width = 1500, height = 1800, res = 300)
phenos_tabs$hs_skf_sum2 %>%
  filter(!Symbol %in% c("NA", "MIR4667;TESK1", "MID2;BC070370")) %>%
  plot_gene_concordance(main = "Skinfold thickness")
dev.off()

skf_genes <- phenos_tabs$hs_skf_sum2 %>%
  filter(Symbol %in% c("TRAK2", "TCP11L2", "SWT1", "SORBS3", "SETBP1", "RGCC", "PSME4", "PMAIP1", 
  "PITHD1", "PGM2L1", "LEPROTL1", "IFIT1B", "DDX18", "CXCL5", "CTAGE5", "CHST10", "CCDC30", "ZRSR2", "XIST", "XG", 
  "TSPAN7", "RPS4X", "P2RY8", "CDK16", "ALAS2") ) %>%
  plot_gene_concordance(main = "Skinfold thickness") +
  theme(legend.position = "none")

png("figures/genes_pheno/heatmap_skf_tissues_selGenes.png", width = 1500, height = 1500, res = 300)
skf_genes
dev.off()


## ADHD
png("figures/genes_pheno/heatmap_adhd_tissues.png", width = 1500, height = 1800, res = 300)
phenos_tabs$hs_ADHD_raw %>%
  filter(!Symbol %in% c("NA", "", "M55536;PANK3", "SCARNA9L", "LOC389906", "HDHD1")) %>%
  plot_gene_concordance(main = "ADHD")
dev.off()

adhd_genes <- phenos_tabs$hs_ADHD_raw %>%
  filter(Symbol %in% c("TTL", "KIAA0368", "FANCC", "ZFX", "XIST", "UBA1", "SYAP1", 
  "KDM6A", "KDM5C", "JPX", "EIF1AX", "DDX3X")) %>%
  plot_gene_concordance(main = "ADHD index")

png("figures/genes_pheno/heatmap_adhd_tissues_selGenes.png", width = 1500, height = 1500, res = 300)
adhd_genes
dev.off()


png("figures/genes_pheno/heatmap_innatention_tissues.png", width = 1500, height = 1800, res = 300)
phenos_tabs$hs_Cognit_raw %>%
  filter(!Symbol %in% c("NA", "", "M55536;PANK3", "TXLNG;AX746622", "D28359;INTS6-AS1", "LOC389906", "HDHD1")) %>%
  plot_gene_concordance(main = "Inattention index")
dev.off()

innatention_genes <- phenos_tabs$hs_Cognit_raw %>%
  filter(Symbol %in% c("KIAA0368", "FANCC", "ZFX", "XIST",  "SYAP1", 
  "KDM6A", "KDM5C", "JPX", "EIF1AX", "DDX3X")) %>%
  plot_gene_concordance(main = "Inattention index")

png("figures/genes_pheno/heatmap_innatention_tissues_selGenes.png", width = 1500, height = 1500, res = 300)
innatention_genes
dev.off()

png("figures/genes_pheno/heatmap_hyper_tissues.png", width = 1500, height = 1800, res = 300)
phenos_tabs$hs_Hyper_raw %>%
  filter(!Symbol %in% c("NA", "", "TXLNG;AX746622", "LOC389906", "HDHD1")) %>%
  plot_gene_concordance(main = "Hyperactivity index")
dev.off()


hyper_genes <- phenos_tabs$hs_Hyper_raw %>%
  filter(Symbol %in% c("TTL", "ZRSR2", "XIST", "UBA1", "SYAP1", 
  "KDM6A", "KDM5C", "JPX", "EIF1AX", "DDX3X", "CA5BP1", "CA5B")) %>%
  plot_gene_concordance(main = "Hyperactivity index")

png("figures/genes_pheno/heatmap_hyper_tissues_selGenes.png", width = 1500, height = 1500, res = 300)
hyper_genes
dev.off()


png("figures/genes_pheno/heatmap_adhd_second_tissues_selGenes.png", width = 3000, height = 1500, res = 300)
plot_grid(innatention_genes + theme(legend.position = "none"), hyper_genes, ncol = 2, labels = c("A", "B"))
dev.off()

png("figures/genes_pheno/heatmap_panel.png", width = 3000, height = 2000, res = 300)
plot_grid(weight_genes, skf_genes, adhd_genes, ncol = 3, 
labels = "AUTO", rel_widths = c(1, 1, 1.4))
dev.off()

## Check correlation with GTEx phenotypes
adipose <- gtex_tissues_sva[["Adipose - Visceral (Omentum)"]]

gtex_phenotype <- read_delim("data/phenosGTex.csv")

pheno_adipose <- left_join(colData(adipose) %>% data.frame() %>% 
  dplyr::select(gtex.subjid, gtex.smtsisch, gtex.smrin, starts_with("SV")), 
  gtex_phenotype, by = c("gtex.subjid" = "SUBJID"))

pheno_adipose$BMI <- as.numeric(pheno_adipose$BMI)
pheno_adipose$WGHT <- as.numeric(pheno_adipose$WGHT)

sv_cols_adip <- grep("^SV", colnames(pheno_adipose), value = TRUE)
mod_adipose <- model.matrix(formula(paste("~ BMI + AGE + gtex.smtsisch + gtex.smrin +", paste(sv_cols_adip, collapse = "+"))),
                                data = pheno_adipose)
adipose_noY <- adipose[as.character(seqnames(adipose)) %in% paste0("chr", c(1:22, "X")), ]

## Run Voom and limma
adipose_voom <- voom(assay(adipose_noY), mod_adipose)
lm.adipose <- lmFit(adipose_voom, mod_adipose) %>%
  eBayes()
tab.adipose <- topTable(lm.adipose, coef = 2, n = Inf, confint = TRUE)
tab.adipose$Symbol <- rowData(adipose_noY)[rownames(tab.adipose), "gene_name"]


top_bmi <- subset(tab.adipose, Symbol %in%  c("ZC3H12D", "TRAK2", "TCP11L2", 
  "PSME4", "PGM2L1", "PCBP2", "LEPROTL1", "IFIT1B", "ESPN", "EBF1", "DDX18", "CXCL5", "CTSE",
  "AP2A1", "UBA1", "STS", "MED14", "EIF2S3", "EIF1AX", "EFHC2", "CA5BP1", "ALAS2" ))


mod_weight <- model.matrix(formula(paste("~ WGHT + AGE + gtex.smtsisch + gtex.smrin +", paste(sv_cols_adip, collapse = "+"))),
                                data = pheno_adipose)
## Run Voom and limma
weight_voom <- voom(assay(adipose_noY[, !is.na(pheno_adipose$WGHT)]), mod_weight)
lm.weight <- lmFit(weight_voom, mod_weight) %>%
  eBayes()
tab.weight <- topTable(lm.weight, coef = 2, n = Inf, confint = TRUE)
tab.weight$Symbol <- rowData(adipose_noY)[rownames(tab.weight), "gene_name"]


top_weight <- subset(tab.weight, Symbol %in%  c("ZC3H12D", "TRAK2", "TCP11L2", 
  "PSME4", "PGM2L1", "PCBP2", "LEPROTL1", "IFIT1B", "ESPN", "EBF1", "DDX18", "CXCL5", "CTSE",
  "AP2A1", "UBA1", "STS", "MED14", "EIF2S3", "EIF1AX", "EFHC2", "CA5BP1", "ALAS2" ))

pheno_adipose2 <- pheno_adipose[!is.na(pheno_adipose$WGHT), ]

library(edgeR)
adipose_cpm <- cpm(assay(adipose_noY), log = TRUE, prior.count = 1)
adipose_cpm <- adipose_cpm[, !is.na(pheno_adipose$WGHT)]

df_gene <- tibble(ebf1 = adipose_cpm["ENSG00000164330.16", ],
  sex = ifelse(pheno_adipose2$SEX == 1, "Men", "Women"),
  ifit1b = adipose_cpm["ENSG00000204010.3", ],
  cxcl5 = adipose_cpm["ENSG00000163735.6", ],
  sts = adipose_cpm["ENSG00000101846.6", ],
  ddx18 = adipose_cpm["ENSG00000088205.12", ],
  pcbp2 = adipose_cpm["ENSG00000197111.15", ],
  ap2a1 = adipose_cpm["ENSG00000196961.12", ],
  age = pheno_adipose2$AGE, bmi = pheno_adipose2$BMI, weight = pheno_adipose2$WGHT, 
  ischemic = pheno_adipose2$gtex.smtsisch, rin = pheno_adipose2$gtex.smrin) %>%
  cbind(pheno_adipose2 %>% dplyr::select(starts_with("SV")) %>% data.frame())



summary(lm(bmi ~ sex + age , df_gene))
summary(lm(weight ~ sex + age , df_gene))


summary(lm(bmi ~ sex + age + ebf1, df_gene))
summary(lm(weight ~ sex + age + ebf1, df_gene))

summary(lm(bmi ~ sex + age + ifit1b, df_gene))
summary(lm(weight ~ sex + age + ifit1b, df_gene))

summary(lm(bmi ~ sex + age + cxcl5, df_gene))
summary(lm(weight ~ sex + age + cxcl5, df_gene))

summary(lm(bmi ~ sex + age + ddx18, df_gene))
summary(lm(weight ~ sex + age + ddx18, df_gene))
summary(lm(weight ~ sex*ddx18 + age, df_gene))
summary(lm(ddx18 ~ sex*age, df_gene))

summary(lm(bmi ~ sex + age + sts, df_gene))
summary(lm(weight ~ sex + age + sts , df_gene))

summary(lm(bmi ~ sex + age + pcbp2, df_gene))
summary(lm(weight ~ sex + age + pcbp2 , df_gene))

summary(lm(bmi ~ sex + age + ap2a1, df_gene))
summary(lm(weight ~ sex + age + ap2a1 , df_gene))
summary(lm(weight ~ sex*ap2a1 + age, df_gene))

cor(df_gene$ebf1, df_gene$weight)
cor(df_gene$ifit1b, df_gene$weight)
cor(df_gene$cxcl5, df_gene$weight)
cor(df_gene$sts, df_gene$weight)
cor(df_gene$ddx18, df_gene$weight)
cor(df_gene$ap2a1, df_gene$weight)

cor(df_gene$ebf1, df_gene$bmi)
cor(df_gene$ifit1b, df_gene$bmi)
cor(df_gene$cxcl5, df_gene$bmi)
cor(df_gene$sts, df_gene$bmi)
cor(df_gene$ddx18, df_gene$bmi)
cor(df_gene$ap2a1, df_gene$bmi)

weight_adipose <- ggplot(df_gene[df_gene$ddx18 > 4.5,], aes(x = ddx18, y = weight, color = sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(~sex) +
  scale_color_manual(values = c("Men" = "blue", "Women" = "maroon")) +
  ylab("Weight (lb)") +
  xlab("DDX18 expression") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

png("figures/genes_pheno/weight_ddx18_adipose.png", width = 1500, height = 1500, res = 300)
weight_adipose
dev.off()

weight_adipose2 <- ggplot(df_gene, aes(x = pcbp2, y = weight, color = sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(~sex) +
  scale_color_manual(values = c("Men" = "blue", "Women" = "maroon")) +
  ylab("Weight (lb)") +
  xlab("PCBP2 expression") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

png("figures/genes_pheno/weight_pcbp2_adipose.png", width = 1500, height = 1500, res = 300)
weight_adipose2
dev.off()



sex_adipose <- ggplot(df_gene[df_gene$ddx18 > 4.5,], aes(y = ddx18, x = sex, color = sex)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = c("Men" = "blue", "Women" = "maroon")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggtitle("Adipose Visceral (Omentum)") +
  ylab("DDX18 expression") +
  xlab("Sex") 

png("figures/genes_pheno/sex_ddx18_adipose.png", width = 1500, height = 1500, res = 300)
sex_adipose
dev.off()

sex_adipose2 <- ggplot(df_gene, aes(y = ap2a1, x = sex, color = sex)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = c("Men" = "blue", "Women" = "maroon")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggtitle("Adipose Visceral (Omentum)") +
  ylab("AP2A1 expression") +
  xlab("Sex") 

png("figures/genes_pheno/sex_ap2a1_adipose.png", width = 1500, height = 1500, res = 300)
sex_adipose2
dev.off()

png("figures/genes_pheno/assocs_panel.png", width = 3000, height = 3000, res = 300)
plot_grid(
  plot_grid(weight_genes, 
    plot_grid(sex_adipose, weight_adipose, ncol = 1, labels = c("B", "C")),
  ncol = 2, labels = c("A", "")),
  plot_grid(skf_genes, adhd_genes, ncol = 2, labels = c("D", "E")),
  nrow = 2
)
dev.off()

# runMediations <- function(pheno, tc){
#   pheno <- as.character(pheno)
#   tc <- as.character(tc)
#   df <- data.frame(pheno = all_phenos[, pheno],
#                      tc = as.numeric(assay(se.noY[tc, ])))
#   df <- cbind(df, colData(se.noY)) %>%
#     filter(!is.na(pheno))

#   mod.med <- lm(formula(paste("tc ~ e3_sex + cohort +  age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6 +",
#   paste(paste0("SV", 1:19), collapse = "+"))), df)
    
#   if (pheno %in% c("hs_ADHD_raw", "hs_Hyper_raw", "hs_Cognit_raw")){
#      mod.out <-  glm.nb(pheno ~ tc + e3_sex + cohort + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, 
#               df)
#   } else {
#     mod.out <- lm(pheno ~ tc + e3_sex + cohort + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, df)
#   }
#   med <- mediate(mod.med, mod.out, treat = "e3_sex", mediator = "tc", sims = 800)
#   med
# }

# pheno_gene.tab2$coef_tc_sex <- tab.gexp[pheno_gene.tab2$TC, "logFC"]
# pheno_gene.tab2 <- left_join(pheno_gene.tab2, select(df_pheno_assocs, var, coef), 
#   by = c("pheno" = "var"))
# pheno_gene.tab2 <- mutate(pheno_gene.tab2, 
#                          coherent = ifelse(sign(coef_tc_sex) * sign(coef_pheno_tc) == sign(coef ), "Coherent", "Discordant")
# )

# pheno_gene.coherent2 <- subset(pheno_gene.tab2, coherent == "Coherent")
# mediation.list2 <- mclapply(seq_len(nrow(pheno_gene.coherent2)), function(i) {
#   runMediations(pheno = pheno_gene.coherent2[i, ]$pheno, tc = pheno_gene.coherent2[i, ]$TC)
# }, mc.cores = 15)


# mediation.df2 <- as_tibble(pheno_gene.coherent2) %>%
#   mutate(med.prop = sapply(mediation.list2, function(x) x$n.avg),
#          med.pvalue = sapply(mediation.list2, function(x) x$n.avg.p),
#          med.fdr = p.adjust(med.pvalue, "fdr"))
# mediation.df2$Symbol <- rowData(se.noY[mediation.df2$TC, ])$GeneSymbolDB2
# save(mediation.list2, mediation.df2, file = "results/mediation/mediation_phenos_metachild.Rdata")


# n_gene <- table(mediation.df2$TC)
# subset(mediation.df2, TC %in% names(n_gene[n_gene > 1]))  %>%
#   arrange(TC) %>% data.frame() %>% select(TC, pheno, Symbol)

# ## Sup Table 4
# mediation.df2_out <- mutate(mediation.df2, 
#   Phenotype = pheno_vec_name[pheno],
#   coef_pheno_sex = coef, 
#   med_prop = med.prop,
#   med_pvalue = med.pvalue,
#   med_fdr = med.fdr) %>%
#   select(TC, Symbol, Phenotype, starts_with("coef"), starts_with("med"))
# write.table(mediation.df2_out, file = "results/mediation/mediation_coding.txt", 
#   quote = FALSE, row.names = FALSE, sep = "\t")


# mediation_sig2 <- subset(mediation.df2, med.fdr < 0.05)

# write.table(mediation_sig2, file = "results/mediation/mediation_coding_sig.txt", 
#   quote = FALSE, row.names = FALSE)

# subset(child_comb, Symbol %in% mediation_sig2$Symbol) %>% data.frame() %>% 
#   select(Symbol, starts_with(c("logFC", "P.Value")))

# subset(tab.adult, Symbol %in% mediation_sig2$Symbol) %>% data.frame() %>% 
#   select(Symbol, starts_with(c("logFC", "P.Value")))

# comb_pheno <- t(assay(se.noY[mediation_sig2$TC,])) %>%  
#   as_tibble() %>%
#   mutate(HelixID = colnames(se.noY)) %>%
#   inner_join(all_phenos, by = "HelixID")

# raw_mod <- summary(lm(hs_skf_sum2 ~ e3_sex + cohort + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, comb_pheno))
# med_mod <- summary(lm(hs_skf_sum2 ~ TC0X000980.hg.1 + TC04001283.hg.1 + TC06004128.hg.1  + TC0X001077.hg.1 + e3_sex + cohort + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, comb_pheno))
# 1 - med_mod$coefficients[6, 1]/ raw_mod$coefficients[2, 1]

# ### Gene descriptives ####
# mediator_genes <- unique(mediation_sig2$TC)

# mediator_genes_df <- data.frame(TC = mediator_genes) %>%
#   mutate(chr = as.character(seqnames(se.noY[TC, ])),
#          Coding = ifelse(TC %in% codingTCs, "Coding", "non-Coding"), 
#          Symbol = rowData(se.noY)[TC, ]$GeneSymbolDB2, 
#          Symbol = ifelse(Symbol == "NA" | Symbol == "", TC, Symbol)) %>% 
#   left_join(xci_genes[, c("Symbol", "Type", "XCI")],
#             by = "Symbol")


# ## N genes
# length(mediator_genes)

# ## N genes autosomic
# table(mediator_genes_df$chr == "chrX")


# ### XCI
# table(mediator_genes_df$XCI)


# ## Coding
# table(mediator_genes_df$Coding)

# ### Plot mediations ####
# ## Sup Figures 4-7 
# plot_mediations <- function(pheno,pheno_name,  TC, gene_name){

#   comb_pheno <- tibble(tc = as.numeric(assay(se.noY[TC,])),
#                       HelixID = colnames(se.noY)) %>%
#     inner_join(all_phenos, by = "HelixID") %>%
#     mutate(Sex = ifelse(e3_sex == "male", "Boys", "Girls"))
#   comb_pheno$Pheno <- comb_pheno[[pheno]]

#   plot_sex_tc <- ggplot(comb_pheno, aes(x = Sex, y = tc)) +
#     geom_boxplot() +
#     xlab("Sex") +
#     geom_label(
#       aes(x = Inf, y = Inf),  
#       label = sprintf("P = %.1e", tab.gexp[TC, ]$P.Value), 
#       hjust = 1, vjust = 1
#     ) +
#     ylab(paste(gene_name, "Expression")) +  
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5))


#   plot_sex_pheno <- ggplot(comb_pheno, aes(x = Sex, y = Pheno)) +
#     geom_boxplot() +
#     xlab("Sex") +
#     ylab(pheno_name) +
#     geom_label(
#       aes(x = Inf, y = Inf),  
#       label = sprintf("P = %.1e", (subset(df_pheno_assocs, var == pheno)$pval)), 
#       hjust = 1, vjust = 1
#     ) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5))

#   plot_tc_pheno <- ggplot(comb_pheno, aes(x = tc, y = Pheno)) +
#     geom_point() +
#     geom_smooth(method = "lm") +
#     xlab(paste(gene_name, "Expression")) +
#     ylab(pheno_name) +
#     geom_label(
#       aes(x = Inf, y = Inf),  
#       label = sprintf("P = %.1e", pheno_assocs.dim_cod[[pheno]][TC, ]$P.Value), 
#       hjust = 1, vjust = 1
#     ) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5))


#     title <- ggdraw() + 
#     draw_label(paste("Sex -", gene_name, "-", pheno_name),
#                fontface = 'bold',  hjust = 0.5  ) 
#   plot_grid(title, 
#     plot_grid(plot_sex_pheno, plot_sex_tc, plot_tc_pheno, nrow = 1, labels = LETTERS[1:3]),
#   ncol = 1, rel_heights = c(1, 10))
# }

# mediation_plots <- lapply(seq_len(nrow(mediation_sig2)), function(i){

#   pheno <- mediation_sig2$pheno[i]
#   pheno_name <- pheno_vec_name[pheno]
#   tc <- mediation_sig2$TC[i]
#   symbol <- mediation_sig2$Symbol[i]

#   plot_mediations(pheno, pheno_name, tc, symbol)
# })


# png("figures/mediation/efhc2_sk.png", height = 900, width = 3000, res = 300)
# mediation_plots[[1]]
# dev.off()

# png("figures/mediation/cxcl5_sk.png", height = 900, width = 3000, res = 300)
# mediation_plots[[2]]
# dev.off()

# png("figures/mediation/hla_sk.png", height = 900, width = 3000, res = 300)
# mediation_plots[[3]]
# dev.off()

# png("figures/mediation/alas2_sk.png", height = 900, width = 3000, res = 300)
# mediation_plots[[4]]
# dev.off()
