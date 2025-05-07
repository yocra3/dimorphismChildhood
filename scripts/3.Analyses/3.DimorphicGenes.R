#'##############################################################################
#' Dimorphic genes 
#' docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.8 R
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(cowplot)
library(SummarizedExperiment)
library(UpSetR)
library(S4Vectors)
library(limma)
library(readxl)
library(ggmanh)
library(ggrepel)
library(ggbreak)
library(enrichR)
library(meta)
library(eulerr)
library(ggpubr)

## Load and preprocess gene expression ####
load("./results/preprocessFiles/Expression_SE_raw_filtered_SV.RData")
load("data/gexpAnnotation.Rdata")

se.filt$cohort <- droplevels(se.filt$cohort)

xci_genes <- read_xlsx("data/XCI_genes.PMID29022598.xlsx", skip = 1) 
xci_genes <- xci_genes[, 1:7]
colnames(xci_genes) <- c("Symbol", "Ensembl", "Chr", "Start", "End", "Type", "XCI")

## Dimorphic TCs ####
### Remove genes in chrY and genes in pseudo autosomal regions
se.noY <- se.filt[seqnames(se.filt) != "chrY", ]
# pseudo_reg <- GRanges(c("chrX:60001-2699520", "chrX:154931044-155260560")) ## hg19
# pseudo_genes <- subsetByOverlaps(rowRanges(se.noY), pseudo_reg)
# se.noY <- se.noY[!rownames(se.noY) %in% names(pseudo_genes), ]
dim(se.noY)

#[1] 30013   734
mod.helix <- model.matrix(formula(paste("~ e3_sex + cohort + age_sample_years + NK_6 +
                                        Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6 +", 
                                              paste(paste0("SV", 1:19), collapse = "+"))),
                            colData(se.noY ))

lm.gexp <- lmFit(t(scale(t(assay(se.noY)))), mod.helix) %>%
  eBayes()
tab.gexp <- topTable(lm.gexp, coef = 2, n = Inf, confint = TRUE)
tab.gexp$chr <- as.character(seqnames(se.noY[rownames(tab.gexp), ]))
tab.gexp$Symbol <-  rowData(se.noY)[rownames(tab.gexp), ]$GeneSymbolDB2

### Coding TCs ####
codingTCs <- subset(expAnnot, Coding == "coding" & locus.type == "Coding" & transcript_cluster_id %in% rownames(se.noY))$transcript_cluster_id
tab.gexp$Coding <- ifelse(rownames(tab.gexp) %in% codingTCs, "Coding", "Non-coding")

save(tab.gexp, file = "results/genes/HELIX_dimorphic_genes.Rdata")

## Sup File 1
write.table(tab.gexp %>%
  mutate(TC = rownames(.)) %>%
  select(TC, Symbol, chr, logFC, CI.L, CI.R, P.Value, adj.P.Val), file = "results/genes/helix_analysis.txt", sep = "\t", 
  quote = FALSE, row.names = FALSE)


## Comparison with children ####
### Controls of T1D study
load("results/preprocessFiles/gexp_GSE43488_sex.Rdata")
### Run model ####
mod.gse43488 <- model.matrix(formula(paste("~ sex + age +", 
                                              paste(paste0("SV", 1:2), collapse = "+"))),
                            colData(gse43488_controls ))
lm.gse43488 <- lmFit(t(scale(t(assay(gse43488_controls)))), mod.gse43488) %>%
  eBayes()
tab.gse43488 <- topTable(lm.gse43488, coef = 2, n = Inf, confint = TRUE)
tab.gse43488$Symbol <- rownames(tab.gse43488)
save(tab.gse43488, file = "results/genes/gse43488_analysis.Rdata")

## Comparison with Generation R ####
load("results/GenR/gtex_analysis.Rdata")
tab.genR <- tab.gtex
rm(tab.gtex)

## Sup File 2
write.table(tab.genR %>%
  mutate(Chromosome = chromosome) %>%
  select(Symbol, Chromosome, logFC, CI.L, CI.R, P.Value, adj.P.Val), 
  file = "results/genes/GenR_analysis.txt", sep = "\t", 
  quote = FALSE, row.names = FALSE)



### Combine with HELIX ####
child_comb <- tab.gexp %>% 
  as_tibble() %>%
  mutate(TC = rownames(tab.gexp)) %>%
  filter(!is.na(Symbol) & Symbol != "" & Symbol != "NA") %>%
  dplyr::select(Symbol, TC, logFC, P.Value, adj.P.Val,chr, CI.L, CI.R, chr) %>%
  left_join(dplyr::select(tab.gse43488, Symbol, logFC,  P.Value, adj.P.Val, CI.L, CI.R) %>%
              group_by(Symbol) %>% filter(P.Value == min(P.Value)), by = "Symbol",
            suffix = c(".helix", ".gse43488")) %>%
  left_join(dplyr::select(tab.genR, Symbol, logFC,  P.Value, adj.P.Val, CI.L, CI.R) %>%
              group_by(Symbol) %>% filter(P.Value == min(P.Value)), by = "Symbol",
            suffix = c("", ".GenR")) %>%
  filter(!is.na(P.Value.gse43488) & !is.na(P.Value))

com_genes_child <- nrow(child_comb)

#### Run meta-analysis in common genes ####
meta_list <- lapply(seq_len(nrow(child_comb)), function(i) {
  metagen(TE = c(child_comb$logFC.helix[i], child_comb$logFC.gse43488[i], child_comb$logFC[i]), 
          upper = c(child_comb$CI.R.helix[i], child_comb$CI.R.gse43488[i], child_comb$CI.R[i]), 
          lower = c(child_comb$CI.L.helix[i], child_comb$CI.L.gse43488[i], child_comb$CI.L[i]))
})
child_comb$logFC.com <- sapply(meta_list, function(x) ifelse(x$pval.Q < 0.05/com_genes_child, x$TE.random,x$TE.common))
child_comb$P.Value.com <- sapply(meta_list, function(x) ifelse(x$pval.Q < 0.05/com_genes_child, x$pval.random, x$pval.common))
child_comb$adj.P.Val.com <- p.adjust(child_comb$P.Value.com, method = "BH")
save(child_comb, file = "results/genes/children_analysis.Rdata")


## Sup File 3
write.table(child_comb %>%
  mutate(logFC.T1D = logFC.gse43488,
          P.Value.T1D = P.Value.gse43488,
          adj.P.Val.T1D = adj.P.Val.gse43488,
          CI.L.T1D = CI.L.gse43488,
          CI.R.T1D = CI.R.gse43488,
          logFC.GenR = logFC,
          P.Value.GenR = P.Value,
          adj.P.Val.GenR = adj.P.Val,
          CI.L.GenR = CI.L,
          CI.R.GenR = CI.R,
          logFC.Meta = logFC.com,
          P.Value.Meta = P.Value.com,
          adj.P.Val.Meta = adj.P.Val.com) %>%
  select(Symbol, chr, TC, ends_with("helix"), ends_with("T1D"), ends_with("GenR"), ends_with("Meta")), file = "results/genes/children_metaanalysis.txt", sep = "\t", 
  quote = FALSE, row.names = FALSE)

child_comb_sig <- subset(child_comb, adj.P.Val.com < 0.05)
dim.gexp <- subset(child_comb_sig, adj.P.Val < 0.05)$TC


### All TCs ####
## N genes
sum(child_comb$adj.P.Val.com < 0.05)
mean(child_comb$adj.P.Val.com < 0.05)

## X Chromosome prop
xchr_tab <- table(sig = child_comb$adj.P.Val.com < 0.05, chr = child_comb$chr == "chrX")
prop.table(xchr_tab, margin = 1)
prop.table(xchr_tab, margin = 2)

fisher.test(xchr_tab)

### XCI
child_comb_chrX <- subset(child_comb, chr == "chrX")
child_comb_chrX <- left_join(child_comb_chrX, xci_genes[, c("Symbol", "Type", "XCI")],
                            by = "Symbol")
child_comb_chrX$XCI2 <- ifelse(child_comb_chrX$XCI %in% c("inactive", "variable"), "Inactive", child_comb_chrX$XCI)

xci_tab <- table(child_comb_chrX$XCI, sig = child_comb_chrX$adj.P.Val.com < 0.05)
chisq.test(xci_tab)
sum(!is.na(child_comb_chrX$XCI) & child_comb_chrX$adj.P.Val.com < 0.05)

fisher.test(xci_tab[c(2, 1), ])

child_comb_chrX_sig <- subset(child_comb_chrX, adj.P.Val.com < 0.05)
mean(child_comb_chrX_sig$XCI == "escape", na.rm = TRUE)


# Proportion of genes DE by sex
prop.table(table(ifelse(child_comb_sig$logFC.com < 0, "Girls", "Boys"), 
  ifelse(child_comb_sig$chr == "chrX", "chrX", "Autosomal")), margin = 2)

prop.table(
  table(ifelse(child_comb_chrX_sig$logFC.com < 0, "Girls", "Boys"), 
   child_comb_chrX_sig$XCI2)
, margin = 2)
child_comb_sig %>%
  arrange(P.Value.com) %>%
  head(40) %>%
  data.frame() %>%
  select(Symbol, chr,logFC.com, adj.P.Val.com )

## Enrichments ####
dbs <- listEnrichrDbs()

db.vec <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "DisGeNET")


enrich_genes_comb <- child_comb_sig$Symbol
background_child <-  child_comb$Symbol
enrich_res_child <- enrichr(enrich_genes_comb, db.vec, background = background_child,
  include_overlap = TRUE)
names(enrich_res_child) <- db.vec
enrich_child_df <-  Reduce(rbind,enrich_res_child) %>%
  mutate(DB = rep(names(enrich_res_comb), sapply(enrich_res_comb, nrow))) %>%
  filter(Adjusted.P.value < 0.05) 


enrich_genesX_comb <- subset(child_comb_sig, chr == "chrX")$Symbol
background_childX <-  subset(child_comb, chr == "chrX")$Symbol
enrich_resX_comb <- enrichr(enrich_genesX_comb, db.vec,
  background = background_childX,
  include_overlap = TRUE)
enrich_combX_df <-  Reduce(rbind,enrich_resX_comb) %>%
  mutate(DB = rep(names(enrich_resX_comb), sapply(enrich_resX_comb, nrow))) %>%
  filter(Adjusted.P.value < 0.05) 

enrich_genes.children <- enrich_comb_df %>% 
  mutate(Type = ifelse(DB == "DisGeNET", "Phenotypes and Diseases", "Transcription Factor"),
          Overlap = Overlap, 
          OR = Odds.Ratio, 
          pvalue = P.value, 
          adj_pvalue = Adjusted.P.value) %>%
  select(Term, DB, Type, Overlap, OR, pvalue, adj_pvalue)
write.table(enrich_genes.comb_children, file = "results/genes/dimorphic_gene_enrichr_children.txt",
            sep = "\t", quote = FALSE, row.names = FALSE) 

child_TFs <- subset(enrich_comb_df, DB == "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")$Term  
child_TFs <- unique(gsub(" ENCODE", "", gsub(" CHEA", "", child_TFs)))

subset(child_comb, Symbol %in% child_TFs) %>% 
  dplyr::select(Symbol, adj.P.Val.com, starts_with(c("logFC", "P.Value"))) %>%
  arrange(P.Value.com) %>%
  data.frame()

subset(child_comb, Symbol %in% child_TFs) %>% 
  dplyr::select(Symbol, adj.P.Val.com, adj.P.Val.helix, starts_with(c("logFC", "P.Value"))) %>%
  arrange(P.Value.com) %>%
  data.frame() %>%
  filter(adj.P.Val.com < 0.05)


### Manhattan ####
## Figure 3A
man_tab_child <- child_comb
man_tab_child$Position <- start(rowRanges(se.filt[man_tab_child$TC, ]))
man_tab_child$chr <- factor(man_tab_child$chr, levels = paste0("chr", c(1:22, "X")))
man_tab_child$chrNum <- factor(gsub("chr", "", man_tab_child$chr), levels = c(1:22, "X"))
man_tab_child$Genes <- man_tab_child$Symbol
man_tab_child$Genes <- ifelse(man_tab_child$P.Value.com < 2e-26, man_tab_child$Genes, "")
man_tab_child$Genes[man_tab_child$Genes == "NA"] <- ""
man_tab_child$pval <- ifelse(man_tab_child$P.Value.com < 1e-200, 1e-200, man_tab_child$P.Value.com)


man_plot_children <- manhattan_plot(x = man_tab_child, pval.colname = "pval", 
                           signif = c(0.05/nrow(man_tab_child)),
                           chr.colname = "chrNum", pos.colname = "Position",
                           plot.title = "Children (1,009 individuals)",
                           rescale = FALSE,
                           label.colname = "Genes") +
  theme(plot.title = element_text(hjust = 0.5)) 


man_plot_children$layers[[2]] <- NULL ## Remove significance line
png("figures/manhattan_children.png", res = 300, height = 1200, width = 2000)
man_plot_children
dev.off()

### Volcano ####
### Figure 3B
volcano_plot_children <- man_tab_child %>%
  mutate(Signif = ifelse(adj.P.Val.com < 0.05, ifelse(logFC.com > 0, "Boys", "Girls"), "None")) %>%
  ggplot(aes(x = logFC.com, y = -log10(pval), color = Signif, label = Genes)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_label_repel(max.overlaps = 15) +
  scale_color_manual(values = c("blue", "maroon", "black")) +
  theme_bw() +
  xlab("Normalized Effect size") +
  ylab(expression(-log[10]("p"))) + 
  annotate("text", x = -2, y = -30, label = "Higher Girls") + 
  annotate("text", x = 0.5, y = -30, label = "Higher Boys") + 
  coord_cartesian(ylim = c(0, 230), clip="off") +
  ggtitle("Children") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

png("figures/volcano_children.png", res = 300, height = 1200, width = 2000)
volcano_plot_children
dev.off()


# Comparison with adults ####
## GSE36382 ####
load("results/preprocessFiles/gexp_gse36382_sex.Rdata")

### Run model ####
### Select probes mapped to chr1-22 and X and probes with symbol
gse36382_noY <- gse36382_se[as.character(seqnames(gse36382_se)) != "chrY", ]
mod_gse36382 <- model.matrix(~ sex, data = colData(gse36382_se))

assay(gse36382_noY)$E <- t(scale(t(assay(gse36382_noY)$E)))
lm.gse36382 <- lmFit(assay(gse36382_noY), mod_gse36382) %>%
  eBayes()
tab.gse36382 <- topTable(lm.gse36382, coef = 2, n = Inf, confint = TRUE)
tab.gse36382$Symbol <- rownames(tab.gse36382)
tab.gse36382$chromosome <- as.character(seqnames(gse36382_noY[rownames(tab.gse36382), ]))
save(tab.gse36382, file = "results/genes/gse36382_analysis.Rdata")

## GSE33828 ####
load("results/preprocessFiles/gexp_gse33828_sex.Rdata")
### Run model ####
gse33828_noY <- gse33828_se[as.character(seqnames(gse33828_se)) != "chrY", ]
mod_gse33828 <- model.matrix(~ sex + age, data = colData(gse33828_noY))

lm.gse33828 <- lmFit(t(scale(t(assay(gse33828_noY)))), mod_gse33828) %>%
  eBayes()
tab.gse33828 <- topTable(lm.gse33828, coef = 2, n = Inf, confint = TRUE)
tab.gse33828$Symbol <- rownames(tab.gse33828)
tab.gse33828$chromosome <- as.character(seqnames(gse33828_noY[rownames(tab.gse33828), ]))
save(tab.gse33828, file = "results/genes/gse33828_analysis.Rdata")

## GTEx ####
load("results/preprocessFiles/gexp_gtex_sva_sex.Rdata")
mod_gtex <- model.matrix(formula(paste("~ sex + gtex.age +", 
                                              paste(paste0("SV", 1:43), collapse = "+"))),
                                data = colData(blood_gtex_filt))

blood_gtex_noY <- blood_gtex_filt[as.character(seqnames(blood_gtex_filt)) %in% paste0("chr", c(1:22, "X")), ]

gtex_voom <- voom(assay(blood_gtex_noY), mod_gtex)
gtex_voom$E <- t(scale(t(gtex_voom$E)))

lm.gtex <- lmFit(gtex_voom, mod_gtex) %>%
  eBayes()
tab.gtex <- topTable(lm.gtex, coef = 2, n = Inf, confint = TRUE)
tab.gtex$Symbol <- rowData(blood_gtex_noY)[rownames(tab.gtex), "gene_name"]
tab.gtex$chromosome <- as.character(seqnames(blood_gtex_noY[rownames(tab.gtex)]))
save(tab.gtex, file = "results/genes/gtex_analysis.Rdata")

## Combine adult results ####
### Sup Table 4
tab.adult <- inner_join(tab.gse33828, tab.gse36382, by = "Symbol", 
                        suffix = c(".gse33828", ".gse36382")) %>%
  inner_join(tab.gtex, by = "Symbol", suffix = c(".array", ".gtex")) %>%
  as_tibble()

com_genes_adult <- nrow(tab.adult)

meta_list_adult <- lapply(seq_len(nrow(tab.adult)), function(i) {
  metagen(TE = c(tab.adult$logFC.gse33828[i], tab.adult$logFC.gse36382[i], tab.adult$logFC[i]), 
          upper = c(tab.adult$CI.R.gse33828[i], tab.adult$CI.R.gse36382[i], tab.adult$CI.R[i]), 
          lower = c(tab.adult$CI.L.gse33828[i], tab.adult$CI.L.gse36382[i], tab.adult$CI.L[i]))
})
tab.adult$logFC.com <- sapply(meta_list_adult, function(x) ifelse(x$pval.Q < 0.05/com_genes_adult, x$TE.random,x$TE.common))
tab.adult$P.Value.com <- sapply(meta_list_adult, function(x) ifelse(x$pval.Q < 0.05/com_genes_adult, x$pval.random, x$pval.common))
tab.adult$adj.P.Val.com <- p.adjust(tab.adult$P.Value.com, method = "BH")
save(tab.adult, file = "results/genes/adult_analysis.Rdata")



write.table(tab.adult %>%
  mutate(logFC.RS3 = logFC.gse33828,
          P.Value.RS3 = P.Value.gse33828,
          adj.P.Val.RS3 = adj.P.Val.gse33828,
          CI.L.RS3 = CI.L.gse33828,
          CI.R.RS3 = CI.R.gse33828,
          logFC.ship = logFC.gse36382,
          P.Value.ship = P.Value.gse36382,
          adj.P.Val.ship = adj.P.Val.gse36382,
          CI.L.ship = CI.L.gse36382,
          CI.R.ship = CI.R.gse36382,
          logFC.gtex = logFC,
          P.Value.gtex = P.Value,
          adj.P.Val.gtex = adj.P.Val,
          CI.L.gtex = CI.L,
          CI.R.gtex = CI.R,
          logFC.Meta = logFC.com,
          P.Value.Meta = P.Value.com,
          adj.P.Val.Meta = adj.P.Val.com,
          Chromosome = chromosome) %>%
  select(Symbol, Chromosome,  ends_with("RS3"), ends_with("ship"), ends_with("gtex"), ends_with("Meta")), 
  file = "results/genes/adult_metaanalysis.txt", sep = "\t", 
  quote = FALSE, row.names = FALSE)


### Make plots #####
tab.adult$chr <- factor(tab.adult$chromosome, levels = paste0("chr", c(1:22, "X")))
tab.adult$chrNum <- factor(gsub("chr", "", tab.adult$chr), levels = c(1:22, "X"))
adult_genes <- subset(rowRanges(blood_gtex_filt), gene_name %in% tab.adult$Symbol)
adult_genes <- adult_genes[!duplicated(adult_genes$gene_name)]
names(adult_genes) <- adult_genes$gene_name
tab.adult$Position <- start(adult_genes[tab.adult$Symbol])
tab.adult$pval <- ifelse(tab.adult$P.Value.com < 1e-200, 1e-200, tab.adult$P.Value.com)
tab.adult$Genes <- tab.adult$Symbol
tab.adult$Genes <- ifelse(tab.adult$pval < 5e-20, tab.adult$Genes, "")
tab.adult$Genes[tab.adult$Genes == "NA"] <- ""

## Figure 3C
man_adult_plot <- manhattan_plot(x = tab.adult, pval.colname = "pval", 
                           signif = 0.05/nrow(tab.adult),
                           chr.colname = "chrNum", pos.colname = "Position",
                           plot.title = "Adults (2,723 individuals)",
                           rescale = FALSE,
                           label.colname = "Genes") +
  theme(plot.title = element_text(hjust = 0.5))


man_adult_plot$layers[[2]] <- NULL ## Remove significance line

png("figures/manhattan_adult.png", res = 300, height = 1200, width = 2000)
man_adult_plot
dev.off()

### Figure 3D
volcano_adult_plot <- tab.adult %>%
  mutate(Signif = ifelse(adj.P.Val.com < 0.05, ifelse(logFC > 0, "Boys", "Girls"), "None")) %>%
  ggplot(aes(x = logFC, y = -log10(pval), color = Signif, label = Genes)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_label_repel() +
  scale_color_manual(values = c("blue", "maroon", "black")) +
  theme_bw() +
  xlab("Normalized Effect size") +
  ggtitle("Adults") +
  ylab(expression(-log[10]("p"))) + 
  annotate("text", x = -2, y = -45, label = "Higher Women") + 
  annotate("text", x = 0.5, y = -45, label = "Higher Men") + 
  coord_cartesian(ylim = c(0, 250), clip="off") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) 

png("figures/volcano_adult.png", res = 300, height = 1200, width = 2000)
volcano_adult_plot
dev.off()

tab.adult.sig <- subset(tab.adult, adj.P.Val.com < 0.05 )
prop.table(table(chr = tab.adult$chromosome == "chrX", sig = tab.adult$adj.P.Val.com < 0.05), margin = 1)


cot <- child_comb %>%
                        mutate(Dimorphic = ifelse(adj.P.Val.com < 0.05 , "Dimorphic", "Non-dimorphic"),
                               Direction = ifelse(logFC.com < 0, "Girls", "Boys"))
prop.table(table(chr = cot$chr == "chrX", sig = cot$Dimorphic), margin = 1)

## Genes enrichment adults ####
background_adult <- tab.adult$Symbol
enrich_res_adult <- enrichr(tab.adult.sig$Symbol, db.vec, background = background_adult,
  include_overlap = TRUE)
names(enrich_res_adult) <- db.vec

enrich_res_df.comb_adult <-  Reduce(rbind,enrich_res_adult) %>%
  mutate(DB = rep(names(enrich_res_adult), sapply(enrich_res_adult, nrow))) %>%
  filter(Adjusted.P.value < 0.05) 

enrich_genes.comb_adult <- enrich_res_df.comb_adult %>% 
  mutate(Type = ifelse(DB == "DisGeNET", "Phenotypes and Diseases", "Transcription Factor"),
          Overlap = Overlap, 
          OR = Odds.Ratio, 
          pvalue = P.value, 
          adj_pvalue = Adjusted.P.value) %>%
  select(Term, DB, Type, Overlap, OR, pvalue, adj_pvalue)


### Sup Table 3
enrich_genes.comb <- 
  left_join(
    Reduce(rbind,enrich_res_child) %>%
      mutate(DB = rep(names(enrich_res_child), sapply(enrich_res_child, nrow))) %>%
      mutate(Overlap_child = Overlap,
              OR_child = Odds.Ratio,
              pvalue_child = P.value,
              adj_pvalue_child = Adjusted.P.value)  %>%
      dplyr::select(Term, ends_with("child"), DB),
    Reduce(rbind, enrich_res_adult) %>%
      mutate(DB = rep(names(enrich_res_adult), sapply(enrich_res_adult, nrow))) %>%
      mutate(Overlap_adult = Overlap,
              OR_adult = Odds.Ratio,
              pvalue_adult = P.value,
              adj_pvalue_adult = Adjusted.P.value)  %>%
      dplyr::select(Term, ends_with("adult"), DB),
  by = c("Term", "DB")) %>%
  filter(adj_pvalue_child < 0.05 | adj_pvalue_adult < 0.05)


enrich_genes.comb_out <- enrich_genes.comb %>%
  mutate(Type = ifelse(DB == "DisGeNET", "Phenotypes and Diseases", "Transcription Factor")) %>%
  select(Term, DB, Type, ends_with("child"), ends_with("adult"))

write.table(enrich_genes.comb_out, file = "results/genes/dimorphic_gene_enrichr_comb.txt",
            sep = "\t", quote = FALSE, row.names = FALSE) 


eul_GS_df <- enrich_genes.comb_out %>%
 mutate(Adults = adj_pvalue_adult < 0.05,
        Children = adj_pvalue_child < 0.05)

## Figure 3F
overlap_TF_plot <- plot(euler(eul_GS_df[eul_GS_df$Type == "Transcription Factor", c("Adults", "Children")], shape = "ellipse"),
                     fills = list(fill = c("blue", "orange"), alpha = 0.5),
                     quantities = list(fontsize = 12),
                     main = "Transcription Factors")

## Figure 3E
overlap_Disease_plot <- plot(euler(eul_GS_df[eul_GS_df$Type != "Transcription Factor", c("Adults", "Children")], shape = "ellipse"),
                             fills = list(fill = c("blue", "orange"), alpha = 0.5),
                             quantities = list(fontsize = 12),
                             main = "Phenotypes and Diseases")

png("figures/adult_vs_children_GS_overlap.png")
plot_grid(overlap_TF_plot, overlap_Disease_plot)
dev.off()


enrich_genes_age2 <- rbind(Reduce(rbind,enrich_res_child) %>%
                            mutate(DB = rep(names(enrich_res_child), sapply(enrich_res_child, nrow)),
                            Dataset = "Children"),
                          Reduce(rbind,enrich_res_adult) %>%
                            mutate(DB = rep(names(enrich_res_adult), sapply(enrich_res_adult, nrow)),
                              Dataset = "Adults")
) %>%
as_tibble() %>%
select(-Genes)


sel_terms <- enrich_genes_age2 %>%
  group_by(Dataset, DB) %>%
  filter(Adjusted.P.value < 0.05) %>%
  slice_min(P.value, n = 4) %>%
  arrange(DB, Dataset) %>%
  pull(Term) %>%
  unique()


## Figure 3H
topGeneSets <- enrich_genes_age2 %>%
  filter(Term %in% c(sel_terms, "Autistic Disorder")) %>%
  mutate(Term = factor(Term, levels = unique(c("Autistic Disorder", sel_terms))),
        logP = -log10(P.value),
        logP = ifelse(Term == "GATA1 CHEA" & Dataset == "Children", logP + 0.3, logP),
        OR = Odds.Ratio) %>%
  ggplot(aes(x = logP, y = Term, color = Dataset, size = OR)) +
  geom_point() +
  xlab(expression(-log[10]("p"))) + 
  scale_color_manual(name = "", values = c("blue", "orange")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))

png("figures/enrichment_terms_combined.png", width = 268, height = 156, units = "mm", res = 300)
topGeneSets
dev.off()





### Descriptives ####
comb_descrip <- rbind(child_comb %>%
                        mutate(Dimorphic = ifelse(adj.P.Val.com < 0.05 , "Dimorphic", "Non-dimorphic"),
                               Direction = ifelse(logFC.com < 0, "Girls", "Boys")) %>%
                        dplyr::select(Symbol, Dimorphic, Direction, chr) %>%
                        mutate(Dataset = "Children"),
                      tab.adult %>%
                        mutate(Dimorphic = ifelse(adj.P.Val.com < 0.05, "Dimorphic", "Non-dimorphic"),
                               Direction = ifelse(logFC.com < 0, "Girls", "Boys"),
                              chr = chromosome,
                               Dataset = "Adults") %>%
                        dplyr::select(Symbol, Dimorphic, Direction, chr, Dataset)) %>%
                        as_tibble() %>%
  mutate(Dataset = factor(Dataset, levels = c("Children", "Adults")))
## Figure 3E
comb_sig_prop <- comb_descrip %>%
  mutate(Autosomal = ifelse(chr == "chrX", "X-chromosome", "Autosomal")) %>%
  group_by(Dataset, Autosomal) %>%
  summarize(Dimorphic = mean(Dimorphic == "Dimorphic")) %>%
  ggplot(aes(x = Autosomal, y = Dimorphic*100, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Dimorphic genes (%)") +
  xlab("Gene location") +
  scale_fill_manual(name = "", values = c( "orange", "blue")) +
  theme_bw()
png("figures/adult_comparison.png", res = 300, height = 1200, width = 2000)
comb_sig_prop
dev.off()

comb_descrip2 <- inner_join(child_comb %>%
                        mutate(Dimorphic = ifelse(adj.P.Val.com < 0.05, "Dimorphic", "Non-dimorphic"),
                               Direction = ifelse(logFC.com < 0, "Girls", "Boys")) %>%
                        dplyr::select(Symbol, Dimorphic, Direction, chr),
                      tab.adult %>%
                        mutate(Dimorphic = ifelse(adj.P.Val.com < 0.05, "Dimorphic", "Non-dimorphic"),
                               Direction = ifelse(logFC.com < 0, "Girls", "Boys")) %>%
                        dplyr::select(Symbol, Dimorphic, Direction),
                      by = "Symbol", suffix = c(".Children", ".Adults") )  %>%
    mutate(Dimorphic_comb = ifelse(Dimorphic.Children == "Dimorphic", 
                                   ifelse(Dimorphic.Adults == "Dimorphic", "Shared", "Children"),
                                   ifelse(Dimorphic.Adults == "Dimorphic", "Adults", "None"))) %>%
  filter(Dimorphic_comb != "None")

### Sex comparison ####
## Figure 3F
comb_girl_prop <- comb_descrip %>%
  filter(Dimorphic == "Dimorphic") %>%
  mutate(Autosomal = ifelse(chr == "chrX", "X-chromosome", "Autosomal")) %>%
  group_by(Dataset, Autosomal) %>%
  summarize(girls_prop = mean(Direction == "Girls")) %>%
  ggplot(aes(x = Autosomal, y = girls_prop*100, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Women Overexpression (%)") +
  xlab("Dimorphic genes locus") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  scale_fill_manual(name = "", values = c( "orange", "blue")) +
  theme_bw()
png("figures/adult_comparison_girlsprop.png", res = 300, height = 1200, width = 2000)
comb_girl_prop
dev.off()


## Figure 3
png("figures/dimorphic_gexp_panel2.png", res = 300, height = 3600, width = 4000)
plot_grid(
  plot_grid(man_plot_children, volcano_plot_children, ncol = 2, rel_widths = c(1.8, 1), labels = c("A", "B")),
  plot_grid(man_adult_plot, volcano_adult_plot,  ncol = 2, rel_widths = c(1.8, 1), labels = c("C", "D")),
  plot_grid(
    plot_grid(comb_sig_prop, comb_girl_prop, overlap_Disease_plot, overlap_TF_plot, overlap_Disease_plot,
      ncol = 2, nrow = 2, labels = c("E", "F", "G", "")),
    topGeneSets, ncol = 2, rel_widths = c(1.5, 1), labels = c("", "H")),
  nrow = 3, rel_heights = c(1, 1, 2)
)
dev.off()

## Effect plots
## Sup Fig 2
effects_comb <- inner_join(child_comb %>%
                        mutate(Dimorphic = ifelse(adj.P.Val.com < 0.05, "Dimorphic", "Non-dimorphic")) %>%
                        dplyr::select(Symbol, Dimorphic, chr, logFC.com, P.Value.com),
                      tab.adult %>%
                        mutate(Dimorphic = ifelse(adj.P.Val.com < 0.05, "Dimorphic", "Non-dimorphic")) %>%
                      dplyr::select(Symbol, Dimorphic, logFC.com, P.Value.com),
                      by = "Symbol", suffix = c(".Children", ".Adults") )  %>%
    mutate(Dimorphic_comb = ifelse(Dimorphic.Children == "Dimorphic", 
                                   ifelse(Dimorphic.Adults == "Dimorphic" | P.Value.com.Adults < 0.05, "Shared", "Children"),
                                   ifelse(Dimorphic.Adults == "Dimorphic", 
                                    ifelse(P.Value.com.Children < 0.05, "Shared", "Adults"), "None"))) 

png("figures/adult_vs_children_effect.png", width = 1500, height = 1500, res = 300)
ggplot(effects_comb, 
  aes(x = logFC.com.Children, y = logFC.com.Adults, color = Dimorphic_comb)) +
  geom_point(alpha = 0.7) +
  xlab("logFC Children") +
  ylab("logFC Adults") +
  scale_color_manual(values = c("blue", "orange", "grey", "pink")) +
  facet_wrap(~Dimorphic_comb) +
  theme_bw() +
  stat_cor(method = "pearson", label.x = -1.2, label.y = -1.8, 
  color = "black",
  label.background = element_rect(fill = "white", color = NA)) +
  theme(legend.position = "none")
dev.off()