#'##############################################################################
#' Dimorphic genes 
#' docker run -it -v $PWD:$PWD -w $PWD dimorphism_r:1.5 R
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
library(sva)
library(parallel)
library(compositions)

## Load and preprocess gene expression ####
load("./results/preprocessFiles/Expression_SE_raw_filtered_SV.RData")
load("data/gexpAnnotation.Rdata")

# se.filt$cohort <- droplevels(se.filt$cohort)

xci_genes <- read_xlsx("data/XCI_genes.PMID29022598.xlsx", skip = 1) 
xci_genes <- xci_genes[, 1:7]
colnames(xci_genes) <- c("Symbol", "Ensembl", "Chr", "Start", "End", "Type", "XCI")

## Dimorphic TCs ####
### Remove genes in chrY 
se.noY <- se.filt[seqnames(se.filt) != "chrY", ]
# pseudo_reg <- GRanges(c("chrX:60001-2699520", "chrX:154931044-155260560")) ## hg19
# pseudo_genes <- subsetByOverlaps(rowRanges(se.noY), pseudo_reg)
# se.noY <- se.noY[!rownames(se.noY) %in% names(pseudo_genes), ]
dim(se.noY)
#[1] 30013   767

mod.helix <- model.matrix(formula(paste("~ e3_sex + cohort + age_sample_years + NK_6 +
                                        Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6 +", 
                                              paste(paste0("SV", 1:3), collapse = "+"))),
                            colData(se.noY))

lm.gexp <- lmFit(assay(se.noY), mod.helix) %>%
  eBayes()
tab.gexp <- topTable(lm.gexp, coef = 2, n = Inf, confint = TRUE)
tab.gexp$padj_holm <- p.adjust(tab.gexp$P.Value, method = "holm")

## Compute number of effective tests
# ntest_gene_dim <- meff(cor(t(assay(se.noY))), method = "galwey")


tab.gexp$TC <- rownames(tab.gexp)
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

## Check consistency between HELIX cohorts
cohorts <- unique(as.character(se.noY$cohort))
cohorts_name <- recode(cohorts, SAB = "INMA")
names(cohorts_name) <- cohorts

cohort_limma_list <- lapply(cohorts, function(coh) {
  message(coh)
  se.coh <- se.noY[, se.noY$cohort == coh]

  mat <- assay(se.coh)
  pd <- as.data.frame(colData(se.coh))

  mod <- model.matrix(~ e3_sex + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, data = pd)
  mod0 <- model.matrix(~ age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, data = pd)
  n.sv <- num.sv(mat, mod, method = "leek") + 1
  sv.obj <- sva(mat, mod, mod0 = mod0, n.sv = n.sv)

  sv.mat <- sv.obj$sv
  colnames(sv.mat) <- paste0("SV", seq_len(ncol(sv.mat)))
  rownames(sv.mat) <- rownames(pd)
  pd <- cbind(pd, sv.mat)

  mod.coh <- model.matrix(
    formula(paste(
      "~ e3_sex + age_sample_years + NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6 +",
      paste(colnames(sv.mat), collapse = " + ")
    )),
    data = pd
  )

  fit.coh <- lmFit(mat, mod.coh) %>%
    eBayes()

  tab.coh <- topTable(fit.coh, coef = 2, n = Inf, confint = TRUE)  

  list(
    table = tab.coh,
    n.sv = n.sv,
    N = ncol(se.coh)
    )
})
names(cohort_limma_list) <- cohorts

cohort_limma_tbl <- tibble(TC = rownames(cohort_limma_list[[1]]$table))
for (coh in cohorts) {
  cohort_limma_tbl[[coh]] <- cohort_limma_list[[coh]]$table[cohort_limma_tbl$TC, "logFC"]
}
cohort_limma_tbl$coherent <- apply(cohort_limma_tbl[, cohorts], 1, function(x) {
  ref_sign <- sign(sum(as.numeric(x), na.rm = TRUE))
  sum(sign(as.numeric(x)) == ref_sign, na.rm = TRUE)
}) 
cohort_limma_tbl <- mutate(cohort_limma_tbl, coherent = pmax(coherent, 6 - coherent))


common_genes <- Reduce(intersect, lapply(cohort_limma_list, function(x) rownames(x$table)))

meta_gene_list <- mclapply(common_genes, function(gene) {
  te <- sapply(cohort_limma_list, function(x) x$table[gene, "logFC"])
  low <- sapply(cohort_limma_list, function(x) x$table[gene, "CI.L"])
  up <- sapply(cohort_limma_list, function(x) x$table[gene, "CI.R"])
  n <- sapply(cohort_limma_list, function(x) x$N)

  tryCatch(
    metagen(
      TE = te,
      lower = low,
      upper = up,
      n.e = n,
      studlab = as.character(cohorts_name[names(te)]),
      sm = "MD",
      method.tau = "REML"
    ),
    error = function(e) {
      list(
        TE.random = NA_real_,
        TE.common = NA_real_,
        pval.random = NA_real_,
        pval.common = NA_real_
      )
    }
  )
}, mc.cores = 5)
names(meta_gene_list) <- common_genes
save(meta_gene_list, cohort_limma_tbl, file = "results/genes/HELIX_cohort_genes.Rdata")


meta_gene_summary <- tibble(
  TC = common_genes,
  TE_fixed = sapply(meta_gene_list, function(x) x$TE.common),
  P_fixed = sapply(meta_gene_list, function(x) x$pval.common),
  TE_random = sapply(meta_gene_list, function(x) x$TE.random),
  P_random = sapply(meta_gene_list, function(x) x$pval.random)
)


gene_cohort_summary <- left_join(tab.gexp, meta_gene_summary, by = "TC") %>%
  left_join(select(cohort_limma_tbl, TC, coherent), by = "TC") %>%
  mutate(
    var_TE_fixed = ifelse(logFC == 0, NA_real_, 100 * (TE_fixed - logFC) / abs(logFC)),
    var_TE_random = ifelse(logFC == 0, NA_real_, 100 * (TE_random - logFC) / abs(logFC)),
    sign_conc_fixed = sign(logFC) == sign(TE_fixed),
    sign_conc_random = sign(logFC) == sign(TE_random)
  ) %>%
  as_tibble() 

## Define dimorphic genes
gene_cohort_summary$dimorphic <- ifelse(gene_cohort_summary$adj.P.Val < 0.01 & gene_cohort_summary$coherent >= 5, "Dimorphic", "Non-Dimorphic")

save(gene_cohort_summary, file = "results/genes/HELIX_cohort_genes_summary.Rdata")

write.table(gene_cohort_summary %>%
  select(TC, Symbol, Coding, chr, logFC, CI.L, CI.R, P.Value, adj.P.Val, TE_fixed, P_fixed, TE_random, P_random, coherent, dimorphic), 
  file = "results/genes/HELIX_cohort_genes_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)


png("figures/Genes/EffectSizesMeta.png", width = 2000, height = 1200, res = 300)
gene_cohort_summary %>% select(logFC, TE_fixed, TE_random, dimorphic) %>%
  filter(logFC > -2) %>%
  pivot_longer(cols = c(TE_fixed, TE_random), names_to = "Method", values_to = "coef") %>%
  mutate(Method = recode(Method, TE_fixed = "Fixed-effects", TE_random = "Random-effects"),
    Effect = ifelse(dimorphic == "Dimorphic", ifelse(logFC > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
  ggplot(aes(x = logFC, y = coef, color = Effect)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("Higher Boys" = "blue", "Higher Girls" = "maroon", "No differences" = "black")) +
  stat_cor(inherit.aes = FALSE, aes(x = logFC, y = coef), method = "pearson", label.x = -1, label.y = 0.35) +
  facet_grid(dimorphic~Method) +
  xlab("Pooled logFC") +
  ylab("Meta-analysis logFC")
dev.off()


# mod.tech <- model.matrix(formula(paste("~ e3_sex + cohort + age_sample_years + NK_6 +
#                                         Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6 + bio_conc_denat_ngul +
#                                         bio_rin_denat +", 
#                                               paste(paste0("SV", 1:3), collapse = "+"))),
#                             colData(se.noY))

# lm.gexp_tech <- lmFit(assay(se.noY), mod.tech) %>%
#   eBayes()
# tab.gexp_tech <- topTable(lm.gexp_tech, coef = 2, n = Inf, confint = TRUE)
# tab.gexp_tech$TC <- rownames(tab.gexp_tech)
# tab.gexp_tech$Dimorphic <- ifelse(tab.gexp_tech$TC %in% subset(gene_cohort_summary, dimorphic == "Dimorphic")$TC, "Dimorphic", "Non-Dimorphic")
# save(tab.gexp_tech, file = "results/genes/HELIX_dimorphic_genes_tech.Rdata")

## Results exploration
dim.gexp <- subset(gene_cohort_summary, dimorphic == "Dimorphic")$TC

### All TCs ####
## N genes
sum(gene_cohort_summary$dimorphic == "Dimorphic")
mean(gene_cohort_summary$dimorphic == "Dimorphic")


## Coding
coding_tab <- table(sig = gene_cohort_summary$dimorphic == "Dimorphic", coding = gene_cohort_summary$Coding)
prop.table(coding_tab, margin = 1)
prop.table(coding_tab, margin = 2)

fisher.test(xchr_tab)

## X Chromosome prop
xchr_tab <- table(sig = gene_cohort_summary$dimorphic == "Dimorphic", chr = gene_cohort_summary$chr == "chrX")
prop.table(xchr_tab, margin = 1)
prop.table(xchr_tab, margin = 2)

fisher.test(xchr_tab)

### XCI
child_comb_chrX <- subset(gene_cohort_summary, chr == "chrX")
child_comb_chrX <- left_join(child_comb_chrX, xci_genes[, c("Symbol", "Type", "XCI")],
                            by = "Symbol")
child_comb_chrX$XCI2 <- ifelse(child_comb_chrX$XCI %in% c("inactive", "variable"), "Inactive", child_comb_chrX$XCI)

xci_tab <- table(child_comb_chrX$XCI, sig = child_comb_chrX$dimorphic == "Dimorphic")
chisq.test(xci_tab)
sum(!is.na(child_comb_chrX$XCI) & child_comb_chrX$dimorphic == "Dimorphic")

fisher.test(xci_tab[c(2, 1), ])

child_comb_chrX_sig <- subset(child_comb_chrX, dimorphic == "Dimorphic")
mean(child_comb_chrX_sig$XCI == "escape", na.rm = TRUE)


# Proportion of genes DE by sex
gene_cohort_sig <- subset(gene_cohort_summary, dimorphic == "Dimorphic")
mean(gene_cohort_sig$logFC < 0)
prop.table(table(ifelse(gene_cohort_sig$logFC < 0, "Girls", "Boys"), 
  ifelse(gene_cohort_sig$chr == "chrX", "chrX", "Autosomal")), margin = 2)

prop.table(
  table(ifelse(child_comb_chrX_sig$logFC < 0, "Girls", "Boys"), 
   child_comb_chrX_sig$XCI2)
, margin = 2)
gene_cohort_summary %>%
  arrange(P.Value) %>%
  head(40) %>%
  data.frame() %>%
  select(Symbol, chr,logFC, adj.P.Val )

## Enrichments ####
dbs <- listEnrichrDbs()

db.vec <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "DisGeNET")

enrich_genes_comb <- subset(gene_cohort_summary, dimorphic == "Dimorphic")$Symbol
background_child <-  gene_cohort_summary$Symbol
enrich_res_child <- enrichr(enrich_genes_comb, db.vec, background = background_child,
  include_overlap = TRUE)
names(enrich_res_child) <- db.vec
enrich_child_df <-  Reduce(rbind,enrich_res_child) %>%
  mutate(DB = rep(names(enrich_res_child), sapply(enrich_res_child, nrow))) %>%
  filter(Adjusted.P.value < 0.05) 


# enrich_genesX_comb <- subset(gene_cohort_sig, chr == "chrX")$Symbol
# background_childX <-  subset(gene_cohort_summary, chr == "chrX")$Symbol
# enrich_resX_comb <- enrichr(enrich_genesX_comb, db.vec,
#   background = background_childX,
#   include_overlap = TRUE)
# enrich_combX_df <-  Reduce(rbind,enrich_resX_comb) %>%
#   mutate(DB = rep(names(enrich_resX_comb), sapply(enrich_resX_comb, nrow))) %>%
#   filter(Adjusted.P.value < 0.05) 

enrich_genes.children <- enrich_child_df %>% 
  mutate(Type = ifelse(DB == "DisGeNET", "Phenotypes and Diseases", "Transcription Factor"),
          Overlap = Overlap, 
          OR = Odds.Ratio, 
          pvalue = P.value, 
          adj_pvalue = Adjusted.P.value) %>%
  select(Term, DB, Type, Overlap, OR, pvalue, adj_pvalue)
write.table(enrich_genes.children, file = "results/genes/dimorphic_gene_enrichr_helix.txt",
            sep = "\t", quote = FALSE, row.names = FALSE) 

child_TFs <- subset(enrich_child_df, DB == "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")$Term  
child_TFs <- unique(gsub(" ENCODE", "", gsub(" CHEA", "", child_TFs)))

subset(gene_cohort_summary, Symbol %in% child_TFs) %>% 
  dplyr::select(Symbol, adj.P.Val, logFC, P.Value, dimorphic) %>%
  arrange(P.Value) %>%
  data.frame()


ar_genes <- strsplit(subset(enrich_child_df, Term == "AR CHEA")$Genes, ";")[[1]]
subset(gene_cohort_summary, Symbol %in% ar_genes) %>% 
  dplyr::select(Symbol, adj.P.Val, logFC, P.Value, dimorphic) %>%
  arrange(P.Value) %>%
  data.frame()



### Manhattan ####
## Figure 3A
man_tab_child <- gene_cohort_summary 
man_tab_child$Position <- start(rowRanges(se.filt[man_tab_child$TC, ]))
man_tab_child$chr <- factor(man_tab_child$chr, levels = paste0("chr", c(1:22, "X")))
man_tab_child$chrNum <- factor(gsub("chr", "", man_tab_child$chr), levels = c(1:22, "X"))
man_tab_child$Genes <- man_tab_child$Symbol
man_tab_child$Genes <- ifelse(man_tab_child$P.Value < 2e-26, man_tab_child$Genes, "")
man_tab_child$Genes[man_tab_child$Genes == "NA"] <- ""
man_tab_child$pval <- ifelse(man_tab_child$P.Value < 1e-200, 1e-200, man_tab_child$P.Value)


man_plot_children <- manhattan_plot(x = man_tab_child, pval.colname = "pval", 
                           signif = c(0.05/nrow(man_tab_child)),
                           chr.colname = "chrNum", pos.colname = "Position",
                           plot.title = "Children",
                           rescale = FALSE,
                           label.colname = "Genes") +
  theme(plot.title = element_text(hjust = 0.5)) 


man_plot_children$layers[[2]] <- NULL ## Remove significance line
png("figures/genes/manhattan_children.png", res = 300, height = 1200, width = 2000)
man_plot_children
dev.off()

### Volcano ####
### Figure 3B
volcano_plot_children <- man_tab_child %>%
  mutate(Signif = ifelse(dimorphic == "Dimorphic", ifelse(logFC > 0, "Boys", "Girls"), "None")) %>%
  ggplot(aes(x = logFC, y = -log10(pval), color = Signif, label = Genes)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_label_repel(max.overlaps = 15) +
  scale_color_manual(values = c("blue", "maroon", "black")) +
  theme_bw() +
  xlab("Effect size") +
  ylab(expression(-log[10]("p"))) + 
  annotate("text", x = -1, y = -30, label = "Higher Girls") + 
  annotate("text", x = 1, y = -30, label = "Higher Boys") + 
  coord_cartesian(ylim = c(0, 230), xlim = c(-1.5, 1.5), clip = "on") +
  ggtitle("Children") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

png("figures/genes/volcano_children.png", res = 300, height = 1200, width = 2000)
volcano_plot_children
dev.off()


# ## Comparison with children ####
# ### Controls of T1D study
# load("results/preprocessFiles/gexp_GSE43488_sex.Rdata")
# ### Run model ####
# mod.gse43488 <- model.matrix(formula(paste("~ sex + age +", 
#                                               paste(paste0("SV", 1:2), collapse = "+"))),
#                             colData(gse43488_controls ))
# lm.gse43488 <- lmFit(t(scale(t(assay(gse43488_controls)))), mod.gse43488) %>%
#   eBayes()
# tab.gse43488 <- topTable(lm.gse43488, coef = 2, n = Inf, confint = TRUE)
# tab.gse43488$Symbol <- rownames(tab.gse43488)
# save(tab.gse43488, file = "results/genes/gse43488_analysis.Rdata")

## Comparison with Generation R ####
load("results/genes/GenerationR_analysis_2026.Rdata")



## Sup File 2
write.table(tab.genR_main %>%
  mutate(Chromosome = chromosome) %>%
  select(Symbol, Chromosome, logFC, CI.L, CI.R, P.Value, adj.P.Val), 
  file = "results/genes/GenR_analysis.txt", sep = "\t", 
  quote = FALSE, row.names = FALSE)

### Combine with HELIX ####
child_comb <- gene_cohort_summary  %>% 
  filter(!is.na(Symbol) & Symbol != "" & Symbol != "NA") %>%
  left_join(dplyr::select(tab.genR_main, Symbol, logFC,  P.Value, adj.P.Val) %>%
              group_by(Symbol) %>% filter(P.Value == min(P.Value)), by = "Symbol",
            suffix = c(".helix", ".GenR")) %>%
  filter(!is.na(P.Value.helix) & !is.na(P.Value.GenR))

com_genes_child <- nrow(child_comb)


child_comb_plot <- child_comb %>%
  mutate(Effect = ifelse(dimorphic == "Dimorphic", ifelse(logFC.helix > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
  ggplot(aes(x = logFC.helix, y = logFC.GenR, color = Effect)) +
  geom_point() +
  stat_cor(data = child_comb, inherit.aes = FALSE, aes(x = logFC.helix, y = logFC.GenR), method = "pearson", label.x = -1.45, label.y = 1.85) +
  theme_bw() +
  scale_color_manual(values = c("Higher Boys" = "blue", "Higher Girls" = "maroon", "No differences" = "black")) +
  facet_grid(dimorphic ~.) +
  xlab("HELIX logFC") +
  ylab("GenR logFC") +
  coord_cartesian(ylim = c(-4, 2), xlim = c(-1.5, 1.5), clip = "off") 

png("figures/Genes/EffectSizesChild.png", width = 2000, height = 1200, res = 300)
child_comb_plot
dev.off()




# #### Run meta-analysis in common genes ####
# meta_list <- lapply(seq_len(nrow(child_comb)), function(i) {
#   metagen(TE = c(child_comb$logFC.helix[i], child_comb$logFC.gse43488[i], child_comb$logFC[i]), 
#           upper = c(child_comb$CI.R.helix[i], child_comb$CI.R.gse43488[i], child_comb$CI.R[i]), 
#           lower = c(child_comb$CI.L.helix[i], child_comb$CI.L.gse43488[i], child_comb$CI.L[i]))
# })
# child_comb$logFC.com <- sapply(meta_list, function(x) ifelse(x$pval.Q < 0.05/com_genes_child, x$TE.random,x$TE.common))
# child_comb$P.Value.com <- sapply(meta_list, function(x) ifelse(x$pval.Q < 0.05/com_genes_child, x$pval.random, x$pval.common))
# child_comb$adj.P.Val.com <- p.adjust(child_comb$P.Value.com, method = "BH")
# save(child_comb, file = "results/genes/children_analysis.Rdata")


# ## Sup File 3
# write.table(child_comb %>%
#   mutate(logFC.T1D = logFC.gse43488,
#           P.Value.T1D = P.Value.gse43488,
#           adj.P.Val.T1D = adj.P.Val.gse43488,
#           CI.L.T1D = CI.L.gse43488,
#           CI.R.T1D = CI.R.gse43488,
#           logFC.GenR = logFC,
#           P.Value.GenR = P.Value,
#           adj.P.Val.GenR = adj.P.Val,
#           CI.L.GenR = CI.L,
#           CI.R.GenR = CI.R,
#           logFC.Meta = logFC.com,
#           P.Value.Meta = P.Value.com,
#           adj.P.Val.Meta = adj.P.Val.com) %>%
#   select(Symbol, chr, TC, ends_with("helix"), ends_with("T1D"), ends_with("GenR"), ends_with("Meta")), file = "results/genes/children_metaanalysis.txt", sep = "\t", 
#   quote = FALSE, row.names = FALSE)


# Comparison with adults ####


## SHIP
load("results/genes/ship_analysis.Rdata")

## RSIII
rs3_raw <- read_xlsx("results/genes/RSIII_TWAS_results_raw_v2.xlsx")
# rs3_scaled <- read_xlsx("results/genes/RSIII_TWAS_results_scaled.xlsx")


# ## GSE36382 ####
# load("results/preprocessFiles/gexp_gse36382_sex.Rdata")

# ### Run model ####
# ### Select probes mapped to chr1-22 and X and probes with symbol
# gse36382_noY <- gse36382_se[as.character(seqnames(gse36382_se)) != "chrY", ]
# mod_gse36382 <- model.matrix(~ sex, data = colData(gse36382_se))

# assay(gse36382_noY)$E <- t(scale(t(assay(gse36382_noY)$E)))
# lm.gse36382 <- lmFit(assay(gse36382_noY), mod_gse36382) %>%
#   eBayes()
# tab.gse36382 <- topTable(lm.gse36382, coef = 2, n = Inf, confint = TRUE)
# tab.gse36382$Symbol <- rownames(tab.gse36382)
# tab.gse36382$chromosome <- as.character(seqnames(gse36382_noY[rownames(tab.gse36382), ]))
# save(tab.gse36382, file = "results/genes/gse36382_analysis.Rdata")

# ## GSE33828 ####
# load("results/preprocessFiles/gexp_gse33828_sex.Rdata")
# ### Run model ####
# gse33828_noY <- gse33828_se[as.character(seqnames(gse33828_se)) != "chrY", ]
# mod_gse33828 <- model.matrix(~ sex + age, data = colData(gse33828_noY))

# lm.gse33828 <- lmFit(t(scale(t(assay(gse33828_noY)))), mod_gse33828) %>%
#   eBayes()
# tab.gse33828 <- topTable(lm.gse33828, coef = 2, n = Inf, confint = TRUE)
# tab.gse33828$Symbol <- rownames(tab.gse33828)
# tab.gse33828$chromosome <- as.character(seqnames(gse33828_noY[rownames(tab.gse33828), ]))
# save(tab.gse33828, file = "results/genes/gse33828_analysis.Rdata")

## Combine adult results ####
### Sup Table 3
tab.array <- inner_join(rs3_raw, tab.ship.main, by = c("GeneSymbol" = "ILMN_Gene"), 
                        suffix = c(".ship", ".rs3")) 

meta_list_array <- mclapply(seq_len(nrow(tab.array)), function(i) {
  metagen(TE = c(tab.array$logFC[i],  tab.array$EffectSize[i]), 
          upper = c(tab.array$CI.R[i],  tab.array$CI_high[i]), 
          lower = c(tab.array$CI.L[i], tab.array$CI_low[i]),
          n.e = c(987, 765),
          studlab = c("SHIP",  "RSIII"), 
          sm = "MD", method.tau = "REML")
}, mc.cores = 3)
tab.array$logFC.fixed    <- sapply(meta_list_array, function(x) x$TE.common)
tab.array$P.Value.fixed  <- sapply(meta_list_array, function(x) x$pval.common)
tab.array$adj.P.Val.fixed <- p.adjust(tab.array$P.Value.fixed, method = "BH")

tab.array$logFC.random   <- sapply(meta_list_array, function(x) x$TE.random)
tab.array$P.Value.random <- sapply(meta_list_array, function(x) x$pval.random)
tab.array$adj.P.Val.random <- p.adjust(tab.array$P.Value.random, method = "BH")

tab.array$I2             <- sapply(meta_list_array, function(x) x$I2)
tab.array$Q.pval         <- sapply(meta_list_array, function(x) x$pval.Q)

tab.array$logFC.com    <-  ifelse(tab.array$I2 < 0.4, tab.array$logFC.fixed, tab.array$logFC.random)
tab.array$P.Value.com    <-  ifelse(tab.array$I2 < 0.4, tab.array$P.Value.fixed, tab.array$P.Value.random)

tab.array$Dimorphic <- ifelse(tab.array$adj.P.Val.random < 0.01 | (tab.array$adj.P.Val.fixed < 0.01 & tab.array$I2 < 0.4), "Dimorphic", "Non-Dimorphic")
save(tab.array, meta_list_array, file = "results/genes/array_analysis.Rdata")
nrow(tab.array)

group_by(tab.array, Dimorphic == "Dimorphic") %>%
  summarize(cor = cor(logFC, EffectSize, use = "complete.obs"))

table(tab.array$Dimorphic, tab.array$logFC > 0)

# tab.adult <- inner_join(tab.ship.scaled, tab.gtex_scaled, by = c("ILMN_Gene" = "Symbol"), 
#                         suffix = c(".ship", ".gtex")) %>%
#   inner_join(rs3_scaled, by = c("ILMN_Gene" = "GeneSymbol"), suffix = c(".array", ".rs3")) %>%
#   as_tibble()

# com_genes_adult <- nrow(tab.adult)

# meta_list_adult <- mclapply(seq_len(nrow(tab.adult)), function(i) {
#   metagen(TE = c(tab.adult$logFC.ship[i], tab.adult$logFC.gtex[i], tab.adult$EffectSize[i]), 
#           upper = c(tab.adult$CI.R.ship[i], tab.adult$CI.R.gtex[i], tab.adult$CI_high[i]), 
#           lower = c(tab.adult$CI.L.ship[i], tab.adult$CI.L.gtex[i], tab.adult$CI_low[i]),
#           n.e = c(987, 881, 765),
#           studlab = c("SHIP", "GTEx", "RSIII"), 
#           sm = "MD", method.tau = "REML")
# }, mc.cores = 3)
# tab.adult$logFC.fixed    <- sapply(meta_list_adult, function(x) x$TE.common)
# tab.adult$P.Value.fixed  <- sapply(meta_list_adult, function(x) x$pval.common)
# tab.adult$adj.P.Val.fixed <- p.adjust(tab.adult$P.Value.fixed, method = "BH")

# tab.adult$logFC.random   <- sapply(meta_list_adult, function(x) x$TE.random)
# tab.adult$P.Value.random <- sapply(meta_list_adult, function(x) x$pval.random)
# tab.adult$adj.P.Val.random <- p.adjust(tab.adult$P.Value.random, method = "BH")

# tab.adult$I2             <- sapply(meta_list_adult, function(x) x$I2)
# tab.adult$Q.pval         <- sapply(meta_list_adult, function(x) x$pval.Q)

# tab.adult$logFC.com    <-  ifelse(tab.adult$I2 < 0.4, tab.adult$logFC.fixed, tab.adult$logFC.random)
# tab.adult$P.Value.com    <-  ifelse(tab.adult$I2 < 0.4, tab.adult$P.Value.fixed, tab.adult$P.Value.random)

# tab.adult$Dimorphic <- ifelse(tab.adult$adj.P.Val.random < 0.01 | (tab.adult$adj.P.Val.fixed < 0.01 & tab.adult$I2 < 0.4), "Dimorphic", "Non-Dimorphic")
# save(tab.adult, meta_list_adult, file = "results/genes/adult_analysis.Rdata")

# group_by(tab.adult, adj.P.Val.fixed < 0.01) %>%
#   summarize(cor_helix_ship = cor(logFC.ship, logFC.gtex, use = "complete.obs"),
#             cor_helix_rs3 = cor(EffectSize, logFC.gtex, use = "complete.obs"),
#             cor_ship_rs3 = cor(logFC.ship, EffectSize, use = "complete.obs"))

# group_by(tab.adult, adj.P.Val.random < 0.01) %>%
#   summarize(cor_helix_ship = cor(logFC.ship, logFC.gtex, use = "complete.obs"),
#             cor_helix_rs3 = cor(EffectSize, logFC.gtex, use = "complete.obs"),
#             cor_ship_rs3 = cor(logFC.ship, EffectSize, use = "complete.obs"))

# group_by(tab.adult, Dimorphic == "Dimorphic") %>%
#   summarize(cor_helix_ship = cor(logFC.ship, logFC.gtex, use = "complete.obs"),
#             cor_helix_rs3 = cor(EffectSize, logFC.gtex, use = "complete.obs"),
#             cor_ship_rs3 = cor(logFC.ship, EffectSize, use = "complete.obs"))

## Gene annotation
annot <- read_delim("data/GPL10558-50081.txt", comment = "#", delim = "\t")
annot_short <- select(annot, Symbol, Chromosome, Probe_Coordinates) %>%
  group_by(Symbol) %>%
  arrange(Symbol) %>%
  slice_head(n = 1) %>%
  mutate(Chromosome = paste("chr", Chromosome, sep = ""))
### Make plots #####
tab.array_chr <- left_join(tab.array, annot_short, by = c("GeneSymbol" = "Symbol"))

write.table(tab.array_chr %>%
  mutate(logFC.RS3 = EffectSize,
          P.Value.RS3 = PValue,
          adj.P.Val.RS3 = FDR,
          CI.L.RS3 = CI_low,
          CI.R.RS3 = CI_high,
          logFC.ship = logFC,
          P.Value.ship = P.Value,
          adj.P.Val.ship = adj.P.Val,
          CI.L.ship = CI.L,
          CI.R.ship = CI.R,
          Symbol = GeneSymbol) %>%
  select(Symbol, Chromosome, Probe_Coordinates, ends_with("RS3"), ends_with("ship"), ends_with("fixed"), 
  ends_with("random"), I2, Dimorphic), 
  file = "results/genes/adult_metaanalysis.txt", sep = "\t", 
  quote = FALSE, row.names = FALSE)


png("figures/Genes/EffectSizesAdult.png", width = 1600, height = 1200, res = 300)
tab.array %>%
  mutate(Effect = ifelse(Dimorphic == "Dimorphic", ifelse(logFC.com > 0, "Higher Men", "Higher Women"), "No differences")) %>%
  ggplot(aes(x = logFC, y = EffectSize, color = Effect)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("Higher Men" = "blue", "Higher Women" = "maroon", "No differences" = "black")) +
  xlab("SHIP logFC") +
  ylab("RSIII logFC") +
  facet_grid(Dimorphic ~.) +
  stat_cor(data = tab.array, inherit.aes = FALSE, aes(x = logFC, y = EffectSize), method = "pearson", label.x = -0.5, label.y = 0.7) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Manhattan and volcano plots ####
man.adult <- tab.array_chr
man.adult$Chromosome <- gsub("chr", "", man.adult$Chromosome)
man.adult$Chromosome <- ifelse(man.adult$Chromosome %in% c("XY", "Y", ""), "X", man.adult$Chromosome)
man.adult$chrNum <- factor(man.adult$Chromosome, levels = c(1:22, "X"))

man.adult$Position <- sapply(strsplit(man.adult$Probe_Coordinates, "-"), function(x) as.numeric(x[1]))
man.adult$pval <- ifelse(man.adult$P.Value.com < 1e-200, 1e-200, man.adult$P.Value.com)
man.adult$Genes <- man.adult$GeneSymbol
man.adult$Genes <- ifelse(man.adult$pval < 2e-26, man.adult$Genes, "")
man.adult$Genes[man.adult$Genes == "NA"] <- ""

## Figure 3C
man_adult_plot <- manhattan_plot(x = man.adult, pval.colname = "pval", 
                           signif = 0.05/nrow(man.adult),
                           chr.colname = "chrNum", pos.colname = "Position",
                           plot.title = "Adults",
                           rescale = FALSE,
                           label.colname = "Genes") +
  theme(plot.title = element_text(hjust = 0.5))


man_adult_plot$layers[[2]] <- NULL ## Remove significance line

png("figures/genes/manhattan_adult.png", res = 300, height = 1200, width = 2000)
man_adult_plot
dev.off()

### Figure 3D
volcano_adult_plot <- man.adult %>%
  mutate(Signif = ifelse(Dimorphic == "Dimorphic", ifelse(logFC.com > 0, "Boys", "Girls"), "None")) %>%
  ggplot(aes(x = logFC.com, y = -log10(pval), color = Signif, label = Genes)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_label_repel() +
  scale_color_manual(values = c("blue", "maroon", "black")) +
  theme_bw() +
  xlab("Combined logFC") +
  ggtitle("Adults") +
  ylab(expression(-log[10]("p"))) + 
  annotate("text", x = -0.5, y = -45, label = "Higher Women") + 
  annotate("text", x = 0.5, y = -45, label = "Higher Men") + 
  coord_cartesian(ylim = c(0, 250), xlim = c(-1.5, 1.5), clip = "on") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) 

png("figures/genes/volcano_adult.png", res = 300, height = 1200, width = 2000)
volcano_adult_plot
dev.off()

man.adult.sig <- subset(man.adult, Dimorphic == "Dimorphic" )
prop.table(table(chr = man.adult$Chromosome == "X", sig = man.adult$Dimorphic == "Dimorphic"), margin = 1)


# cot <- child_comb %>%
#                         mutate(Dimorphic = ifelse(adj.P.Val.com < 0.05 , "Dimorphic", "Non-dimorphic"),
#                                Direction = ifelse(logFC.com < 0, "Girls", "Boys"))
# prop.table(table(chr = cot$chr == "chrX", sig = cot$Dimorphic), margin = 1)

## Genes enrichment adults ####
background_adult <- tab.array$GeneSymbol
enrich_res_adult <- enrichr(man.adult.sig$GeneSymbol, db.vec, background = background_adult,
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
      Adults = ifelse(is.na(Adults), FALSE, Adults),
        Children = adj_pvalue_child < 0.05)

## Figure 3F
overlap_TF_plot <- plot(euler(eul_GS_df[eul_GS_df$Type == "Transcription Factor", c("Adults", "Children")], shape = "ellipse"),
                     fills = list(fill = c("forestgreen", "orange"), alpha = 0.5),
                     quantities = list(fontsize = 12),
                     main = "Transcription Factors")

## Figure 3E
overlap_Disease_plot <- plot(euler(eul_GS_df[eul_GS_df$Type != "Transcription Factor", c("Adults", "Children")], shape = "ellipse"),
                             fills = list(fill = c("forestgreen", "orange"), alpha = 0.5),
                             quantities = list(fontsize = 12),
                             main = "Phenotypes and Diseases")

png("figures/genes/adult_vs_children_GS_overlap.png")
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
        OR = Odds.Ratio,
        Dataset = factor(Dataset, levels = c("Children", "Adults"))) %>%
  ggplot(aes(x = logP, y = Term, color = Dataset, size = OR)) +
  geom_point() +
  xlab(expression(-log[10]("p"))) + 
  scale_color_manual(name = "", values = c( "orange", "forestgreen")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))

png("figures/genes/enrichment_terms_combined.png", width = 268, height = 156, units = "mm", res = 300)
topGeneSets
dev.off()


### Descriptives ####
comb_descrip <- rbind(gene_cohort_summary %>%
                        mutate(Direction = ifelse(logFC < 0, "Girls", "Boys"),
                        Dimorphic = dimorphic) %>%
                        dplyr::select(Symbol, Dimorphic, Direction, chr) %>%
                        mutate(Dataset = "Children"),
                      man.adult %>%
                        mutate(Direction = ifelse(logFC.com < 0, "Girls", "Boys"),
                              chr = Chromosome,
                               Dataset = "Adults",
                               Symbol = GeneSymbol) %>%
                        dplyr::select(Symbol, Dimorphic, Direction, chr, Dataset)) %>%
                        as_tibble() %>%
  mutate(Dataset = factor(Dataset, levels = c("Children", "Adults")))
## Figure 3E
comb_sig_prop <- comb_descrip %>%
  mutate(Autosomal = ifelse(chr %in% c("X","chrX"), "X-chromosome", "Autosomal")) %>%
  group_by(Dataset, Autosomal) %>%
  summarize(Dimorphic = mean(Dimorphic == "Dimorphic")) %>%
  ggplot(aes(x = Autosomal, y = Dimorphic*100, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Dimorphic genes (%)") +
  xlab("Gene location") +
  scale_fill_manual(name = "", values = c( "orange", "forestgreen")) +
  theme_bw()
png("figures/genes/adult_comparison.png", res = 300, height = 1200, width = 2000)
comb_sig_prop
dev.off()

comb_descrip2 <- inner_join(gene_cohort_summary %>%
                        mutate(Dimorphic = dimorphic,
                               Direction = ifelse(logFC < 0, "Girls", "Boys")) %>%
                        dplyr::select(Symbol, Dimorphic, Direction, chr),
                      tab.array_chr %>%
                        mutate(Symbol = GeneSymbol,
                               Direction = ifelse(logFC < 0, "Girls", "Boys")) %>%
                        dplyr::select(Symbol, Dimorphic, Direction),
                      by = "Symbol", suffix = c(".Children", ".Adults") )  %>%
    mutate(Dimorphic_comb = ifelse(Dimorphic.Children == "Dimorphic", 
                                   ifelse(Dimorphic.Adults == "Dimorphic", "Shared", "Children"),
                                   ifelse(Dimorphic.Adults == "Dimorphic", "Adults", "None"))) %>%
  filter(Dimorphic_comb != "None") %>%
  mutate(Children = ifelse(Dimorphic.Children == "Dimorphic", TRUE, FALSE),
         Adults = ifelse(Dimorphic.Adults == "Dimorphic", TRUE, FALSE))

overlap_genes_plot <- plot(euler(comb_descrip2[, c("Children", "Adults")], shape = "ellipse"),
                     fills = list(fill = c("orange", "forestgreen" ), alpha = 0.5),
                     quantities = list(fontsize = 12),
                     main = "Genes")

png("figures/genes/overlap_genes.png", res = 300, height = 1200, width = 2000)
overlap_genes_plot
dev.off()

### Sex comparison ####
## Figure 3F
comb_girl_prop <- comb_descrip %>%
  filter(Dimorphic == "Dimorphic") %>%
  mutate(Autosomal = ifelse(chr %in% c("X","chrX"), "X-chromosome", "Autosomal")) %>%
  group_by(Dataset, Autosomal) %>%
  summarize(girls_prop = mean(Direction == "Girls")) %>%
  ggplot(aes(x = Autosomal, y = girls_prop*100, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Women Overexpression (%)") +
  xlab("Dimorphic genes locus") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  scale_fill_manual(name = "", values = c( "orange", "forestgreen")) +
  theme_bw()
png("figures/genes/adult_comparison_girlsprop.png", res = 300, height = 1200, width = 2000)
comb_girl_prop
dev.off()




## Figure 3
png("figures/genes/dimorphic_gexp_panel2.png", res = 300, height = 3600, width = 4400)
plot_grid(
  plot_grid(man_plot_children, volcano_plot_children, ncol = 2, rel_widths = c(1.8, 1), labels = c("A", "B")),
  plot_grid(man_adult_plot, volcano_adult_plot,  ncol = 2, rel_widths = c(1.8, 1), labels = c("C", "D")),
  plot_grid(
    plot_grid(
      plot_grid(comb_sig_prop, comb_girl_prop, ncol = 2, nrow = 1, labels = c("E", "F",  "")),
      plot_grid(overlap_genes_plot, NULL, overlap_TF_plot, NULL, overlap_Disease_plot, 
      nrow = 1, labels = c("G", "H", ""), rel_widths = c(1, 0.05, 0.7, 0.05, 1)),
      nrow = 2),
  topGeneSets, ncol = 2, rel_widths = c(1.5, 1), labels = c("", "I")),
  nrow = 3, rel_heights = c(1, 1, 2)
)
dev.off()


## Pre vs Post-menopause
tab.ship_joint <- left_join(tab.ship.pre_main, tab.ship.post_main, by = "ILMN_Gene", suffix = c(".pre", ".post")) %>%
  left_join(tab.ship.main, by = "ILMN_Gene") %>%
  left_join(select(tab.array, GeneSymbol, Dimorphic), by = c("ILMN_Gene" = "GeneSymbol")) %>%
  as_tibble() %>%
  mutate(Dimorphic = ifelse(is.na(Dimorphic), "Not in Meta-analysis", Dimorphic))

group_by(tab.ship_joint, Dimorphic) %>%
  summarize(pre_post = cor(logFC.pre, logFC.post, use = "complete.obs"),
            main_post = cor(logFC.post, logFC, use = "complete.obs"),
            main_pre = cor(logFC.pre, logFC, use = "complete.obs"))


png("figures/genes/ship_pre_post.png", width = 1800, height = 1500, res = 300)
tab.ship_joint %>%
  mutate(Effect = ifelse(Dimorphic == "Dimorphic", ifelse(logFC > 0, "Higher Men", "Higher Women"), "No differences")) %>%
  ggplot(aes(x = logFC.pre, y = logFC.post, color = Effect)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("Higher Men" = "blue", "Higher Women" = "maroon", "No differences" = "black")) +
  xlab("Pre-menopause logFC") +
  ylab("Post-menopause logFC") +
  facet_grid(Dimorphic ~.) +
  stat_cor(data = tab.ship_joint, inherit.aes = FALSE, aes(x = logFC.pre, y = logFC.post), method = "pearson", label.x = -0.5, label.y = 0.5)
dev.off()

write.table(tab.ship_joint %>%
select(ILMN_Gene, chr, ends_with(".pre"), ends_with(".post"), Dimorphic), 
file = "results/genes/SHIP_menopause_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Comparison with GTEx
## GTEx ####
load("results/preprocessFiles/gexp_gtex_sva_sex.Rdata")
ciber_cells_gtex <- colnames(colData(gtex_cells_filt))[390:411]
ciber_cells_callrate <- sapply(ciber_cells_gtex, function(x) mean(colData(gtex_cells_filt)[, x] > 0.01))
ciber_cells <- ciber_cells_gtex[ciber_cells_callrate > 0.05]


mod_gtex <- model.matrix(formula(paste("~ sex + AGE + ischemic + gtex.smrin +", paste(ciber_cells, collapse = "+"), "+",
                                              paste(paste0("SV", 1:56), collapse = "+"))),
                                data = colData(gtex_cells_filt))

blood_gtex_noY <- gtex_cells_filt[as.character(seqnames(gtex_cells_filt)) %in% paste0("chr", c(1:22, "X")), ]

gtex_voom <- voom(assay(blood_gtex_noY), mod_gtex)
lm.gtex <- lmFit(gtex_voom, mod_gtex) %>%
  eBayes()
tab.gtex <- topTable(lm.gtex, coef = 2, n = Inf, confint = TRUE)
tab.gtex$Symbol <- rowData(blood_gtex_noY)[rownames(tab.gtex), "gene_name"]
tab.gtex$chromosome <- as.character(seqnames(blood_gtex_noY[rownames(tab.gtex)]))
save(tab.gtex, file = "results/genes/gtex_analysis.Rdata")
write.table(tab.gtex %>%
  mutate(ENSG = rownames(tab.gtex)) %>%
  select(ENSG, Symbol, chromosome, logFC,CI.L, CI.R, P.Value, adj.P.Val), 
  file = "results/genes/gtex_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# ### Run analysis stratified by ventilator
# load("results/preprocessFiles/gexp_gtex_vent_stratified.Rdata")

# ciber_cells_filt <- ciber_cells[ciber_cells != "Dendritic.cells.resting"]

# mod_gtex_vent <- model.matrix(formula(paste("~ sex + AGE + ischemic + gtex.smrin +", paste(ciber_cells_filt, collapse = "+"), "+",
#                                               paste(paste0("SV", 1:40), collapse = "+"))),
#                                 data = colData(gtex_vent))

# vent_gtex_noY <- gtex_vent[as.character(seqnames(gtex_vent)) %in% paste0("chr", c(1:22, "X")), ]

# vent_voom <- voom(assay(vent_gtex_noY), mod_gtex_vent)
# lm.vent <- lmFit(vent_voom, mod_gtex_vent) %>%
#   eBayes()
# tab.vent <- topTable(lm.vent, coef = 2, n = Inf, confint = TRUE)
# tab.vent$Symbol <- rowData(vent_gtex_noY)[rownames(tab.vent), "gene_name"]
# tab.vent$chromosome <- as.character(seqnames(vent_gtex_noY[rownames(tab.vent)]))
# tab.vent$ENSG <- rownames(tab.vent)


# mod_gtex_novent <- model.matrix(formula(paste("~ sex + AGE + ischemic + gtex.smrin +", paste(ciber_cells_filt, collapse = "+"), "+",
#                                               paste(paste0("SV", 1:30), collapse = "+"))),
#                                 data = colData(gtex_novent))

# novent_gtex_noY <- gtex_novent[as.character(seqnames(gtex_novent)) %in% paste0("chr", c(1:22, "X")), ]

# novent_voom <- voom(assay(novent_gtex_noY), mod_gtex_novent)
# lm.novent <- lmFit(novent_voom, mod_gtex_novent) %>%
#   eBayes()
# tab.novent <- topTable(lm.novent, coef = 2, n = Inf, confint = TRUE)
# tab.novent$Symbol <- rowData(novent_gtex_noY)[rownames(tab.novent), "gene_name"]
# tab.novent$chromosome <- as.character(seqnames(novent_gtex_noY[rownames(tab.novent)]))
# tab.novent$ENSG <- rownames(tab.novent)
# save(tab.vent, tab.novent, file = "results/genes/gtex_analysis_ventilator.Rdata")

# tab.strat <- left_join(select(tab.vent, ENSG, Symbol, logFC, P.Value), select(tab.novent, ENSG, logFC, P.Value), 
#   by = "ENSG", suffix = c(".vent", ".novent")) %>%
#   as_tibble()


# tab.adult_strats <- left_join(tab.array, tab.strat, by = c("GeneSymbol" = "Symbol"), suffix = c(".array", ".gtex")) %>%
#   as_tibble()

# group_by(tab.adult_strats, Dimorphic, logFC < 0) %>%
#   summarize(vent = cor(logFC.com, logFC.vent, use = "complete.obs"),
#   no_vent = cor(logFC.com, logFC.novent, use = "complete.obs"),
#   vent_no_vent = cor(logFC.vent, logFC.novent, use = "complete.obs"))


# png("figures/genes/array_vs_vent.png", width = 1500, height = 1500, res = 300)
# tab.adult_strats %>% mutate(Effect = ifelse(Dimorphic == "Dimorphic", ifelse(logFC.com > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
#   ggplot( aes(x = logFC.com, y = logFC.vent, color = Effect)) +
#   geom_point() +
#   xlab("logFC Meta-analysis") +
#   ylab("logFC GTEx") +
#   scale_color_manual(values = c("Higher Boys" = "blue", "Higher Girls" = "maroon", "No differences" = "black")) +
#   facet_grid(Dimorphic ~.) +
#   theme_bw() +
#   stat_cor(method = "pearson", label.x = -1.2, label.y = -0.7, 
#   color = "black") +
#   theme(legend.position = "none")

# dev.off()

# png("figures/genes/array_vs_novent.png", width = 1500, height = 1500, res = 300)
# tab.adult_strats %>% mutate(Effect = ifelse(Dimorphic == "Dimorphic", ifelse(logFC.com > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
#   ggplot( aes(x = logFC.com, y = logFC.novent, color = Effect)) +
#   geom_point() +
#   xlab("logFC Meta-analysis") +
#   ylab("logFC GTEx") +
#   scale_color_manual(values = c("Higher Boys" = "blue", "Higher Girls" = "maroon", "No differences" = "black")) +
#   facet_grid(Dimorphic ~.) +
#   theme_bw() +
#   stat_cor(method = "pearson", label.x = -1.2, label.y = -0.7, 
#   color = "black") +
#   theme(legend.position = "none")
# dev.off()

# png("figures/genes/vent_vs_novent.png", width = 1500, height = 1500, res = 300)
# tab.strat %>% mutate(Dim_vent = ifelse(P.Value.vent  < 1e-4, "Ventilator Signif.", "Ventilator Not Signif."),
# Dim_novent = ifelse(P.Value.novent  < 1e-4, "Non-Ventilator Signif.", "Non-Ventilator Not Signif.")) %>%
#   ggplot( aes(x = logFC.vent, y = logFC.novent)) +
#   geom_point() +
#   xlab("logFC Ventilator") +
#   ylab("logFC Non-Ventilator") +
#   facet_grid(Dim_vent ~Dim_novent) +
  
#   theme_bw() +
#   stat_cor(method = "pearson", label.x = -10, label.y = 5, 
#   color = "black") +
#   theme(legend.position = "none")
# dev.off()



# df <- tibble(gexp = gtex_voom$E["ENSG00000206047.2",], 
#              sex = gtex_cells_filt$sex,
#              age = gtex_cells_filt$AGE,
#               B.cells.naive = gtex_cells_filt$B.cells.naive,
#               ischem = gtex_cells_filt$gtex.smtsisch,
#               hardy = factor(gtex_cells_filt$DTHHRDY))

# summary(lm(gexp ~ sex + age, df))
# summary(lm(gexp ~ sex + age + B.cells.naive, df))
# summary(lm(gexp ~ sex*B.cells.naive, df))
# summary(lm(gexp ~ sex*ischem + age + B.cells.naive, df))
# summary(lm(gexp ~ sex*hardy + age + B.cells.naive, df))

# summary(lm(gexp ~ sex + age + B.cells.naive, df, subset = ischem < 0))
# summary(lm(gexp ~ sex + age + B.cells.naive, df, subset = ischem > 0))
# summary(lm(gexp ~ sex + age + B.cells.naive, df, subset = hardy != 0))
# summary(lm(gexp ~ sex + age + B.cells.naive, df, subset = hardy == 0))

## CIBERSORT-based cell type proportion
ciber_mat <- acomp(data.matrix(colData(gtex_cells_filt)[, ciber_cells]))
ciber_d <- apply(ciber_mat, 2, function(x) min(x[x > 0], na.rm = TRUE))
ciber_mat_zero <- zeroreplace(ciber_mat, d = ciber_d)


gtex_clr <- clr(ciber_mat_zero) %>%
  as.data.frame() %>%
  mutate(age = gtex_cells_filt$AGE,
  sex = gtex_cells_filt$sex,
  rin = gtex_cells_filt$gtex.smrin,
  ischem = gtex_cells_filt$ischemic, 
  hardy = factor(gtex_cells_filt$DTHHRDY))  %>%
  mutate(ventilator = ifelse(hardy == 0, "Ventilator", "Other")) 


clr_gtex_tests <- lapply(ciber_cells, function(cell) {
  lm_mod <- robustbase::lmrob(as.formula(paste(cell, "~ sex + age + rin + ischem")), data = gtex_clr) %>%
    summary() 
  tibble(Cell = cell, 
         Estimate = coef(lm_mod)[2, 1], 
         Estimate_exp = exp(coef(lm_mod)[2, 1]),
         Pvalue = coef(lm_mod)[2, 4])
}) %>% Reduce(rbind, .)

clr_gtex_ischem_tests <- lapply(ciber_cells, function(cell) {
  lm_mod <- robustbase::lmrob(as.formula(paste(cell, "~ ischem + age + rin + hardy")), data = gtex_clr) %>%
    summary() 
  tibble(Cell = cell, 
         Estimate = coef(lm_mod)[2, 1], 
         Estimate_exp = exp(coef(lm_mod)[2, 1]),
         Pvalue = coef(lm_mod)[2, 4])
}) %>% Reduce(rbind, .)

clr_gtex_ventilator_tests <- lapply(ciber_cells, function(cell) {
  lm_mod <- robustbase::lmrob(as.formula(paste(cell, "~ ventilator + ischem + age + rin")), data = gtex_clr) %>%
    summary() 
  tibble(Cell = cell, 
         Estimate = coef(lm_mod)[2, 1], 
         Estimate_exp = exp(coef(lm_mod)[2, 1]),
         Pvalue = coef(lm_mod)[2, 4])
}) %>% Reduce(rbind, .)



clr_gtex_matrix_tests <- lapply(ciber_cells, function(cell) {
  lm_mod <- lm(as.formula(paste(cell, "~ sex + age + ventilator + ischem + rin")), data = gtex_clr) %>%
    summary() 
  coef(lm_mod)[-1, 4]
}) %>% Reduce(rbind, .)
rownames(clr_gtex_matrix_tests) <- ciber_cells
colnames(clr_gtex_matrix_tests) <- c("Sex", "Age", "Ventilator", "Ischemic Time", "RIN")


clr_df_tests <- clr_gtex_matrix_tests %>%
  as.data.frame() %>%
  rownames_to_column("Cell_type") %>%
  pivot_longer(
    -Cell_type,
    names_to = "Covariate",
    values_to = "Pvalue"
  ) %>%
  mutate(
    Cell_type = gsub("\\.", " ", Cell_type),
    minus_log10_p = -log10(Pvalue),
    label = sprintf("%.2e", Pvalue),
    Covariates = factor(Covariate, levels = c("Sex", "Age", "Ventilator", "Ischemic Time", "RIN"))
  )

png("figures/genes/gtex_cibersort_celltype_pvalues.png", width = 2000, height = 1500, res = 300)
ggplot(clr_df_tests, aes(x = Covariates, y = Cell_type, fill = minus_log10_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = expression(-log[10](P))
  ) +
  labs(
    title = "Covariates vs Cell type proportions",
    x = "Covariates",
    y = "Cell types"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# clr_gtex_hardy_tests <- lapply(ciber_cells, function(cell) {
#   lm_mod <- lm(as.formula(paste(cell, "~ hardy + ischem + age + rin")), data = gtex_clr) %>%
#     summary() 
# })


# clr_gtex_tests2 <- lapply(ciber_cells, function(cell) {
#   lm_mod <- lm(as.formula(paste(cell, "~ sex + age")), data = gtex_clr) %>%
#     summary() 
#   tibble(Cell = cell, 
#          Estimate = coef(lm_mod)[2, 1], 
#          Estimate_exp = exp(coef(lm_mod)[2, 1]),
#          Pvalue = coef(lm_mod)[2, 4])
# }) %>% Reduce(rbind, .)

## Check association between sex and covariates
sex_covar_tests <- lapply(c("AGE", "gtex.smrin", "ischemic"), function(covar) {
  lm_mod <- lm(as.formula(paste(covar, "~ sex")), data = colData(gtex_cells_filt)) %>%
    summary()
  tibble(Covariate = covar,
         Estimate = coef(lm_mod)[2, 1],
         Pvalue = coef(lm_mod)[2, 4])
}) %>% Reduce(rbind, .)


# gtex_voom_scaled <- gtex_voom
# gtex_voom_scaled$E <- t(scale(t(gtex_voom_scaled$E)))
# lm.gtex_scaled <- lmFit(gtex_voom_scaled, mod_gtex) %>%
#   eBayes()
# tab.gtex_scaled <- topTable(lm.gtex_scaled, coef = 2, n = Inf, confint = TRUE)
# tab.gtex_scaled$Symbol <- rowData(blood_gtex_noY)[rownames(tab.gtex_scaled), "gene_name"]
# tab.gtex_scaled$chromosome <- as.character(seqnames(blood_gtex_noY[rownames(tab.gtex_scaled)]))

tab.adults <- left_join(tab.array, tab.gtex, by = c("GeneSymbol" = "Symbol"), suffix = c(".array", ".gtex")) %>%
  as_tibble()


adult_comb_plot <- tab.adults %>% mutate(Effect = ifelse(Dimorphic == "Dimorphic", ifelse(logFC.com > 0, "Higher Men", "Higher Women"), "No differences")) %>%
  ggplot( aes(x = logFC.com, y = logFC.gtex, color = Effect)) +
  geom_point() +
  xlab("logFC Meta-analysis") +
  ylab("logFC GTEx") +
  facet_grid(Dimorphic ~ .) +
  theme_bw() +
  stat_cor(aes(group = Effect, color = Effect), method = "pearson", label.x = -1.2, label.y = c(0.5, 0.4, 0.5),
    p.accuracy = 1e-7, show.legend = FALSE) +
    scale_color_manual(values = c("Higher Men" = "blue", "Higher Women" = "maroon", "No differences" = "black"))

png("figures/genes/array_vs_gtex.png", width = 1500, height = 1200, res = 300)
adult_comb_plot
dev.off()

group_by(tab.adults, Dimorphic, logFC.com < 0) %>%
  summarize(cor = cor(logFC.com, logFC.gtex, use = "complete.obs"),
  cor2 = cor(logFC.array, EffectSize, use = "complete.obs"))

group_by(tab.adults, adj.P.Val.random < 0.01) %>%
    summarize(cor = cor(logFC.com, logFC.gtex, use = "complete.obs"),
    cor2 = cor(logFC.array, logFC.gtex, use = "complete.obs"))

table(tab.adults$Dimorphic, sign(tab.adults$logFC.com) == sign(tab.adults$logFC.gtex))


png("figures/Genes/EffectSizesChild.png", width = 2000, height = 1200, res = 300)
child_comb %>%
  mutate(Effect = ifelse(dimorphic == "Dimorphic", ifelse(logFC.helix > 0, "Higher Boys", "Higher Girls"), "No differences")) %>%
  ggplot(aes(x = logFC.helix, y = logFC.GenR, color = Effect)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point() +
  stat_cor(data = child_comb, inherit.aes = FALSE, aes(x = logFC.helix, y = logFC.GenR), method = "pearson", label.x = -1.45, label.y = 1.85) +
  theme_bw() +
  scale_color_manual(values = c("Higher Boys" = "blue", "Higher Girls" = "maroon", "No differences" = "black")) +
  facet_grid(dimorphic ~.) +
  xlab("Helix logFC") +
  ylab("GenR logFC") +
  coord_cartesian(ylim = c(-4, 2), xlim = c(-1.5, 1.5), clip = "off") 
dev.off()




## Children vs adults
## Sup Fig 2
effects_comb <- inner_join(gene_cohort_summary %>%
                        mutate(Dimorphic = dimorphic) %>%
                        dplyr::select(Symbol, Dimorphic, chr, logFC, P.Value),
                      tab.array %>%
                        mutate(Symbol = GeneSymbol) %>%
                      dplyr::select(Symbol, Dimorphic, logFC.com, P.Value.com),
                      by = "Symbol", suffix = c(".Children", ".Adults") )  %>%
    mutate(Dimorphic_comb = ifelse(Dimorphic.Children == "Dimorphic", 
                                   ifelse(Dimorphic.Adults == "Dimorphic", "Shared", "Children"),
                                   ifelse(Dimorphic.Adults == "Dimorphic", "Adults", "None")))

png("figures/genes/adult_vs_children_effect.png", width = 1500, height = 1500, res = 300)
ggplot(effects_comb, 
  aes(x = logFC, y = logFC.com, color = Dimorphic_comb)) +
  geom_point(alpha = 0.7) +
  xlab("logFC Children") +
  ylab("logFC Adults") +
  scale_color_manual(values = c("forestgreen", "orange", "black", "#CFCD8B")) +
  facet_wrap(~Dimorphic_comb) +
  theme_bw() +
  stat_cor(method = "pearson", label.x = -1.2, label.y = -0.7, 
  color = "black") +
  theme(legend.position = "none")
dev.off()

