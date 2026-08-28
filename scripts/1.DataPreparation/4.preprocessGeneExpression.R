#'##############################################################################
#' Preprocess Gene expression data
#' docker run -it -v $PWD:$PWD -w $PWD dimorphic_r:1.4 R
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(limma)
# library(GEOquery)
library(SummarizedExperiment)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
library(sva)
# library(affy)


#'##############################################################################
#' HELIX
#'##############################################################################
load("results/preprocessFiles/Expression_SE_raw_filtered.Rdata")

# Compute required SVs
mat <- assay(se.filt)
pd <- colData(se.filt)

## Add covariates to the model
pd$cohort <- droplevels(pd$cohort)
mod <- model.matrix(~ cohort + e3_sex + age_sample_years +
                      NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                      Mono_6, data = pd)
mod0 <- model.matrix(~ cohort + age_sample_years +
                      NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                      Mono_6, data = pd)
n.sv <- num.sv(mat, mod, method="leek") + 1
print(n.sv)
# 3 SVs

sv.obj <- sva(mat, mod, mod0 = NULL, n.sv = n.sv)

## Add SVs to SE
helix_sv <- sv.obj$sv
colnames(helix_sv) <- paste0("SV", seq_len(ncol(helix_sv)))
rownames(helix_sv) <- colnames(se.filt)

colData(se.filt) <- cbind(colData(se.filt), helix_sv)

# Save
save(se.filt, file = "results/preprocessFiles/Expression_SE_raw_filtered_SV.RData")

# ## Prepare gene annotaion ####
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# all_genes <- genes(txdb, "gene_id")
# ENTREZID2SYMBOL <- select(org.Hs.eg.db, mcols(all_genes)$gene_id, c("ENTREZID", "SYMBOL"))
# mcols(all_genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
# names(all_genes) <- mcols(all_genes)$SYMBOL

# #'##############################################################################
# #' GSE36382
# #'##############################################################################

# ## Download data ####
# gse36382_geo <- getGEO("GSE36382")[[1]]
# gse36382_probes <- read.ilmn("data/GSE36382_non_normalized.txt", expr = "Sample",
#                              probeid = "ID_REF")
# gse36382_exprs <- neqc(gse36382_probes)

# ## Preprocess data ####
# ## Combine probes for the same gene
# gse36382_gene <- avereps(gse36382_exprs, ID = fData(gse36382_geo)[rownames(gse36382_exprs), "Symbol"])

# ## Discard probes not annotated to a Symbol
# gse36382_gene <- gse36382_gene[rownames(gse36382_gene) != "", ]

# ## Create SE and remove genes not present in TxDb
# gse36382_se <- SummarizedExperiment(gse36382_gene)
# gse36382_se <- gse36382_se[rownames(gse36382_se) %in% names(all_genes), ]
# rowRanges(gse36382_se) <- all_genes[rownames(gse36382_se)]

# ### Infer sex from chrY expression ####
# chry_exprs <- assay(gse36382_se[as.character(seqnames(gse36382_se)) == "chrY",])$E
# pc_chry <- prcomp(t(chry_exprs))
# gse36382_se$sex <- ifelse(pc_chry$x[, 1] > 0, "male", "female")

# ### Compute SVA ####
# gse36382_mat <- assay(gse36382_se) %>% data.matrix()
# Y.r_gse36382 <- t(resid(lm(t(gse36382_mat) ~ sex, 
#   data = colData(gse36382_se))))

# mod_gse36382 <- model.matrix(~ sex, data = colData(gse36382_se))
# mod0_gse36382 <- model.matrix(~ 1, data = colData(gse36382_se))

# n.sv_gse36382 <- num.sv(assay(gse36382_se), mod_gse36382, method="leek")
# print(n.sv_gse36382)
# # 0
# save(gse36382_se, file = "results/preprocessFiles/gexp_gse36382_sex.Rdata")

# #'##############################################################################
# #' GSE33828
# #'##############################################################################

# ### Download data ####
# gse33828_geo <- GEOquery::getGEO("GSE33828")[[1]]

# gse33828_probes <- read_delim("data/GSE33828_RS3-ExpressionMatrix-881-raw-data-Ordered.txt.gz")
# gse33828_mat <- data.matrix(gse33828_probes[, -1])
# rownames(gse33828_mat) <- gse33828_probes$PROBE_ID
# gse33828_exprs <- neqc(gse33828_mat, 
#                        detection.p = matrix(runif(nrow(gse33828_mat)*ncol(gse33828_mat), 1e-10, 1e-3),
#                                                   nrow(gse33828_mat), ncol(gse33828_mat)))
# ## Preprocess data ####
# ## Combine probes for the same gene
# gse33828_gene <- avereps(gse33828_exprs, ID = fData(gse33828_geo)[rownames(gse33828_exprs), "Symbol"])

# ## Discard probes not annotated to a Symbol
# gse33828_gene <- gse33828_gene[!is.na(rownames(gse33828_gene)), ]

# ## Create SE and remove genes not present in TxDb
# gse33828_se <- SummarizedExperiment(gse33828_gene)
# gse33828_se <- gse33828_se[rownames(gse33828_se) %in% names(all_genes), ]
# rowRanges(gse33828_se) <- all_genes[rownames(gse33828_se)]

# gse33828_se$age <- as.numeric(gse33828_geo$`age (y):ch1`)
# gse33828_se$sex <- gse33828_geo$`gender:ch1`


# ## Remove sample with negative age
# gse33828_se <- gse33828_se[,gse33828_se$age > 0]

# ### Compute SVA ####
# mod_gse33828 <- model.matrix(~ sex + age, data = colData(gse33828_se))
# mod0_gse33828 <- model.matrix(~ age, data = colData(gse33828_se))

# n.sv_gse33828 <- num.sv(assay(gse33828_se), mod_gse33828, method="leek")
# print(n.sv_gse33828)
# # 0
# save(gse33828_se, file = "results/preprocessFiles/gexp_gse33828_sex.Rdata")

# #'##############################################################################
# #' GSE43488
# #'##############################################################################

# ### Download data
# gse43488_geo <- getGEO("GSE43488")[[1]]
# gse43488_gene <- avereps(exprs(gse43488_geo), ID = fData(gse43488_geo)$"Gene Symbol")

# gse43488_se <- SummarizedExperiment(gse43488_gene)
# gse43488_se$age <- as.numeric(gse43488_geo$`age at sample (months):ch1`)/12
# gse43488_se$sex <- gse43488_geo$`gender:ch1`
# gse43488_se$t1d_diagnosis <- gse43488_geo$`time from t1d diagnosis (months):ch1`
# gse43488_se$seroconversion <- gse43488_geo$`time from seroconversion (months):ch1`

# ### Select controls
# gse43488_controls <- gse43488_se[, gse43488_se$seroconversion == "no seroconversion" & gse43488_se$t1d_diagnosis == "no T1D diagnosis"]
# ## Discard samples > 10yo
# gse43488_controls <- gse43488_controls[, gse43488_controls$age < 10]

# ### Compute SVA ####
# mod_gse43488 <- model.matrix(~ sex + age, data = colData(gse43488_controls))
# mod0_gse43488 <- model.matrix(~ age, data = colData(gse43488_controls))

# n.sv_gse43488 <- num.sv(assay(gse43488_controls), mod_gse43488, method="leek") + 1
# sv.obj_gse43488 <- sva(assay(gse43488_controls), mod_gse43488, mod0 = mod0_gse43488, n.sv = n.sv_gse43488)
# svs_gse43488 <- sv.obj_gse43488$sv
# rownames(svs_gse43488) <- rownames(mod_gse43488)
# colnames(svs_gse43488) <- paste0("SV", seq_len(n.sv_gse43488))

# colData(gse43488_controls) <- cbind(colData(gse43488_controls), svs_gse43488)
# save(gse43488_controls, file = "results/preprocessFiles/gexp_GSE43488_sex.Rdata")

# GTEx ####
## Whole blood samples
load("results/preprocessFiles/gexp_gtex_EUR_sex.Rdata")

gtex_cells <- read_tsv("results/CIBERSORT/CIBERSORTx_Job10_Results.txt", name_repair = "universal") %>%
    mutate(ID = gsub("\\.1", "", Mixture)) 

colData_gtex <- left_join(colData(rse_gtex_filt) %>% as.data.frame(), gtex_cells, by = c("gtex.sampid" = "ID")) %>%
  as_tibble()


gtex_cells <- SummarizedExperiment(assay(rse_gtex_filt), colData = colData_gtex, rowRanges = rowRanges(rse_gtex_filt))
gtex_cells$ischemic <- ifelse(gtex_cells$gtex.smtsisch < 0, 0, gtex_cells$gtex.smtsisch)
ciber_cells <- colnames(colData(gtex_cells))[390:410]

##
gtex_cells_filt <- gtex_cells[, !is.na(gtex_cells$ischemic)]
mod1_gtex <- model.matrix(formula(paste("~sex + AGE + ischemic + gtex.smrin + ", paste(ciber_cells, collapse = " + "))), colData(gtex_cells_filt))
mod0_gtex <- model.matrix(formula(paste("~ AGE + ischemic + gtex.smrin +", paste(ciber_cells, collapse = " + "))), colData(gtex_cells_filt))
svseq <- svaseq(assay(gtex_cells_filt), mod1_gtex, mod0_gtex)$sv
## n.sv = 56
colnames(svseq) <- paste0("SV", 1:56)
rownames(svseq) <- colnames(gtex_cells_filt)
colData(gtex_cells_filt) <- cbind(colData(gtex_cells_filt), svseq)
save(gtex_cells_filt, file = "results/preprocessFiles/gexp_gtex_sva_sex.Rdata")

## Compute SVAs for ventilator and non-ventilator samples
gtex_cells_hardy <- gtex_cells_filt[, !is.na(gtex_cells_filt$DTHHRDY)]
gtex_vent <- gtex_cells_hardy[, gtex_cells_hardy$DTHHRDY == 0]

ciber_cells_filt <- ciber_cells[ciber_cells != "Dendritic.cells.resting"]

mod1_gtex_vent <- model.matrix(formula(paste("~sex + AGE + ischemic + gtex.smrin + ", paste(ciber_cells_filt, collapse = " + "))), colData(gtex_vent))
mod0_gtex_vent <- model.matrix(formula(paste("~ AGE + ischemic + gtex.smrin +", paste(ciber_cells_filt, collapse = " + "))), colData(gtex_vent))
svseq_vent <- svaseq(assay(gtex_vent), mod1_gtex_vent, mod0_gtex_vent)$sv
## n.sv = 40
colnames(svseq_vent) <- paste0("SV", 1:40)
rownames(svseq_vent) <- colnames(gtex_vent)
colData(gtex_vent) <- cbind(colData(gtex_vent), svseq_vent)


gtex_novent <- gtex_cells_hardy[, gtex_cells_hardy$DTHHRDY != 0]
mod1_gtex_novent <- model.matrix(formula(paste("~sex + AGE + ischemic + gtex.smrin + ", paste(ciber_cells_filt, collapse = " + "))), colData(gtex_novent))
mod0_gtex_novent <- model.matrix(formula(paste("~ AGE + ischemic + gtex.smrin +", paste(ciber_cells_filt, collapse = " + "))), colData(gtex_novent))
svseq_novent <- svaseq(assay(gtex_novent), mod1_gtex_novent, mod0_gtex_novent)$sv
## n.sv = 30
colnames(svseq_novent) <- paste0("SV", 1:30)
rownames(svseq_novent) <- colnames(gtex_novent)
colData(gtex_novent) <- cbind(colData(gtex_novent), svseq_novent)


save(gtex_vent, gtex_novent, file = "results/preprocessFiles/gexp_gtex_vent_stratified.Rdata")

# Other tissues ####
load("results/preprocessFiles/gtex_tissues.Rdata")


runSVAseq <- function(rse){

  rse_filt <- rse[, !is.na(rse$gtex.smtsisch)]
  mod1_gtex <- model.matrix(~sex + AGE + gtex.smtsisch + gtex.smrin, colData(rse_filt))
  mod0_gtex <- model.matrix(~AGE + gtex.smtsisch + gtex.smrin, colData(rse_filt))
  
  svseq <- svaseq(assay(rse_filt), mod1_gtex, mod0_gtex)$sv
  colnames(svseq) <- paste0("SV", seq_len(ncol(svseq)))
  rownames(svseq) <- colnames(rse_filt)
  colData(rse_filt) <- cbind(colData(rse_filt), svseq)
  rse_filt
}
gtex_tissues_sva <- lapply(gtex_tissues_proc, runSVAseq)

save(gtex_tissues_sva, file = "results/preprocessFiles/gtex_tissues_sva_sex.Rdata")

