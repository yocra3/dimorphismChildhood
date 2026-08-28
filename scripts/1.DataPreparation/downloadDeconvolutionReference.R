#'##############################################################################
#' Download dataset for deconvolution reference
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(GEOquery)
library(SummarizedExperiment)
library(edgeR)

GSE107011_geo <- getGEO("GSE107011")[[1]]
counts <- read_tsv("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE107011&format=file&file=GSE107011_raw_counts_GRCh38.p13_NCBI.tsv.gz")

GSE107011_se <- SummarizedExperiment(assay = as.matrix(counts[, -1]), 
                               rowData = DataFrame(gene_id = counts$GeneID),
                               colData = DataFrame(GSE107011_geo))
annot <- read_tsv("data/Human.GRCh38.p13.annot.tsv")

## Mapping between cell types and main groups
cell_map <- c(
  "Naive CD4 T cells" = "CD4_T_cells",
  "Terminal effector CD4 T cells" = "CD4_T_cells",
  "Th1 cells" = "CD4_T_cells",
  "Th1/Th17 cells" = "CD4_T_cells",
  "Th17 cells" = "CD4_T_cells",
  "Th2 cells" = "CD4_T_cells",
  "Follicular helper T cells" = "CD4_T_cells",
  "T regulatory cells" = "CD4_T_cells",

  "Naive CD8 T cells" = "CD8_T_cells",
  "Central memory CD8 T cell" = "CD8_T_cells",
  "Effector memory CD8 T cells" = "CD8_T_cells",
  "Terminal effector CD8 T cells" = "CD8_T_cells",

  "Natural killer cells" = "NK_cells",
  "MAIT cells" = "NK_cells",
  "Vd2 gd T cells" = "NK_cells",
  "Non-Vd2 gd T cells" = "NK_cells",

  "Naive B cells" = "B_cells",
  "Exhausted B cells" = "B_cells",
  "Non-switched memory B cells" = "B_cells",
  "Switched memory B cells" = "B_cells",
  "Plasmablasts" = "B_cells",

  "Classical monocytes" = "Monocytes",
  "Intermediate monocytes" = "Monocytes",
  "Non classical monocytes" = "Monocytes",

  "Low-density neutrophils" = "Granulocytes",
  "Low-density basophils" = "Granulocytes",

  "Myeloid dendritic cells" = NA,
  "Plasmacytoid dendritic cells" = NA,
  "Progenitor cells" = NA,
  "PBMCs" = NA
)
GSE107011_se$cell_type <- cell_map[GSE107011_se$cell.type.ch1]

## Select only main cell types
GSE107011_filt <- GSE107011_se[, !is.na(GSE107011_se$cell_type)]

## Compute CPM
GSE107011_dge <- DGEList(counts = assay(GSE107011_filt), genes = rowData(GSE107011_filt)$gene_name)
GSE107011_dge <- calcNormFactors(GSE107011_dge)
GSE107011_cpm <- cpm(GSE107011_dge, normalized.lib.sizes = TRUE)

keep <- filterByExpr(assay(GSE107011_filt), group = GSE107011_filt$cell_type)

## Define matrix
gene_names <- annot$Symbol
names(gene_names) <- annot$GeneID

signature_matrix <- GSE107011_cpm
rownames(signature_matrix) <- gene_names[as.character(rowData(GSE107011_filt)$gene_id)]
colnames(signature_matrix) <- GSE107011_filt$cell_type

signature_matrix_filt <- signature_matrix[keep, ]
write.table(cbind(data.frame(GeneID = rownames(signature_matrix_filt)), as.data.frame(signature_matrix_filt)), 
    file = "results/CIBERSORT/DeconvolutionReference.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cell_types <- unique(GSE107011_filt$cell_type)
pheno_mat <- sapply(cell_types, function(cell) {
    as.numeric(cell == colnames(signature_matrix_filt))
}) %>% t()
pheno_mat[pheno_mat == 0] <- 2
colnames(pheno_mat) <- colnames(signature_matrix_filt)
write.table(cbind(data.frame(Cell = rownames(pheno_mat)), as.data.frame(pheno_mat)), 
    file = "results/CIBERSORT/DeconvolutionPhenotype.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


## Use all cell types
cell_map2 <- c(
  "Naive CD4 T cells" = "CD4_T_naive",
  "Th1 cells" = "CD4_T_effector",
  "Th2 cells" = "CD4_T_effector",
  "Th17 cells" = "CD4_T_effector",
  "Th1/Th17 cells" = "CD4_T_effector",
  "Follicular helper T cells" = "CD4_T_effector",
  "T regulatory cells" = "CD4_Treg",
  "Terminal effector CD4 T cells" = "CD4_T_effector",

  "Naive CD8 T cells" = "CD8_T_naive",
  "Central memory CD8 T cell" = "CD8_T_memory",
  "Effector memory CD8 T cells" = "CD8_T_memory",
  "Terminal effector CD8 T cells" = "CD8_T_effector",

  "Natural killer cells" = "NK_cells",
  "MAIT cells" = "MAIT",
  "Vd2 gd T cells" = "Gamma_delta_T",
  "Non-Vd2 gd T cells" = "Gamma_delta_T",

  "Naive B cells" = "B_naive",
  "Non-switched memory B cells" = "B_memory",
  "Switched memory B cells" = "B_memory",
  "Exhausted B cells" = "B_memory",
  "Plasmablasts" = "Plasmablasts",

  "Classical monocytes" = "Classical_monocytes",
  "Intermediate monocytes" = "Intermediate_monocytes",
  "Non classical monocytes" = "Nonclassical_monocytes",

  "Myeloid dendritic cells" = "mDC",
  "Plasmacytoid dendritic cells" = "pDC",

  "Low-density neutrophils" = "Granulocytes",
  "Low-density basophils" = "Granulocytes",

  "PBMCs" = NA,
  "Progenitor cells" = NA
)

GSE107011_se$cell_type2 <- cell_map2[GSE107011_se$cell.type.ch1]

## Select rpesent cell types
GSE107011_filt2 <- GSE107011_se[, !is.na(GSE107011_se$cell_type2)]

## Compute CPM
GSE107011_dge2 <- DGEList(counts = assay(GSE107011_filt2), genes = rowData(GSE107011_filt2)$gene_name)
GSE107011_dge2 <- calcNormFactors(GSE107011_dge2)
GSE107011_cpm2 <- cpm(GSE107011_dge2, normalized.lib.sizes = TRUE)

keep2 <- filterByExpr(assay(GSE107011_filt2), group = GSE107011_filt2$cell_type2)

## Define matrix
signature_matrix2 <- GSE107011_cpm2
rownames(signature_matrix2) <- gene_names[as.character(rowData(GSE107011_filt2)$gene_id)]
colnames(signature_matrix2) <- GSE107011_filt2$cell_type2



## Define matrix
signature_matrix_filt2 <- signature_matrix2[keep2, ]
colnames(signature_matrix_filt2) <- GSE107011_filt2$cell_type2
write.table(cbind(data.frame(GeneID = rownames(signature_matrix_filt2)), as.data.frame(signature_matrix_filt2)), 
    file = "results/CIBERSORT/DeconvolutionReference_allcells.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cell_types2 <- unique(GSE107011_filt2$cell_type2)
pheno_mat2 <- sapply(cell_types2, function(cell) {
    as.numeric(cell == colnames(signature_matrix_filt2))
}) %>% t()
pheno_mat2[pheno_mat2 == 0] <- 2
colnames(pheno_mat2) <- colnames(signature_matrix_filt2)
write.table(cbind(data.frame(Cell = rownames(pheno_mat2)), as.data.frame(pheno_mat2)), 
    file = "results/CIBERSORT/DeconvolutionPhenotype_allcells.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
