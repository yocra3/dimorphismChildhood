#'##############################################################################
#' Preprocess data for replication
#'##############################################################################

## Load libraries ####
library(tidyverse)
library(limma)
library(GEOquery)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(sva)
library(affy)

## Prepare gene annotaion ####
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- genes(txdb, "gene_id")
ENTREZID2SYMBOL <- select(org.Hs.eg.db, mcols(all_genes)$gene_id, c("ENTREZID", "SYMBOL"))
mcols(all_genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
names(all_genes) <- mcols(all_genes)$SYMBOL

# GSE36382 ####
## Download data ####
gse36382_geo <- getGEO("GSE36382")[[1]]
gse36382_probes <- read.ilmn("data/GSE36382_non_normalized.txt", expr = "Sample",
                             probeid = "ID_REF")
gse36382_exprs <- neqc(gse36382_probes)

## Preprocess data ####
## Combine probes for the same gene
gse36382_gene <- avereps(gse36382_exprs, ID = fData(gse36382_geo)[rownames(gse36382_exprs), "Symbol"])

## Discard probes not annotated to a Symbol
gse36382_gene <- gse36382_gene[rownames(gse36382_gene) != "", ]

## Create SE and remove genes not present in TxDb
gse36382_se <- SummarizedExperiment(gse36382_gene)
gse36382_se <- gse36382_se[rownames(gse36382_se) %in% names(all_genes), ]
rowRanges(gse36382_se) <- all_genes[rownames(gse36382_se)]

### Infer sex from chrY expression ####
chry_exprs <- assay(gse36382_se[as.character(seqnames(gse36382_se)) == "chrY",])$E
pc_chry <- prcomp(t(chry_exprs))
gse36382_se$sex <- ifelse(pc_chry$x[, 1] > 0, "male", "female")

### Compute SVA ####
mod_gse36382 <- model.matrix(~ sex, data = colData(gse36382_se))
mod0_gse36382 <- model.matrix(~ 1, data = colData(gse36382_se))

n.sv_gse36382 <- num.sv(assay(gse36382_se), mod_gse36382, method="leek")
print(n.sv_gse36382)
# 0
save(gse36382_se, file = "results/preprocessFiles/gexp_gse36382_sex.Rdata")


## GSE33828 ####
### Download data ####
gse33828_geo <- GEOquery::getGEO("GSE33828")[[1]]

gse33828_probes <- read_delim("data/GSE33828_RS3-ExpressionMatrix-881-raw-data-Ordered.txt")
gse33828_mat <- data.matrix(gse33828_probes[, -1])
rownames(gse33828_mat) <- gse33828_probes$PROBE_ID
gse33828_exprs <- neqc(gse33828_mat, 
                       detection.p = matrix(runif(nrow(gse33828_mat)*ncol(gse33828_mat), 1e-10, 1e-3),
                                                  nrow(gse33828_mat), ncol(gse33828_mat)))
## Preprocess data ####
## Combine probes for the same gene
gse33828_gene <- avereps(gse33828_exprs, ID = fData(gse33828_geo)[rownames(gse33828_exprs), "Symbol"])

## Discard probes not annotated to a Symbol
gse33828_gene <- gse33828_gene[!is.na(rownames(gse33828_gene)), ]

## Create SE and remove genes not present in TxDb
gse33828_se <- SummarizedExperiment(gse33828_gene)
gse33828_se <- gse33828_se[rownames(gse33828_se) %in% names(all_genes), ]
rowRanges(gse33828_se) <- all_genes[rownames(gse33828_se)]

gse33828_se$age <- as.numeric(gse33828_geo$`age (y):ch1`)
gse33828_se$sex <- gse33828_geo$`gender:ch1`


## Remove sample with negative age
gse33828_se <- gse33828_se[,gse33828_se$age > 0]

### Compute SVA ####
mod_gse33828 <- model.matrix(~ sex + age, data = colData(gse33828_se))
mod0_gse33828 <- model.matrix(~ age, data = colData(gse33828_se))

n.sv_gse33828 <- num.sv(assay(gse33828_se), mod_gse33828, method="leek")
print(n.sv_gse33828)
# 0
save(gse33828_se, file = "results/preprocessFiles/gexp_gse33828_sex.Rdata")

## GSE43488 ####
### Download data
gse43488_geo <- getGEO("GSE43488")[[1]]
gse43488_gene <- avereps(exprs(gse43488_geo), ID = fData(gse43488_geo)$"Gene Symbol")

gse43488_se <- SummarizedExperiment(gse43488_gene)
gse43488_se$age <- as.numeric(gse43488_geo$`age at sample (months):ch1`)/12
gse43488_se$sex <- gse43488_geo$`gender:ch1`
gse43488_se$t1d_diagnosis <- gse43488_geo$`time from t1d diagnosis (months):ch1`
gse43488_se$seroconversion <- gse43488_geo$`time from seroconversion (months):ch1`

### Select controls
gse43488_controls <- gse43488_se[, gse43488_se$seroconversion == "no seroconversion" & gse43488_se$t1d_diagnosis == "no T1D diagnosis"]
## Discard samples > 12yo
gse43488_controls <- gse43488_controls[, gse43488_controls$age < 10]

### Compute SVA ####
mod_rep <- model.matrix(~ sex + age, data = colData(gse43488_controls))
mod_rep0 <- model.matrix(~ age, data = colData(gse43488_controls))

n.sv_rep <- num.sv(assay(gse43488_controls), mod_rep, method="leek")
print(n.sv_rep)
# 0
sv.obj_rep <- sva(assay(gse43488_controls), mod_rep, mod0 = mod_rep0, n.sv = n.sv_rep)
svs_rep <- sv.obj_rep$sv
rownames(svs_rep) <- rownames(mod_rep)
colnames(svs_rep) <- c("SV1")

colData(gse43488_controls) <- cbind(colData(gse43488_controls), svs_rep)
save(gse43488_controls, file = "results/preprocessFiles/gexp_sva_sex.Rdata")


save(gse43488_controls, file = "results/preprocessFiles/gexp_GSE43488_sex.Rdata")

# GTEx ####
### Download data ####
library(recount3)
library(edgeR)
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "BLOOD" 
)
rse_gtex <- create_rse(proj_info)
rse_gtex$sex <- ifelse(rse_gtex$gtex.sex == 1, "male", "female")
keep <- filterByExpr(assay(rse_gtex), group = rse_gtex$sex)
rse_gtex_filt <- rse_gtex[keep, ]

mod1_gtex <- model.matrix(~sex + gtex.age, colData(rse_gtex_filt) )
mod0_gtex <- model.matrix(~ gtex.age, colData(rse_gtex_filt))
svseq <- svaseq(assay(rse_gtex_filt), mod1_gtex, mod0_gtex)$sv
## n.sv = 37
colnames(svseq) <- paste0("SV", 1:37)
rownames(svseq) <- colnames(rse_gtex_filt)
colData(rse_gtex_filt) <- cbind(colData(rse_gtex_filt), svseq)
save(rse_gtex_filt, file = "results/preprocessFiles/gexp_gtex_sva_sex.Rdata")
