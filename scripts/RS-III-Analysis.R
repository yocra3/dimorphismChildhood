#'##############################################################################
#' Main model for sexual dimorphism of gene expression in RS-III
#'##############################################################################

## Libraries ####
library(limma)
library(sva)
library(openxlsx)

## Set working directory ####
setwd("~/2ExtraProjects/Carlos2/array/")

## Load expression and phenotype data ####
Data <- read.table(
  "RS3-881samples-21238probes-NORMALIZED-QN-L2T.txt",
  header = TRUE,
  check.names = FALSE
)

Phen <- read.table(
  "RS3-881samples-COVARIATES.txt",
  header = TRUE,
  check.names = FALSE
)

sampleQC <- fread("RS3-767samples-with-good-quality-SNP-ARRAYS.txt")[[1]]
Data <- Data[, c("PROBES", intersect(sampleQC, colnames(Data)))]
Phen <- Phen[Phen$ergoid %in% sampleQC,]
Data <- Data[, colnames(Data) != "RS3-8354"]

expr_raw <- Data

#'##############################################################################
#' Prepare expression matrix
#'##############################################################################

## If first column contains probe/gene IDs, move it to rownames ####
if (!all(sapply(expr_raw, is.numeric))) {
  id_col <- which(!sapply(expr_raw, is.numeric))[1]
  rownames(expr_raw) <- expr_raw[[id_col]]
  expr_raw <- expr_raw[, -id_col]
}

## Convert to matrix ####
expr_raw <- as.matrix(expr_raw)
mode(expr_raw) <- "numeric"


#'##############################################################################
#' Match expression columns to phenotype rows
#'##############################################################################

if (!any(colnames(expr_raw) %in% rownames(Phen))) {
  
  possible_id_cols <- names(Phen)[sapply(
    Phen,
    function(x) any(as.character(x) %in% colnames(expr_raw))
  )]
  
  if (length(possible_id_cols) == 0) {
    stop("No matching sample IDs found between expression columns and Phen.")
  }
  
  sample_col <- possible_id_cols[1]
  rownames(Phen) <- as.character(Phen[[sample_col]])
}

common_samples <- intersect(colnames(expr_raw), rownames(Phen))

expr_raw <- expr_raw[, common_samples]
Phen <- Phen[common_samples, ]

stopifnot(all(colnames(expr_raw) == rownames(Phen)))


#'##############################################################################
#' Phenotype formatting
#'##############################################################################

## Sex: female should be reference
Phen$sex <- factor(
  Phen$sex,
  levels = c(1, 0),
  labels = c("Female", "Male")
)

Phen$fastingstate <- factor(
  Phen$fastingstate,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

Phen$current_smoking <- factor(
  Phen$current_smoking,
  levels = c(0, 1),
  labels = c("Not-smoker", "Smoker")
)

Phen$age <- Phen$age

Phen$plateidamplification <- as.factor(Phen$plateidamplification)

#'##############################################################################
#' Gene/probe filtering
#'##############################################################################

## Remove Y chromosome probes/genes if chromosome annotation is in rownames
expr_raw <- expr_raw[
  !grepl("^Y-|^chrY|^Y$|ChrY|chromosomeY", rownames(expr_raw), ignore.case = TRUE),
]

## Remove zero-variance probes
expr_raw <- expr_raw[apply(expr_raw, 1, var, na.rm = TRUE) > 0, ]

## Remove rows with too much missingness
expr_raw <- expr_raw[rowMeans(is.na(expr_raw)) < 0.20, ]

## Impute remaining missing values by probe mean
if (anyNA(expr_raw)) {
  for (i in seq_len(nrow(expr_raw))) {
    expr_raw[i, is.na(expr_raw[i, ])] <- mean(expr_raw[i, ], na.rm = TRUE)
  }
}


#'##############################################################################
#' Create raw and scaled expression matrices
#'##############################################################################

## dataset_raw: normalized log2 microarray values
dataset_raw <- expr_raw

## dataset_scaled: each probe/gene mean 0, SD 1
dataset_scaled <- t(scale(t(dataset_raw)))

## Remove probes that became NA after scaling
keep_scaled <- apply(dataset_scaled, 1, function(x) all(is.finite(x)))

dataset_raw <- dataset_raw[keep_scaled, ]
dataset_scaled <- dataset_scaled[keep_scaled, ]


#'##############################################################################
#' Complete-case filtering for model covariates
#'##############################################################################

model_covariates <- c(
  "sex",
  "age",
  "nrlymphocytes",
  "nrmonocytes",
  "nrgranulocytes",
  "RQS",
  "plateidamplification"
)

## Check missingness before filtering
print(colSums(is.na(Phen[, model_covariates])))

keep_samples <- complete.cases(Phen[, model_covariates])

Phen <- Phen[keep_samples, ]
dataset_raw <- dataset_raw[, rownames(Phen)]
dataset_scaled <- dataset_scaled[, rownames(Phen)]

stopifnot(all(colnames(dataset_scaled) == rownames(Phen)))


#'##############################################################################
#' Design matrix without SVs
#'##############################################################################

design_base <- model.matrix(
  ~ sex +
    age +
    nrlymphocytes +
    nrmonocytes +
    nrgranulocytes +
    RQS +
    plateidamplification,
  data = Phen
)

design_null <- model.matrix(
  ~ age +
    nrlymphocytes +
    nrmonocytes +
    nrgranulocytes +
    RQS +
    plateidamplification,
  data = Phen
)

stopifnot(ncol(dataset_scaled) == nrow(design_base))


#'##############################################################################
#' Estimate surrogate variables using sva
#'##############################################################################

set.seed(123)

n_sv <- num.sv(
  dat = dataset_scaled,
  mod = design_base,
  method = "leek"
)

cat("Estimated number of SVs:", n_sv, "\n")

if (n_sv > 0) {
  
  svobj <- sva(
    dat = dataset_scaled,
    mod = design_base,
    mod0 = design_null,
    n.sv = n_sv
  )
  
  sv_df <- as.data.frame(svobj$sv)
  colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
  
  design <- cbind(design_base, as.matrix(sv_df))
  
} else {
  
  sv_df <- NULL
  design <- design_base
}


#'##############################################################################
#' Function to run limma and extract sex effect
#'##############################################################################

run_limma_sex <- function(expr_mat, design, output_prefix) {
  
  fit <- lmFit(expr_mat, design)
  fit <- eBayes(fit)
  
  sex_coef <- grep("^sexMale$", colnames(design), value = TRUE)
  
  if (length(sex_coef) != 1) {
    stop("Could not uniquely identify the sex coefficient. Check design column names.")
  }
  
  tab <- topTable(
    fit,
    coef = sex_coef,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P",
    confint = TRUE
  )
  
  out <- data.frame(
    GeneSymbol = rownames(tab),
    EffectSize = tab$logFC,
    SE = abs(tab$logFC / tab$t),
    PValue = tab$P.Value,
    FDR = tab$adj.P.Val,
    CI_low = tab$CI.L,
    CI_high = tab$CI.R,
    stringsAsFactors = FALSE
  )
  
  write.xlsx(
    out,
    file = paste0(output_prefix, ".xlsx"),
    overwrite = TRUE
  )
  
  return(out)
}


#'##############################################################################
#' Run models
#'##############################################################################

## Raw/log2-normalized expression:
## for interpretation, reporting and plotting
twas_raw <- run_limma_sex(
  expr_mat = dataset_raw,
  design = design,
  output_prefix = "BIOS_RSIII_TWAS_results_raw"
)

## Scaled expression:
## for meta-analysis
twas_scaled <- run_limma_sex(
  expr_mat = dataset_scaled,
  design = design,
  output_prefix = "BIOS_RSIII_TWAS_results_scaled"
)


lines <- readLines("HumanHT-12_V4_Illumina_annotation.txt")

start <- grep("^\\[Probes\\]", lines)

annotationIllumina <- read.delim(
  text = lines[(start + 1):length(lines)],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

## Raw results
twas_raw$ProbeID <- twas_raw$GeneSymbol

twas_raw <- merge(
  twas_raw,
  annotationIllumina,
  by.x = "ProbeID",
  by.y = "Probe_Id",
  all.x = TRUE,
  sort = FALSE
)

twas_raw <- twas_raw[, c(
  "ProbeID",
  "Symbol",
  "EffectSize",
  "SE",
  "PValue",
  "FDR",
  "CI_low",
  "CI_high"
)]

names(twas_raw)[2] <- "GeneSymbol"

## Scaled results
twas_scaled$ProbeID <- twas_scaled$GeneSymbol

twas_scaled <- merge(
  twas_scaled,
  annotationIllumina,
  by.x = "ProbeID",
  by.y = "Probe_Id",
  all.x = TRUE,
  sort = FALSE
)

twas_scaled <- twas_scaled[, c(
  "ProbeID",
  "Symbol",
  "EffectSize",
  "SE",
  "PValue",
  "FDR",
  "CI_low",
  "CI_high"
)]

names(twas_scaled)[2] <- "GeneSymbol"

## Overwrite the Excel files with annotated results
write.xlsx(
  twas_raw,
  "BIOS_RSIII_TWAS_results_raw.xlsx",
  overwrite = TRUE
)

write.xlsx(
  twas_scaled,
  "BIOS_RSIII_TWAS_results_scaled.xlsx",
  overwrite = TRUE
)

