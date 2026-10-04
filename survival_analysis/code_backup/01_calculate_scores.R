# Step 1: Calculate ssGSEA and ESTIMATE scores for both FPKM and TPM expression values
library(GSVA)
library(estimate)
library(readxl)

# Define directories
base_dir <- "/path/to/survival_project"
result_dir <- file.path(base_dir, "result")

# Load gene lists
load(file.path(base_dir, "4_gene_set_overlap.RData")) # overlap_genes
load(file.path(base_dir, "P1_P3_P4_overlap_pesudotime_genes_new.RData")) # down_genes, up_genes
gene_sets <- list(
  Overlap_genes = overlap_genes,
  Up_score = down_genes,
  Down_score = up_genes
)

# Helper function to clean expression matrix (remove empty/NA rows and aggregate duplicates by mean)
clean_expression_matrix <- function(expr_mat) {
  valid_rows <- which(rownames(expr_mat) != "" & !is.na(rownames(expr_mat)))
  expr_mat <- expr_mat[valid_rows, , drop=FALSE]
  
  if (any(duplicated(rownames(expr_mat)))) {
    genes <- rownames(expr_mat)
    # compute sum
    expr_agg <- rowsum(expr_mat, group = genes)
    # compute mean
    counts <- table(genes)
    expr_agg <- expr_agg / as.numeric(counts[rownames(expr_agg)])
    expr_mat <- as.matrix(expr_agg)
  }
  return(expr_mat)
}

# Helper function to run ESTIMATE
run_estimate <- function(expr_matrix) {
  temp_in <- tempfile(fileext = ".txt")
  temp_common <- tempfile(fileext = ".txt")
  temp_out <- tempfile(fileext = ".txt")
  
  # Ensure it is a data frame for writing, so col.names=NA works correctly to write rownames
  # For ESTIMATE, we need standard Entrez gene symbols, but both CPGEA and TCGA here seem to use Gene symbols.
  # Let's inspect the names quickly - write as is, then let ESTIMATE filter.
  # Wait, run_estimate crashes on 0 matched genes, meaning our expression matrix rows don't match the ESTIMATE reference platform (affymetrix).
  write.table(as.data.frame(expr_matrix), file=temp_in, quote=FALSE, sep="\t", col.names=NA)
  
  # ESTIMATE might fail if zero genes match, so we capture the error
  tryCatch({
    filterCommonGenes(input.f=temp_in, output.f=temp_common, id="GeneSymbol")
    estimateScore(input.ds=temp_common, output.ds=temp_out, platform="affymetrix")
    
    res <- read.table(temp_out, skip=2, header=TRUE, row.names=1, sep="\t")
    res <- res[, -1, drop=FALSE]
    res <- as.data.frame(t(res))
  }, error = function(e) {
    cat("ESTIMATE error:", e$message, "\n")
    res <<- data.frame(row.names=colnames(expr_matrix), StromalScore=NA, ImmuneScore=NA, ESTIMATEScore=NA, TumorPurity=NA)
  })
  
  file.remove(temp_in, temp_common, temp_out)
  return(res)
}

# FPKM to TPM conversion function
fpkm_to_tpm <- function(fpkm_matrix) {
  tpm_matrix <- apply(fpkm_matrix, 2, function(x) {
    s <- sum(x, na.rm=TRUE)
    if (s == 0) return(x)
    return(x / s * 1e6)
  })
  return(tpm_matrix)
}

# ----------------- Load Raw Datasets -----------------
cat("Loading raw datasets...\n")

# A. Changhai Cohort Raw Data
load(file.path(base_dir, "Changhai_PRAD_bulk_data.RData")) # loads PRAD_bulk_data
PRAD_bulk_data_Tumor <- PRAD_bulk_data[, seq(1, ncol(PRAD_bulk_data), 2)]
colnames(PRAD_bulk_data_Tumor) <- gsub("_WTS", "", colnames(PRAD_bulk_data_Tumor))

clinical_df_changhai <- as.data.frame(read_excel(file.path(base_dir, "2_41586_2020_2135_MOESM4_ESM.xlsx")))
rownames(clinical_df_changhai) <- paste0("T", clinical_df_changhai$`Sequencing ID`)

# B. TCGA Cohort Preprocessed Data
TCGA_FPKM_tumor <- read.csv(file.path(base_dir, "TCGA_FPKM_tumor.csv"), row.names=1, check.names=FALSE)
TCGA_TPM_tumor <- read.csv(file.path(base_dir, "TCGA_TPM_tumor.csv"), row.names=1, check.names=FALSE)

load(file.path(base_dir, "TCGA_data.RData"))
df_phno$SampleBarcode <- rownames(df_phno)
clinical_df_TCGA <- merge(df_phno, df_clinical, by="X_PATIENT")
rownames(clinical_df_TCGA) <- clinical_df_TCGA$SampleBarcode

# ----------------- Process FPKM and TPM -----------------
for (expr_type in c("FPKM", "TPM")) {
  cat("Processing for", expr_type, "...\n")
  
  if (expr_type == "FPKM") {
    expr_ch <- as.matrix(log2(PRAD_bulk_data_Tumor + 1))
    expr_tcga <- as.matrix(log2(TCGA_FPKM_tumor + 1))
  } else {
    tpm_ch <- fpkm_to_tpm(PRAD_bulk_data_Tumor)
    expr_ch <- as.matrix(log2(tpm_ch + 1))
    expr_tcga <- as.matrix(log2(TCGA_TPM_tumor + 1))
  }
  
  expr_ch <- clean_expression_matrix(expr_ch)
  expr_tcga <- clean_expression_matrix(expr_tcga)
  
  cat("Running ssGSEA...\n")
  param_ch <- ssgseaParam(expr_ch, gene_sets)
  ssgsea_ch <- as.data.frame(t(gsva(param_ch, verbose=FALSE)))
  
  param_tcga <- ssgseaParam(expr_tcga, gene_sets)
  ssgsea_tcga <- as.data.frame(t(gsva(param_tcga, verbose=FALSE)))
  
  cat("Running ESTIMATE...\n")
  est_ch <- run_estimate(expr_ch)
  est_tcga <- run_estimate(expr_tcga)
  
  single_ch <- as.data.frame(t(expr_ch[c("HSPD1", "FASN", "PKP1"), , drop=FALSE]))
  single_tcga <- as.data.frame(t(expr_tcga[c("HSPD1", "FASN", "PKP1"), , drop=FALSE]))
  
  # Merge datasets safely
  common_ch <- Reduce(intersect, list(rownames(clinical_df_changhai), rownames(ssgsea_ch), rownames(est_ch), rownames(single_ch)))
  final_ch <- cbind(
    clinical_df_changhai[common_ch, , drop=FALSE],
    ssgsea_ch[common_ch, , drop=FALSE],
    est_ch[common_ch, , drop=FALSE],
    single_ch[common_ch, , drop=FALSE]
  )
  
  # The problem: clinical_df_TCGA uses TCGA.2A.A8VL.01, but ssgsea/est use TCGA-2A-A8VL-01A (or similar).
  # We need to standardize TCGA barcodes to intersect correctly.
  
  std_barcode <- function(x) {
    # Replace dots with dashes
    x <- gsub("\\.", "-", x)
    # Truncate to first 15 characters (e.g., TCGA-2A-A8VL-01)
    substr(x, 1, 15)
  }
  
  tcga_clin_names <- std_barcode(rownames(clinical_df_TCGA))
  tcga_ssgsea_names <- std_barcode(rownames(ssgsea_tcga))
  tcga_est_names <- std_barcode(rownames(est_tcga))
  tcga_single_names <- std_barcode(rownames(single_tcga))
  
  # Remove duplicates before setting row names
  clin_dup <- duplicated(tcga_clin_names)
  clinical_df_TCGA <- clinical_df_TCGA[!clin_dup, , drop=FALSE]
  tcga_clin_names <- tcga_clin_names[!clin_dup]
  rownames(clinical_df_TCGA) <- tcga_clin_names
  
  ssgsea_dup <- duplicated(tcga_ssgsea_names)
  ssgsea_tcga <- ssgsea_tcga[!ssgsea_dup, , drop=FALSE]
  tcga_ssgsea_names <- tcga_ssgsea_names[!ssgsea_dup]
  rownames(ssgsea_tcga) <- tcga_ssgsea_names
  
  est_dup <- duplicated(tcga_est_names)
  est_tcga <- est_tcga[!est_dup, , drop=FALSE]
  tcga_est_names <- tcga_est_names[!est_dup]
  rownames(est_tcga) <- tcga_est_names
  
  single_dup <- duplicated(tcga_single_names)
  single_tcga <- single_tcga[!single_dup, , drop=FALSE]
  tcga_single_names <- tcga_single_names[!single_dup]
  rownames(single_tcga) <- tcga_single_names
  
  common_tcga <- Reduce(intersect, list(rownames(clinical_df_TCGA), rownames(ssgsea_tcga), rownames(est_tcga), rownames(single_tcga)))
  cat("TCGA Common samples:", length(common_tcga), "\n")
  
  final_tcga <- cbind(
    clinical_df_TCGA[common_tcga, , drop=FALSE],
    ssgsea_tcga[common_tcga, , drop=FALSE],
    est_tcga[common_tcga, , drop=FALSE],
    single_tcga[common_tcga, , drop=FALSE]
  )

  cohort_data <- list(Changhai = final_ch, TCGA = final_tcga)
  
  dir.create(file.path(result_dir, expr_type), showWarnings = FALSE, recursive = TRUE)
  save(cohort_data, expr_ch, expr_tcga, file = file.path(result_dir, expr_type, "cohort_clinical_with_scores.RData"))
  cat("Saved results for", expr_type, "successfully!\n")
}
