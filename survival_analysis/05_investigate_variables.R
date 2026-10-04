# 05_investigate_variables.R
# This script explores the clinical columns of the Changhai and TCGA cohorts, summarises how the TNM and margin information are recorded, and discusses how to set up the multivariable Cox regression

# ---- 1. Load the data ----
rdata_path <- "/path/to/survival_project/result/FPKM/cohort_clinical_with_scores.RData"
if (file.exists(rdata_path)) {
  load(rdata_path)
  cat("cohort_clinical_with_scores.RData loaded successfully.\n\n")
} else {
  stop("RData file not found: ", rdata_path)
}

ch_df <- cohort_data$Changhai
tcga_df <- cohort_data$TCGA

# ---- 2. Changhai cohort: column names and value summaries ----
cat("===================== 1. Changhai cohort =====================\n")
cat("Total samples: ", nrow(ch_df), "\n\n")

cat("[2.1 Clinical column names]\n")
ch_target_cols <- c("Pathology T Stage", "Pathology N stage", "Pathology M stage", "Surgical Margin")
for (col in ch_target_cols) {
  cat("  -", col, "\n")
}
cat("\n")

cat("[2.2 Clinical column values]\n")
for (col in ch_target_cols) {
  cat("\n---", col, "distribution ---\n")
  tab <- table(ch_df[[col]], useNA = "always")
  print(tab)
  cat("Non-missing fraction: ", round(sum(!is.na(ch_df[[col]])) / nrow(ch_df) * 100, 2), "%\n")
}
cat("\n")

# ---- 3. TCGA cohort: column names and value summaries ----
cat("===================== 2. TCGA cohort =====================\n")
cat("Total samples: ", nrow(tcga_df), "\n\n")

cat("[3.1 Clinical column names]\n")
tcga_target_cols <- c("pathologic_T", "pathologic_N", "residual_tumor")
for (col in tcga_target_cols) {
  cat("  -", col, "\n")
}
cat("\n")

cat("[3.2 Clinical column values]\n")
for (col in tcga_target_cols) {
  cat("\n---", col, "distribution ---\n")
  tab <- table(tcga_df[[col]], useNA = "always")
  print(tab)
  cat("Non-missing fraction: ", round(sum(!is.na(tcga_df[[col]])) / nrow(tcga_df) * 100, 2), "%\n")
}
cat("\n")

# ---- 4. How to restructure the multivariable Cox regression ----
cat("===================== 3. Restructuring the multivariable Cox regression =====================\n")
cat("Before fitting the multivariable Cox model these clinical variables have to be grouped and recoded sensibly:\n\n")

cat("A. Pathology T stage:\n")
cat("   - Changhai (Pathology T Stage): pT2 (3), pT2a (8), pT2b (7), pT2c (48), pT3a (41), pT3b (22), pT4 (4).\n")
cat("   - TCGA (pathologic_T): T2a (13), T2b (10), T2c (163), T3a (155), T3b (131), T4 (10), plus a few blanks or discrepancies.\n")
cat("   - Suggested treatment: to avoid subgroups that are too small and make Cox regression non-convergent, merge them into two groups (binary):\n")
cat("     * T2 (localised disease): T2/pT2/pT2a/pT2b/pT2c\n")
cat("     * T3_T4 (advanced disease): T3a/T3b/T4/pT3a/pT3b/pT4\n\n")

cat("B. Pathology N stage:\n")
cat("   - Changhai (Pathology N stage): pn0 (1), pN0 (94), pN1 (17), pNx (13), pNX (8).\n")
cat("   - TCGA (pathologic_N): N0 (341), N1 (76), blank (72).\n")
cat("   - Suggested treatment:\n")
cat("     * N0 (no lymph-node metastasis): N0/pN0/pn0\n")
cat("     * N1 (lymph-node metastasis): N1/pN1\n")
cat("     * NX/pNx/pNX/NA/blank: treated as NA (missing) and excluded from the regression\n\n")

cat("C. Pathology M stage - Changhai cohort only:\n")
cat("   - Changhai (Pathology M stage): M0 (11), pM0 (56), M1b (2), pM1a (12), pM1b (1), Mx/pMx/pMX/PMx (many unknown/missing).\n")
cat("     Including pM would cut the complete-case sample size from 112 to 78, which can give very low power or a non-convergent model.\n")
cat("     In radical-prostatectomy cohorts almost no patient has distant metastasis at surgery (i.e. M0), so M stage is usually not forced into the multivariable analysis;\n")
cat("     alternatively it could be dichotomised as M0 (M0/pM0/0) vs M1 (M1b/pM1a/pM1b), but with so few events this has to be handled cautiously (only 78 patients would remain for the multivariable Cox model).\n\n")

cat("D. Surgical margin:\n")
cat("   - Changhai (Surgical Margin): NO (81), YES (51), NA (1).\n")
cat("     * NO maps to 'Negative'\n")
cat("     * YES maps to 'Positive'\n")
cat("   - TCGA (residual_tumor): R0 (311), R1 (143), R2 (5), RX (15), blank (15).\n")
cat("     * R0 maps to 'Negative' (margin negative, no residual tumour)\n")
cat("     * R1 and R2 map to 'Positive' (microscopic/gross residual tumour)\n")
cat("     * RX/blank becomes NA (missing)\n\n")

cat("E. Model-fitting strategy for the multivariable Cox regression:\n")
cat("   The final multivariable model formula should be:\n")
cat("   1. Changhai cohort: Surv(time, status) ~ Expression + PSA + GleasonScore + TumorPurity + pT + pN + Margin + strata(Age_group)\n")
cat("      (Note: whether pM is included depends on whether the loss of sample size is acceptable. At this stage pT, pN and Margin are included as binary variables;\n")
cat("       a dichotomised pM can be tested later if needed, while watching the sample size.)\n")
cat("   2. TCGA cohort: Surv(time, status) ~ Expression + PSA + GleasonScore + TumorPurity + pT + pN + Margin + strata(Age_group)\n")
cat("   Note: pT, pN and Margin are converted into binary factors (with T2, N0 and Negative as the reference levels).\n")
