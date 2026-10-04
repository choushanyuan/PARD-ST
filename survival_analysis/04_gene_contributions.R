# Step 4: Gene Set Contribution Analysis (Univariate Cox)
library(survival)
library(ggplot2)
library(dplyr)

base_dir <- "/path/to/survival_project"
result_dir <- file.path(base_dir, "result")
figure_dir <- file.path(base_dir, "figure")

load(file.path(base_dir, "P1_P3_P4_overlap_pesudotime_genes_new.RData")) # down_genes (Up), up_genes (Down)

run_uni_cox_genes <- function(expr_mat, clin_df, genes) {
  available_genes <- intersect(genes, rownames(expr_mat))
  res_list <- list()
  for (g in available_genes) {
    expr_vals <- as.numeric(expr_mat[g, ])
    fit_data <- data.frame(time = clin_df$BCR_time, status = clin_df$BCR_status, expr = expr_vals)
    fit_data <- na.omit(fit_data)
    if (nrow(fit_data) < 10 || sum(fit_data$status == 1) < 3) next
    
    fit <- tryCatch({ coxph(Surv(time, status) ~ expr, data = fit_data) }, error = function(e) { NULL })
    if (!is.null(fit)) {
      s <- summary(fit)
      res_list[[g]] <- data.frame(
        Gene = g, HR = s$coefficients[1, "exp(coef)"], Pvalue = s$coefficients[1, "Pr(>|z|)"],
        Lower_95 = s$conf.int[1, "lower .95"], Upper_95 = s$conf.int[1, "upper .95"]
      )
    }
  }
  return(bind_rows(res_list))
}

plot_univariate_forest <- function(df, title) {
  df_sorted <- df[order(df$HR, decreasing = TRUE), ]
  df_sorted$Gene <- factor(df_sorted$Gene, levels = rev(df_sorted$Gene))
  
  # Format P-value label
  df_sorted$p_label <- ifelse(df_sorted$Pvalue < 0.001, "P < 0.001", sprintf("P = %.3f", df_sorted$Pvalue))
  
  x_max <- max(df_sorted$Upper_95, na.rm=TRUE)
  x_text <- x_max + (x_max - min(df_sorted$Lower_95, na.rm=TRUE)) * 0.05
  
  p <- ggplot(df_sorted, aes(x = HR, y = Gene, xmin = Lower_95, xmax = Upper_95)) +
    geom_point(color = "#3182BD", size = 1.2) +
    geom_errorbarh(height = 0.2, linewidth = 0.35, color = "#3182BD") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.35) +
    geom_text(aes(x = x_text, label = p_label), hjust = 0, size = 2.5, color = "black") +
    labs(title = title, x = "Hazard Ratio (95% CI)", y = "Gene") +
    theme_classic(base_size = 7.5, base_family = "sans") +
    coord_cartesian(clip = "off") +
    theme(
      axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8.5),
      axis.text.y = element_text(color = "black", size = 5),
      plot.margin = margin(t=10, r=50, b=10, l=10)
    )
  return(p)
}

for (expr_type in c("FPKM", "TPM")) {
  cat("Running Gene Contributions for", expr_type, "...\n")
  load(file.path(result_dir, expr_type, "cohort_clinical_with_scores.RData"))
  
  df_ch <- cohort_data$Changhai
  df_ch$BCR_status <- ifelse(df_ch$BCR == "YES", 1, ifelse(df_ch$BCR == "NO", 0, NA))
  df_ch$BCR_time <- as.numeric(gsub("m", "", df_ch$`Month to surgery/BCR`))
  df_ch_cap <- df_ch[which(df_ch$BCR_time <= 70), ]
  exp_ch_cap <- expr_ch[, rownames(df_ch_cap)]
  
  df_tcga <- cohort_data$TCGA
  df_tcga$BCR_status <- df_tcga$PFI
  df_tcga$BCR_time <- as.numeric(df_tcga$PFI.time) / 30.4
  
  # Standardize both to 15 chars for a guaranteed match ("TCGA-XX-XXXX-01")
  common_tcga_samples <- intersect(rownames(df_tcga), substr(colnames(expr_tcga), 1, 15))
  
  df_tcga <- df_tcga[common_tcga_samples, ]
  
  # Align expression matrix by finding the exact column indices that match the 15-char prefixes
  match_idx <- match(common_tcga_samples, substr(colnames(expr_tcga), 1, 15))
  exp_tcga_aligned <- expr_tcga[, match_idx, drop=FALSE]
  
  # Up-regulated
  uni_ch_up <- run_uni_cox_genes(exp_ch_cap, df_ch_cap, down_genes)
  write.csv(uni_ch_up, file.path(result_dir, expr_type, "Changhai_univariate_cox_up_genes.csv"), row.names = FALSE)
  if (nrow(uni_ch_up) > 0) {
    p <- plot_univariate_forest(uni_ch_up, "Changhai: Pseudotime Up-regulated Genes Univariate Cox")
    ggsave(file.path(figure_dir, expr_type, "Changhai_Up_Genes_Forest.pdf"), p, width = 6, height = 3 + 0.1 * nrow(uni_ch_up), limitsize=FALSE)
  }
  
  uni_tcga_up <- run_uni_cox_genes(exp_tcga_aligned, df_tcga, down_genes)
  if (nrow(uni_tcga_up) > 0) {
    write.csv(uni_tcga_up, file.path(result_dir, expr_type, "TCGA_univariate_cox_up_genes.csv"), row.names = FALSE)
    p <- plot_univariate_forest(uni_tcga_up, "TCGA: Pseudotime Up-regulated Genes Univariate Cox")
    ggsave(file.path(figure_dir, expr_type, "TCGA_Up_Genes_Forest.pdf"), p, width = 6, height = 3 + 0.1 * nrow(uni_tcga_up), limitsize=FALSE)
  } else {
    write.csv(data.frame(Gene=character(), HR=numeric(), Pvalue=numeric(), Lower_95=numeric(), Upper_95=numeric()), 
              file.path(result_dir, expr_type, "TCGA_univariate_cox_up_genes.csv"), row.names = FALSE)
  }
  
  # Down-regulated
  uni_ch_down <- run_uni_cox_genes(exp_ch_cap, df_ch_cap, up_genes)
  write.csv(uni_ch_down, file.path(result_dir, expr_type, "Changhai_univariate_cox_down_genes.csv"), row.names = FALSE)
  if (nrow(uni_ch_down) > 0) {
    p <- plot_univariate_forest(uni_ch_down, "Changhai: Pseudotime Down-regulated Genes Univariate Cox")
    ggsave(file.path(figure_dir, expr_type, "Changhai_Down_Genes_Forest.pdf"), p, width = 6, height = 3 + 0.1 * nrow(uni_ch_down), limitsize=FALSE)
  }
  
  uni_tcga_down <- run_uni_cox_genes(exp_tcga_aligned, df_tcga, up_genes)
  if (nrow(uni_tcga_down) > 0) {
    write.csv(uni_tcga_down, file.path(result_dir, expr_type, "TCGA_univariate_cox_down_genes.csv"), row.names = FALSE)
    p <- plot_univariate_forest(uni_tcga_down, "TCGA: Pseudotime Down-regulated Genes Univariate Cox")
    ggsave(file.path(figure_dir, expr_type, "TCGA_Down_Genes_Forest.pdf"), p, width = 6, height = 3 + 0.1 * nrow(uni_tcga_down), limitsize=FALSE)
  } else {
    write.csv(data.frame(Gene=character(), HR=numeric(), Pvalue=numeric(), Lower_95=numeric(), Upper_95=numeric()), 
              file.path(result_dir, expr_type, "TCGA_univariate_cox_down_genes.csv"), row.names = FALSE)
  }
}
cat("Gene contributions analysis completed successfully!\n")
