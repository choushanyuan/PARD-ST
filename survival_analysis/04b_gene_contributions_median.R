# Re-draw per-gene univariate Cox forest plots using a Median High/Low split.
# Inputs/outputs mirror code/04_gene_contributions.R, but each gene is grouped
# by its own median expression instead of being modeled as a continuous variable.
# Expected location: <project_root>/survival_analysis/04b_gene_contributions_median.R

library(survival)
library(ggplot2)
library(dplyr)

# Locate the project root from this script's own path, so the script is
# portable and free of hard-coded (non-ASCII) paths.
find_base_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg_idx <- grep("^--file=", args)
  if (length(file_arg_idx) > 0) {
    script_path <- sub("^--file=", "", args[file_arg_idx[1]])
    return(dirname(dirname(normalizePath(script_path, winslash = "/"))))
  }
  # RStudio: use the active source document
  ok <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
  if (ok) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(dirname(normalizePath(ctx$path, winslash = "/"))))
    }
  }
  stop("Could not locate the project root. Run this file from <root>/survival_analysis/.")
}

base_dir <- find_base_dir()
result_dir <- file.path(base_dir, "result")
figure_dir <- file.path(base_dir, "figure")

load(file.path(base_dir, "P1_P3_P4_overlap_pesudotime_genes_new.RData")) # down_genes, up_genes

# Median High/Low univariate Cox for every gene
run_uni_cox_genes_median <- function(expr_mat, clin_df, genes) {
  available_genes <- intersect(genes, rownames(expr_mat))
  res_list <- list()
  for (g in available_genes) {
    expr_vals <- as.numeric(expr_mat[g, ])
    fit_data <- data.frame(time = clin_df$BCR_time,
                           status = clin_df$BCR_status,
                           expr = expr_vals)
    fit_data <- na.omit(fit_data)
    if (nrow(fit_data) < 10 || sum(fit_data$status == 1) < 3) next

    med_val <- median(fit_data$expr)
    if (!is.finite(med_val)) next

    fit_data$group <- factor(ifelse(fit_data$expr >= med_val, "High", "Low"),
                             levels = c("Low", "High"))
    if (nlevels(droplevels(fit_data$group)) < 2) next

    fit <- tryCatch({
      coxph(Surv(time, status) ~ group, data = fit_data)
    }, error = function(e) NULL)
    if (!is.null(fit)) {
      s <- summary(fit)
      res_list[[g]] <- data.frame(
        Gene = g,
        Median_cutoff = med_val,
        N_high = sum(fit_data$group == "High"),
        N_low = sum(fit_data$group == "Low"),
        HR = s$coefficients["groupHigh", "exp(coef)"],
        Pvalue = s$coefficients["groupHigh", "Pr(>|z|)"],
        Lower_95 = s$conf.int["groupHigh", "lower .95"],
        Upper_95 = s$conf.int["groupHigh", "upper .95"]
      )
    }
  }
  return(bind_rows(res_list))
}

plot_univariate_forest_median <- function(df, title) {
  df_sorted <- df[order(df$HR, decreasing = TRUE), ]
  df_sorted$Gene <- factor(df_sorted$Gene, levels = rev(df_sorted$Gene))

  df_sorted$p_label <- ifelse(df_sorted$Pvalue < 0.001,
                              "P < 0.001",
                              sprintf("P = %.3f", df_sorted$Pvalue))

  x_max <- max(df_sorted$Upper_95, na.rm = TRUE)
  x_text <- x_max + (x_max - min(df_sorted$Lower_95, na.rm = TRUE)) * 0.05

  p <- ggplot(df_sorted, aes(x = HR, y = Gene, xmin = Lower_95, xmax = Upper_95)) +
    geom_point(color = "#3182BD", size = 1.2) +
    geom_errorbarh(height = 0.2, linewidth = 0.35, color = "#3182BD") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.35) +
    geom_text(aes(x = x_text, label = p_label), hjust = 0, size = 2.5, color = "black") +
    labs(title = title, x = "Hazard Ratio (95% CI)", y = "Gene") +
    theme_classic(base_size = 7.5, base_family = "sans") +
    coord_cartesian(clip = "off") +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8.5),
      axis.text.y = element_text(color = "black", size = 5),
      plot.margin = margin(t = 10, r = 50, b = 10, l = 10)
    )
  return(p)
}

for (expr_type in c("FPKM", "TPM")) {
  cat("Running median-cut gene contributions for", expr_type, "...\n")
  load(file.path(result_dir, expr_type, "cohort_clinical_with_scores.RData"))

  # Changhai
  df_ch <- cohort_data$Changhai
  df_ch$BCR_status <- ifelse(df_ch$BCR == "YES", 1, ifelse(df_ch$BCR == "NO", 0, NA))
  df_ch$BCR_time <- as.numeric(gsub("m", "", df_ch$`Month to surgery/BCR`))
  df_ch_cap <- df_ch[which(df_ch$BCR_time <= 70), ]
  exp_ch_cap <- expr_ch[, rownames(df_ch_cap)]

  # TCGA
  df_tcga <- cohort_data$TCGA
  df_tcga$BCR_status <- df_tcga$PFI
  df_tcga$BCR_time <- as.numeric(df_tcga$PFI.time) / 30.4
  common_tcga_samples <- intersect(rownames(df_tcga), substr(colnames(expr_tcga), 1, 15))
  df_tcga <- df_tcga[common_tcga_samples, ]
  match_idx <- match(common_tcga_samples, substr(colnames(expr_tcga), 1, 15))
  exp_tcga_aligned <- expr_tcga[, match_idx, drop = FALSE]

  # Up-regulated genes (pseudotime-up genes = down_genes in the RData)
  uni_ch_up <- run_uni_cox_genes_median(exp_ch_cap, df_ch_cap, down_genes)
  write.csv(uni_ch_up,
            file.path(result_dir, expr_type, "Changhai_univariate_cox_up_genes_median.csv"),
            row.names = FALSE)
  if (nrow(uni_ch_up) > 0) {
    p <- plot_univariate_forest_median(
      uni_ch_up,
      "Changhai: Pseudotime Up-regulated Genes Univariate Cox (Median High/Low)")
    ggsave(file.path(figure_dir, expr_type, "Changhai_Up_Genes_Forest_Median.pdf"), p,
           # cairo_pdf() embeds the font subset -> PDF stays editable; the
           # ggsave default pdf() device does not embed fonts.
           device = cairo_pdf,
           width = 6, height = 3 + 0.1 * nrow(uni_ch_up), limitsize = FALSE)
  }

  uni_tcga_up <- run_uni_cox_genes_median(exp_tcga_aligned, df_tcga, down_genes)
  write.csv(uni_tcga_up,
            file.path(result_dir, expr_type, "TCGA_univariate_cox_up_genes_median.csv"),
            row.names = FALSE)
  if (nrow(uni_tcga_up) > 0) {
    p <- plot_univariate_forest_median(
      uni_tcga_up,
      "TCGA: Pseudotime Up-regulated Genes Univariate Cox (Median High/Low)")
    ggsave(file.path(figure_dir, expr_type, "TCGA_Up_Genes_Forest_Median.pdf"), p,
           # cairo_pdf() embeds the font subset -> PDF stays editable; the
           # ggsave default pdf() device does not embed fonts.
           device = cairo_pdf,
           width = 6, height = 3 + 0.1 * nrow(uni_tcga_up), limitsize = FALSE)
  }

  # Down-regulated genes (pseudotime-down genes = up_genes in the RData)
  uni_ch_down <- run_uni_cox_genes_median(exp_ch_cap, df_ch_cap, up_genes)
  write.csv(uni_ch_down,
            file.path(result_dir, expr_type, "Changhai_univariate_cox_down_genes_median.csv"),
            row.names = FALSE)
  if (nrow(uni_ch_down) > 0) {
    p <- plot_univariate_forest_median(
      uni_ch_down,
      "Changhai: Pseudotime Down-regulated Genes Univariate Cox (Median High/Low)")
    ggsave(file.path(figure_dir, expr_type, "Changhai_Down_Genes_Forest_Median.pdf"), p,
           # cairo_pdf() embeds the font subset -> PDF stays editable; the
           # ggsave default pdf() device does not embed fonts.
           device = cairo_pdf,
           width = 6, height = 3 + 0.1 * nrow(uni_ch_down), limitsize = FALSE)
  }

  uni_tcga_down <- run_uni_cox_genes_median(exp_tcga_aligned, df_tcga, up_genes)
  write.csv(uni_tcga_down,
            file.path(result_dir, expr_type, "TCGA_univariate_cox_down_genes_median.csv"),
            row.names = FALSE)
  if (nrow(uni_tcga_down) > 0) {
    p <- plot_univariate_forest_median(
      uni_tcga_down,
      "TCGA: Pseudotime Down-regulated Genes Univariate Cox (Median High/Low)")
    ggsave(file.path(figure_dir, expr_type, "TCGA_Down_Genes_Forest_Median.pdf"), p,
           # cairo_pdf() embeds the font subset -> PDF stays editable; the
           # ggsave default pdf() device does not embed fonts.
           device = cairo_pdf,
           width = 6, height = 3 + 0.1 * nrow(uni_tcga_down), limitsize = FALSE)
  }

  cat("Finished median-cut analysis for", expr_type, "\n")
}

cat("All median-cut gene forest plots completed successfully!\n")
