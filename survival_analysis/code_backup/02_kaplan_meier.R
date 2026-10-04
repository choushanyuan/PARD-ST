# Step 2: Kaplan-Meier Survival Analysis for FPKM and TPM
library(survival)
library(survminer)
library(ggplot2)
library(patchwork)
library(cowplot)

base_dir <- "/path/to/survival_project"
figure_dir <- file.path(base_dir, "figure")

palette_nature <- c("#D24B40", "#3182BD")

theme_km_nature <- function(base_size = 7.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
}

run_km_analysis <- function(df, time_col, status_col, score_col, cutoff_method, title_prefix, y_label) {
  df_clean <- df[, c(time_col, status_col, score_col)]
  colnames(df_clean) <- c("time", "status", "score")
  
  df_clean$time <- as.numeric(df_clean$time)
  df_clean$status <- as.numeric(df_clean$status)
  df_clean$score <- as.numeric(df_clean$score)
  
  if (grepl("Changhai", title_prefix)) {
    df_clean <- df_clean[which(df_clean$time <= 70), ]
  }
  
  df_clean <- df_clean[complete.cases(df_clean), ]
  if (nrow(df_clean) < 5 || sum(df_clean$status == 1) < 2) return(NULL)
  
  if (cutoff_method == "Median") {
    med_val <- median(df_clean$score)
    df_clean$group <- ifelse(df_clean$score >= med_val, "High", "Low")
  } else if (cutoff_method == "Q3_Upper") {
    q3_val <- quantile(df_clean$score, 0.75)
    df_clean$group <- ifelse(df_clean$score >= q3_val, "High", "Low")
  } else if (cutoff_method == "Q1_Lower") {
    q1_val <- quantile(df_clean$score, 0.25)
    df_clean$group <- ifelse(df_clean$score >= q1_val, "High", "Low")
  } else if (cutoff_method == "Q1_vs_Q3") {
    q1_val <- quantile(df_clean$score, 0.25)
    q3_val <- quantile(df_clean$score, 0.75)
    df_clean$group <- ifelse(df_clean$score >= q3_val, "High", 
                             ifelse(df_clean$score <= q1_val, "Low", NA))
    df_clean <- df_clean[!is.na(df_clean$group), ]
  } else if (cutoff_method == "Optimal") {
    tryCatch({
      cut <- surv_cutpoint(df_clean, time = "time", event = "status", variables = "score")
      df_clean$group <- ifelse(df_clean$score >= cut$cutpoint$cutpoint, "High", "Low")
    }, error = function(e) {
      med_val <- median(df_clean$score)
      df_clean$group <- ifelse(df_clean$score >= med_val, "High", "Low")
    })
  }
  
  # Calculate sample sizes
  n_high <- sum(df_clean$group == "High")
  n_low <- sum(df_clean$group == "Low")
  lbl_high <- sprintf("High (n=%d)", n_high)
  lbl_low <- sprintf("Low (n=%d)", n_low)
  
  df_clean$group <- factor(df_clean$group, levels = c("High", "Low"), labels = c(lbl_high, lbl_low))
  fit <- survfit(Surv(time, status) ~ group, data = df_clean)
  
  p_surv <- ggsurvplot(
    fit, 
    data = df_clean,
    pval = TRUE, 
    pval.method = TRUE, 
    pval.size = 2.5,
    size = 0.45,
    palette = palette_nature,
    legend.title = score_col,
    legend.labs = c(lbl_high, lbl_low),
    ggtheme = theme_km_nature(),
    ylab = y_label,
    xlab = "Time (months)",
    title = paste0(title_prefix, "
(", cutoff_method, ")")
  )
  
  return(p_surv$plot)
}

analytes <- c("HSPD1", "FASN", "PKP1", "Up_score", "Down_score")
cutoffs <- c("Median", "Q3_Upper", "Q1_Lower", "Q1_vs_Q3", "Optimal")

for (expr_type in c("FPKM", "TPM")) {
  cat("Running KM analysis for", expr_type, "...\n")
  load(file.path(base_dir, "result", expr_type, "cohort_clinical_with_scores.RData"))
  
  plots_ch <- list()
  df_ch <- cohort_data$Changhai
  df_ch$BCR_status <- ifelse(df_ch$BCR == "YES", 1, ifelse(df_ch$BCR == "NO", 0, NA))
  df_ch$BCR_time <- as.numeric(gsub("m", "", df_ch$`Month to surgery/BCR`))
  
  if (!all(analytes %in% colnames(df_ch))) {
    cat("Missing columns in Changhai:", setdiff(analytes, colnames(df_ch)), "\n")
  }
  
  # Ensure time and status are numeric before plotting
  df_ch$BCR_status <- as.numeric(df_ch$BCR_status)
  df_ch$BCR_time <- as.numeric(df_ch$BCR_time)
  
  for (a in analytes) {
    for (c in cutoffs) {
      p <- run_km_analysis(df_ch, "BCR_time", "BCR_status", a, c, paste0("Changhai: ", a), "BCR probability")
      if (!is.null(p)) plots_ch[[paste(a, c, sep="_")]] <- p
    }
  }
  
  plots_tcga <- list()
  df_tcga <- cohort_data$TCGA
  df_tcga$BCR_status <- df_tcga$PFI
  df_tcga$BCR_time <- as.numeric(df_tcga$PFI.time) / 30.4
  
  if (!all(analytes %in% colnames(df_tcga))) {
    cat("Missing columns in TCGA:", setdiff(analytes, colnames(df_tcga)), "\n")
  }
  
  # Ensure time and status are numeric before plotting
  df_tcga$BCR_status <- as.numeric(df_tcga$BCR_status)
  df_tcga$BCR_time <- as.numeric(df_tcga$BCR_time)
  
  for (a in analytes) {
    for (c in cutoffs) {
      p <- run_km_analysis(df_tcga, "BCR_time", "BCR_status", a, c, paste0("TCGA: ", a), "PFS probability")
      if (!is.null(p)) plots_tcga[[paste(a, c, sep="_")]] <- p
    }
  }
  
  dir.create(file.path(figure_dir, expr_type), showWarnings = FALSE, recursive = TRUE)
  pdf(file.path(figure_dir, expr_type, "Cohorts_Survival_KM.pdf"), width = 11, height = 11, family = "sans")
  
  if (length(plots_ch) > 0) {
    grid_ch <- wrap_plots(plots_ch, ncol = 5)
    print(grid_ch)
  }
  
  if (length(plots_tcga) > 0) {
    grid_tcga <- wrap_plots(plots_tcga, ncol = 5)
    print(grid_tcga)
  }
  
  dev.off()
}
cat("KM Survival Curves completed!\n")
