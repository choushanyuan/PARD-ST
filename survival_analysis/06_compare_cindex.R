# 06_compare_cindex.R
# Compute the C-index of the Clinical-only vs the Clinical + Gene/Score model,
# use a likelihood-ratio test (LRT) to assess whether the improvement is significant, loop over all cutoffs,
# draw the comparison bar chart with significance labels and save one PDF per expr_type.

library(survival)
library(dplyr)
library(ggplot2)
library(grDevices)

# ------------------ 1. Setup paths ------------------
base_dir <- "/path/to/survival_project"
result_dir <- file.path(base_dir, "result")
figure_dir <- file.path(base_dir, "figure")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------ 2. Styles and Palette (Nature specifications) ------------------
color_map <- c(
  "Clinical only" = "#8E8E93",
  "Clinical + Score" = "#1F4E9E"
)

theme_nature_cindex <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      # rotate the x-axis text by 45 degrees, right- and top-aligned, so that it does not overlap the axis
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                 size = base_size - 0.5, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.position = "top",
      legend.margin = margin(b = -2, t = 0),
      # shrink the legend keys and their spacing
      legend.key.size = unit(3, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.spacing.x = unit(2, "pt"),
      strip.text = element_text(size = base_size + 0.5, face = "bold", colour = "black"),
      strip.background = element_blank(),
      plot.title = element_text(size = base_size + 1, face = "bold", colour = "black",
                                hjust = 0.5, margin = margin(b = 6)),
      panel.grid = element_blank(),
      panel.spacing = unit(10, "pt"),
      plot.margin = margin(t = 6, r = 8, b = 6, l = 6)
    )
}

# Significance stars
sig_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# ------------------ Stage / margin cleaning functions ------------------
clean_pT <- function(x) {
  x <- toupper(trimws(as.character(x)))
  g <- ifelse(grepl("T2", x), "T2",
              ifelse(grepl("T3|T4", x), "T3-T4", NA))
  factor(g, levels = c("T2", "T3-T4"))
}

clean_pN <- function(x) {
  x <- toupper(trimws(as.character(x)))
  g <- ifelse(grepl("N0", x), "N0",
              ifelse(grepl("N1", x), "N1", NA))
  factor(g, levels = c("N0", "N1"))
}

clean_margin_changhai <- function(x) {
  x <- toupper(trimws(as.character(x)))
  g <- ifelse(x == "NO", "Negative",
              ifelse(x == "YES", "Positive", NA))
  factor(g, levels = c("Negative", "Positive"))
}

clean_margin_tcga <- function(x) {
  x <- toupper(trimws(as.character(x)))
  g <- ifelse(x == "R0", "Negative",
              ifelse(x %in% c("R1", "R2"), "Positive", NA))
  factor(g, levels = c("Negative", "Positive"))
}

# ------------------ 3. Calculate C-indices ------------------
analytes <- c("HSPD1", "FASN", "PKP1", "Up_score", "Down_score")
cutoffs <- c("Median", "Q3_Upper", "Q1_Lower", "Q1_vs_Q3", "Optimal")
covars <- c("PSA", "GleasonScore", "TumorPurity", "pT", "pN", "Margin")

cindex_results_list <- list()

for (expr_type in c("FPKM", "TPM")) {
  cat("\n========== Calculating C-index for", expr_type, "==========\n")
  
  load(file.path(result_dir, expr_type, "cohort_clinical_with_scores.RData"))
  
  # A. Changhai cohort preprocessing
  df_ch <- cohort_data$Changhai
  df_ch <- df_ch[!is.na(df_ch$Age), ]
  df_ch$BCR_status <- ifelse(df_ch$BCR == "YES", 1, ifelse(df_ch$BCR == "NO", 0, NA))
  df_ch$BCR_time <- as.numeric(gsub("m", "", df_ch$`Month to surgery/BCR`))
  df_ch$GleasonScore <- as.numeric(as.character(df_ch$GS))
  df_ch$PSA <- log2(as.numeric(as.character(df_ch$PSA)) + 1)
  df_ch$TumorPurity <- df_ch$TumorPurity * 10
  df_ch$Age_group <- factor(ifelse(df_ch$Age >= 65, ">=65", "<65"), levels = c("<65", ">=65"))
  df_ch$pT     <- clean_pT(df_ch$`Pathology T Stage`)
  df_ch$pN     <- clean_pN(df_ch$`Pathology N stage`)
  df_ch$Margin <- clean_margin_changhai(df_ch$`Surgical Margin`)
  
  # B. TCGA cohort preprocessing
  df_tcga <- cohort_data$TCGA
  df_tcga <- df_tcga[!is.na(df_tcga$age_at_initial_pathologic_diagnosis), ]
  df_tcga$BCR_status <- df_tcga$PFI
  df_tcga$BCR_time <- as.numeric(df_tcga$PFI.time) / 30.4
  df_tcga$Age <- df_tcga$age_at_initial_pathologic_diagnosis
  df_tcga$GleasonScore <- as.numeric(as.character(df_tcga$gleason_score))
  df_tcga$PSA <- log2(as.numeric(as.character(df_tcga$psa_value)) + 1)
  df_tcga$TumorPurity <- df_tcga$TumorPurity * 10
  df_tcga$Age_group <- factor(ifelse(df_tcga$Age >= 65, ">=65", "<65"), levels = c("<65", ">=65"))
  df_tcga$pT     <- clean_pT(df_tcga$pathologic_T)
  df_tcga$pN     <- clean_pN(df_tcga$pathologic_N)
  df_tcga$Margin <- clean_margin_tcga(df_tcga$residual_tumor)

  get_cindex_comparison <- function(df, cohort_name, time_col, status_col, score_col, cutoff_method) {
    df_clean <- df[, c(time_col, status_col, score_col, covars, "Age_group")]
    colnames(df_clean) <- c("time", "status", "score", covars, "Age_group")
    df_clean$time   <- as.numeric(df_clean$time)
    df_clean$status <- as.numeric(df_clean$status)
    df_clean$score  <- as.numeric(df_clean$score)
    df_clean$Age_group <- factor(df_clean$Age_group, levels = c("<65", ">=65"))
    
    # convert only the non-factor covariates to numeric and keep pT/pN/Margin as factors with their reference levels
    for (cov in covars) {
      if (!is.factor(df_clean[[cov]])) df_clean[[cov]] <- as.numeric(df_clean[[cov]])
    }
    
    if (cohort_name == "Changhai") df_clean <- df_clean[which(df_clean$time <= 70), ]
    df_clean <- df_clean[complete.cases(df_clean), ]
    
    # drop the factor covariates that have only one level left after complete.cases, which would make coxph fail
    drop_covars <- c()
    for (cov in covars) {
      if (is.factor(df_clean[[cov]])) {
        df_clean[[cov]] <- droplevels(df_clean[[cov]])
        if (nlevels(df_clean[[cov]]) < 2) drop_covars <- c(drop_covars, cov)
      }
    }
    active_covars <- setdiff(covars, drop_covars)
    
    if (nrow(df_clean) < 10 || sum(df_clean$status == 1) < 3) return(NULL)
    
    # cutoff grouping for score
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
        cut <- survminer::surv_cutpoint(df_clean, time = "time", event = "status", variables = "score")
        df_clean$group <- ifelse(df_clean$score >= cut$cutpoint$cutpoint, "High", "Low")
      }, error = function(e) {
        med_val <- median(df_clean$score)
        df_clean$group <- ifelse(df_clean$score >= med_val, "High", "Low")
      })
    }
    df_clean$group <- factor(df_clean$group, levels = c("Low", "High"))
    
    # after grouping, check that both groups exist; a single-level Expression would make the combined model meaningless
    if (nlevels(droplevels(df_clean$group)) < 2) return(NULL)
    
    # 1. Clinical-only model
    formula_clin <- paste("Surv(time, status) ~ ", paste(active_covars, collapse = " + "),
                          "+ strata(Age_group)")
    fit_clin <- tryCatch(coxph(as.formula(formula_clin), data = df_clean), error = function(e) NULL)
    if (is.null(fit_clin)) return(NULL)
    
    # 2. Clinical + Score model
    df_clean$Expression <- df_clean$group
    formula_comb <- paste("Surv(time, status) ~ Expression +", paste(active_covars, collapse = " + "),
                          "+ strata(Age_group)")
    fit_comb <- tryCatch(coxph(as.formula(formula_comb), data = df_clean), error = function(e) NULL)
    if (is.null(fit_comb)) return(NULL)
    
    # Extract C-index & SE
    sum_clin <- concordance(fit_clin)
    sum_comb <- concordance(fit_comb)
    
    c_clin  <- sum_clin$concordance
    se_clin <- sqrt(sum_clin$var)
    c_comb  <- sum_comb$concordance
    se_comb <- sqrt(sum_comb$var)
    delta   <- c_comb - c_clin
    
    # likelihood-ratio test: nested-model comparison of whether adding the Score significantly improves the goodness of fit
    lrt_p <- tryCatch({
      lrt <- anova(fit_clin, fit_comb, test = "Chisq")
      lrt$`Pr(>|Chi|)`[2]
    }, error = function(e) NA)
    
    return(data.frame(
      Cohort = cohort_name,
      Analyte = score_col,
      Cutoff = cutoff_method,
      ExprType = expr_type,
      N_patients = nrow(df_clean),
      C_Clinical = c_clin,
      SE_Clinical = se_clin,
      C_Combined = c_comb,
      SE_Combined = se_comb,
      Delta_C = delta,
      LRT_Pvalue = lrt_p
    ))
  }

  for (a in analytes) {
    for (c in cutoffs) {
      res_ch <- get_cindex_comparison(df_ch, "Changhai", "BCR_time", "BCR_status", a, c)
      if (!is.null(res_ch)) cindex_results_list[[length(cindex_results_list) + 1]] <- res_ch
      
      res_tcga <- get_cindex_comparison(df_tcga, "TCGA", "BCR_time", "BCR_status", a, c)
      if (!is.null(res_tcga)) cindex_results_list[[length(cindex_results_list) + 1]] <- res_tcga
    }
  }
}

final_cindex_df <- bind_rows(cindex_results_list)
write.csv(final_cindex_df, file.path(result_dir, "cindex_comparison_results.csv"), row.names = FALSE)
cat("Completed C-index calculation. Saved to:",
    file.path(result_dir, "cindex_comparison_results.csv"), "\n")

# ------------------ 4. Build a single comparison figure (returns a ggplot object) ------------------
build_cindex_plot <- function(plot_data, expr_subset, cutoff_subset) {
  if (nrow(plot_data) == 0) return(NULL)
  
  # wide to long
  data_long <- rbind(
    plot_data %>% select(Cohort, Analyte, C_index = C_Clinical, SE = SE_Clinical) %>%
      mutate(Model = "Clinical only"),
    plot_data %>% select(Cohort, Analyte, C_index = C_Combined, SE = SE_Combined) %>%
      mutate(Model = "Clinical + Score")
  )
  data_long$Model   <- factor(data_long$Model, levels = c("Clinical only", "Clinical + Score"))
  data_long$Analyte <- factor(data_long$Analyte, levels = analytes)
  
  p <- ggplot(data_long, aes(x = Analyte, y = C_index, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.58) +
    geom_errorbar(aes(ymin = C_index - SE, ymax = C_index + SE),
                  position = position_dodge(0.7), width = 0.22,
                  linewidth = 0.35, colour = "black") +
    facet_wrap(~Cohort, scales = "free_y") +
    scale_fill_manual(values = color_map) +
    # leave headroom at the top for the delta and significance labels; do not hard-code the upper limit so the labels are not clipped
    scale_y_continuous(expand = expansion(mult = c(0, 0.20)), limits = c(0, NA)) +
    labs(
      x = NULL,
      y = "Concordance Index (C-index)",
      title = paste0(expr_subset, " (", cutoff_subset, " Cutoff)")
    ) +
    theme_nature_cindex()
  
  # labels: the delta C value plus the LRT-based significance stars, centred above each pair of bars
  delta_df <- plot_data %>%
    mutate(
      Analyte = factor(Analyte, levels = analytes),
      y_pos = pmax(C_Clinical + SE_Clinical, C_Combined + SE_Combined) + 0.03,
      star = vapply(LRT_Pvalue, sig_star, character(1)),
      delta_txt = ifelse(Delta_C >= 0, sprintf("+%.3f", Delta_C), sprintf("%.3f", Delta_C)),
      label_text = paste0(delta_txt, "\n", star)
    )
  
  p <- p + geom_text(
    data = delta_df,
    aes(x = Analyte, y = y_pos, label = label_text),
    inherit.aes = FALSE,
    size = 2.0,
    lineheight = 0.85,
    fontface = "bold",
    colour = "#1F4E9E",
    vjust = 0
  )
  
  return(p)
}

# ------------------ 5. Write one separate PDF per expr_type ------------------
# Nature layout: slightly taller to fit the rotated x-axis text and the significance annotations
w_in <- 140 / 25.4
h_in <- 90 / 25.4

for (expr in c("FPKM", "TPM")) {
  pdf_path <- file.path(figure_dir, paste0("CIndex_Comparison_", expr, "_AllCutoffs.pdf"))
  cairo_pdf(pdf_path, width = w_in, height = h_in, family = "sans", onefile = TRUE)
  
  for (cut in cutoffs) {
    plot_data <- final_cindex_df %>% filter(ExprType == expr, Cutoff == cut)
    p <- tryCatch(build_cindex_plot(plot_data, expr, cut),
                  error = function(e) {
                    cat("Plotting failed for", expr, "/", cut, ":", e$message, "\n")
                    NULL
                  })
    if (!is.null(p)) {
      print(p)
      cat("Rendered page:", expr, "-", cut, "\n")
    } else {
      cat("Skipped (no data):", expr, "-", cut, "\n")
    }
  }
  
  dev.off()
  cat("Saved PDF for", expr, "to:\n- ", pdf_path, "\n")
}

cat("\n========== All tasks completed successfully! ==========\n")
