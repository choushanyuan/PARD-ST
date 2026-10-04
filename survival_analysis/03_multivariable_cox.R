# ============================================================
# Multivariable Cox regression + forest plot + model diagnostics
# Data: prostate cancer cohorts (Changhai and TCGA)
# Analytes: HSPD1, FASN, PKP1, Up_score, Down_score
# Grouping methods: Median, Q3_Upper, Q1_Lower, Q1_vs_Q3, Optimal
# New adjustment covariates: pathologic T stage (pT), pathologic N stage (pN), surgical margin (Margin)
# ============================================================
# ------------------ 1. Load the required packages ------------------
library(survival)
library(survminer)      # for surv_cutpoint (optional)
library(dplyr)
library(forestploter)   # publication-quality forest plots
library(grid)           # drawing backend for forestploter

# ------------------ 2. Set the working paths ------------------
base_dir <- "/path/to/survival_project"   # change this to your own path
result_dir <- file.path(base_dir, "result")
figure_dir <- file.path(base_dir, "figure")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------ 3. Linearity test (continuous covariates only) ------------------
test_linearity <- function(fit, var_name) {
  data <- fit$data
  if (!is.numeric(data[[var_name]])) return(NA)
  
  data[[paste0(var_name, "_sq")]] <- data[[var_name]]^2
  formula_orig <- formula(fit)
  new_formula <- update(formula_orig, as.formula(paste0("~ . + ", var_name, "_sq")))
  fit2 <- tryCatch(coxph(new_formula, data = data), error = function(e) NULL)
  if (is.null(fit2)) return(NA)
  
  lrt <- anova(fit, fit2)
  pval <- lrt[2, "Pr(>|Chi|)"]
  return(pval)
}

# ------------------ 4. Core function: fit the multivariable Cox model (age-stratified at 65) ------------------
run_multivariable_cox <- function(df, cohort_name, time_col, status_col, 
                                  score_col, cutoff_method, covar_cols,
                                  perform_diagnostic = TRUE) {
  # Arguments:
  #   df          : must already contain the Age_group column ("<65" or ">=65"),
  #                 and the pT/pN/Margin columns already cleaned into factors
  #   covar_cols  : excludes Age; contains the continuous variables (PSA/GleasonScore/TumorPurity)
  #                 and the factor variables (pT/pN/Margin)
  
  # ---- 4.1 Extract and clean the data ----
  df_clean <- df[, c(time_col, status_col, score_col, covar_cols, "Age_group")]
  colnames(df_clean) <- c("time", "status", "score", covar_cols, "Age_group")
  
  df_clean$time   <- as.numeric(df_clean$time)
  df_clean$status <- as.numeric(df_clean$status)
  df_clean$score  <- as.numeric(df_clean$score)
  
  # Handle each covariate according to its type:
  #   already a factor -> leave it as is (keeps the preset reference level)
  #   character        -> convert to factor
  #   anything else    -> convert to numeric
  for (cov in covar_cols) {
    if (is.factor(df_clean[[cov]])) {
      # keep the factor and the order of its levels unchanged
    } else if (is.character(df_clean[[cov]])) {
      df_clean[[cov]] <- factor(df_clean[[cov]])
    } else {
      df_clean[[cov]] <- as.numeric(df_clean[[cov]])
    }
  }
  
  # make sure Age_group is a factor with "<65" as the first level
  df_clean$Age_group <- factor(df_clean$Age_group, levels = c("<65", ">=65"))
  
  # Special case: truncate the Changhai follow-up at 70 months
  if (cohort_name == "Changhai") {
    df_clean <- df_clean[which(df_clean$time <= 70), ]
  }
  
  df_clean <- df_clean[complete.cases(df_clean), ]
  
  # drop the factor covariates that have only one level left after complete.cases, which would make coxph fail
  drop_covars <- c()
  for (cov in covar_cols) {
    if (is.factor(df_clean[[cov]])) {
      df_clean[[cov]] <- droplevels(df_clean[[cov]])
      if (nlevels(df_clean[[cov]]) < 2) drop_covars <- c(drop_covars, cov)
    }
  }
  active_covars <- setdiff(covar_cols, drop_covars)
  
  if (nrow(df_clean) < 10 || sum(df_clean$status == 1) < 3) return(NULL)
  
  # ---- 4.2 Group the samples according to cutoff_method ----
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
  
  df_clean$group <- factor(df_clean$group, levels = c("Low", "High"))
  df_clean$Expression <- df_clean$group
  
  # ---- 4.3 Fit the Cox model (age-stratified) ----
  # keep only the covariates that still have >= 2 levels in this subset
  formula_str <- paste("Surv(time, status) ~ Expression +", 
                       paste(active_covars, collapse = " + "), 
                       "+ strata(Age_group)")
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df_clean)
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NULL)
  fit$data <- df_clean   # keep the data for plotting and diagnostics
  
  # ---- 4.4 Model diagnostics ----
  diag <- list()
  if (perform_diagnostic) {
    zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if (!is.null(zph)) {
      diag$ph_global_p <- zph$table["GLOBAL", "p"]
      ph_vars <- rownames(zph$table)[rownames(zph$table) != "GLOBAL"]
      diag$ph_variable_p <- zph$table[ph_vars, "p", drop = FALSE]
    } else {
      diag$ph_global_p <- NA
      diag$ph_variable_p <- NA
    }
    
    linear_p <- sapply(active_covars, function(v) {
      if (is.numeric(df_clean[[v]])) test_linearity(fit, v) else NA
    })
    diag$linearity_p <- linear_p
  }
  
  # ---- 4.5 Extract the results ----
  s <- summary(fit)
  coef_matrix <- s$coefficients
  conf_matrix <- s$conf.int
  
  if (!"ExpressionHigh" %in% rownames(coef_matrix)) return(NULL)
  
  res_row <- data.frame(
    Cohort = cohort_name,
    Analyte = score_col,
    Cutoff = cutoff_method,
    N_patients = nrow(df_clean),
    N_events = sum(df_clean$status == 1),
    Group_HR = coef_matrix["ExpressionHigh", "exp(coef)"],
    Group_HR_Lower = conf_matrix["ExpressionHigh", "lower .95"],
    Group_HR_Upper = conf_matrix["ExpressionHigh", "upper .95"],
    Group_Pvalue = coef_matrix["ExpressionHigh", "Pr(>|z|)"]
  )
  
  # extract the HR and P value of every covariate (a column is kept for every requested covar_cols; missing values filled with NA)
  for (cov in covar_cols) {
    match_idx <- which(grepl(paste0("^", cov), rownames(coef_matrix)))
    if (length(match_idx) > 0) {
      res_row[[paste0(cov, "_HR")]]     <- coef_matrix[match_idx[1], "exp(coef)"]
      res_row[[paste0(cov, "_Pvalue")]] <- coef_matrix[match_idx[1], "Pr(>|z|)"]
    } else {
      res_row[[paste0(cov, "_HR")]]     <- NA
      res_row[[paste0(cov, "_Pvalue")]] <- NA
    }
  }
  
  if (perform_diagnostic) {
    res_row$PH_global_p <- diag$ph_global_p
    if (!is.null(diag$ph_variable_p) && is.data.frame(diag$ph_variable_p)) {
      ph_vals <- paste(rownames(diag$ph_variable_p), 
                       round(diag$ph_variable_p[,1], 4), 
                       sep = "=", collapse = "; ")
      res_row$PH_variable_p <- ph_vals
    } else {
      res_row$PH_variable_p <- NA
    }
    res_row$linearity_p <- paste(names(diag$linearity_p), 
                                 round(diag$linearity_p, 4), 
                                 sep = "=", collapse = "; ")
  }
  
  attr(res_row, "fit") <- fit
  attr(res_row, "diagnostics") <- diag
  
  return(res_row)
}

# ------------------ 5. Build the forest-plot table (reference row est = 1) ------------------
build_forest_table <- function(fit) {
  df   <- fit$data
  s    <- summary(fit)
  cf   <- s$coefficients
  ci   <- s$conf.int
  vars <- attr(fit$terms, "term.labels")
  vars <- vars[!grepl("strata", vars)]
  n_total <- nrow(df)
  
  chars <- c(); nums <- c(); hrtext <- c(); ptext <- c()
  est <- c(); low <- c(); up <- c()
  
  push <- function(ch, num, hr, pv, e, l, u) {
    chars  <<- c(chars, ch)
    nums   <<- c(nums, num)
    hrtext <<- c(hrtext, hr)
    ptext  <<- c(ptext, pv)
    est    <<- c(est, e)
    low    <<- c(low, l)
    up     <<- c(up, u)
  }
  
  fmt_p <- function(p) {
    if (is.na(p)) return("")
    stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
    p_text <- ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    paste0(p_text, stars)
  }
  
  fmt_hr <- function(e, l, u) sprintf("%.2f [%.2f, %.2f]", e, l, u)
  
  for (v in vars) {
    x <- df[[v]]
    if (is.factor(x) || is.character(x)) {
      x <- factor(x)
      lv <- levels(x)
      ref <- lv[1]
      push(v, "", "", "", NA, NA, NA)
      for (l in lv) {
        cnt <- sum(x == l, na.rm = TRUE)
        pct <- round(100 * cnt / n_total, 1)
        numstr <- paste0(cnt, " (", pct, ")")
        if (l == ref) {
          push(paste0("  ", l), numstr, "Reference", "", 1, 1, 1)
        } else {
          rn <- paste0(v, l)
          e <- ci[rn, "exp(coef)"]; lo <- ci[rn, "lower .95"]; hi <- ci[rn, "upper .95"]
          p <- cf[rn, "Pr(>|z|)"]
          push(paste0("  ", l), numstr, fmt_hr(e, lo, hi), fmt_p(p), e, lo, hi)
        }
      }
    } else {
      m <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
      numstr <- sprintf("%.2f (%.2f)", m, sdv)
      rn <- v
      e <- ci[rn, "exp(coef)"]; lo <- ci[rn, "lower .95"]; hi <- ci[rn, "upper .95"]
      p <- cf[rn, "Pr(>|z|)"]
      push(v, numstr, fmt_hr(e, lo, hi), fmt_p(p), e, lo, hi)
    }
  }
  
  df_out <- data.frame(
    Characteristics = chars,
    Statistics      = nums,
    `HR (95%CI)`    = hrtext,
    `P.value`       = ptext,
    est = est, low = low, up = up,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  
  # Relabel the variables
  df_out$Characteristics <- ifelse(df_out$Characteristics == "PSA", 
                                   "log2(PSA+1)", df_out$Characteristics)
  df_out$Characteristics <- ifelse(df_out$Characteristics == "TumorPurity", 
                                   "Tumorpurity(per 10%)", df_out$Characteristics)
  df_out$Characteristics <- ifelse(df_out$Characteristics == "pT", 
                                   "Pathologic T", df_out$Characteristics)
  df_out$Characteristics <- ifelse(df_out$Characteristics == "pN", 
                                   "Pathologic N", df_out$Characteristics)
  df_out$Characteristics <- ifelse(df_out$Characteristics == "Margin", 
                                   "Surgical margin", df_out$Characteristics)
  
  return(df_out)
}

# ------------------ 6. Draw the forest plot (forestploter style) ------------------
plot_forest_style <- function(fit, title,
                              xlim = c(0.1, 20),
                              ticks = c(0.1, 0.5, 1, 2, 5, 20)) {
  # note: adding pN (N1) can produce large HRs, so the default xlim is widened on purpose
  dat <- build_forest_table(fit)
  dat$` ` <- paste(rep(" ", 22), collapse = " ")
  show_cols <- c("Characteristics", "Statistics", " ", "HR (95%CI)", "P.value")
  disp <- dat[, show_cols]
  colnames(disp) <- c("Characteristics", "n (%) / Mean (SD)", " ", "HR (95%CI)", "P.value")
  tm <- forest_theme(
    base_size = 10,
    ci_pch = 18,
    ci_col = "#1F4E9E",
    ci_lwd = 1.6,
    ci_Theight = 0.2,
    refline_gp = gpar(lwd = 1, lty = "solid", col = "black"),
    footnote_gp = gpar(cex = 0.8),
    arrow_type = "closed",
    core = list(bg_params = list(fill = "white"))
  )
  
  p <- forest(
    disp,
    est = dat$est,
    lower = dat$low,
    upper = dat$up,
    ci_column = 3,
    ref_line = 1,
    x_trans = "log",
    xlim = xlim,
    ticks_at = ticks,
    xlab = "Hazard Ratio",
    title = title,
    theme = tm
  )
  return(p)
}

# ------------------ 6b. Helper functions to clean the stage / margin variables ------------------
# pathologic T: contains T2 -> "T2" (reference), contains T3/T4 -> "T3-T4", otherwise NA
clean_pT <- function(x) {
  x <- toupper(trimws(as.character(x)))
  g <- ifelse(grepl("T2", x), "T2",
              ifelse(grepl("T3|T4", x), "T3-T4", NA))
  factor(g, levels = c("T2", "T3-T4"))
}
# pathologic N: contains N0 -> "N0" (reference), contains N1 -> "N1", NX/blank -> NA
clean_pN <- function(x) {
  x <- toupper(trimws(as.character(x)))
  g <- ifelse(grepl("N0", x), "N0",
              ifelse(grepl("N1", x), "N1", NA))
  factor(g, levels = c("N0", "N1"))
}
# margin (Changhai): NO -> Negative (reference), YES -> Positive, otherwise NA
clean_margin_changhai <- function(x) {
  x <- toupper(trimws(as.character(x)))
  g <- ifelse(x == "NO", "Negative",
              ifelse(x == "YES", "Positive", NA))
  factor(g, levels = c("Negative", "Positive"))
}
# margin (TCGA): R0 -> Negative (reference), R1/R2 -> Positive, RX/blank -> NA
clean_margin_tcga <- function(x) {
  x <- toupper(trimws(as.character(x)))
  g <- ifelse(x == "R0", "Negative",
              ifelse(x %in% c("R1", "R2"), "Positive", NA))
  factor(g, levels = c("Negative", "Positive"))
}

# ------------------ 7. Main loop over all analytes and cutoffs ------------------
analytes <- c("HSPD1", "FASN", "PKP1", "Up_score", "Down_score")
cutoffs <- c("Median", "Q3_Upper", "Q1_Lower", "Q1_vs_Q3", "Optimal")

# the covariate list now includes pT / pN / Margin
covars <- c("PSA", "GleasonScore", "TumorPurity", "pT", "pN", "Margin")

for (expr_type in c("FPKM", "TPM")) {
  cat("\n========== Processing ", expr_type, " ==========\n")
  
  load(file.path(result_dir, expr_type, "cohort_clinical_with_scores.RData"))
  
  all_results <- list()
  fit_list <- list()
  
  # ------------------ 7.1 Changhai cohort pre-processing ------------------
  df_ch <- cohort_data$Changhai
  df_ch <- df_ch[!is.na(df_ch$Age), ]
  df_ch$BCR_status <- ifelse(df_ch$BCR == "YES", 1, ifelse(df_ch$BCR == "NO", 0, NA))
  df_ch$BCR_time <- as.numeric(gsub("m", "", df_ch$`Month to surgery/BCR`))
  df_ch$GleasonScore <- as.numeric(as.character(df_ch$GS))
  df_ch$PSA <- as.numeric(as.character(df_ch$PSA))
  df_ch$PSA <- log2(df_ch$PSA + 1)
  df_ch$TumorPurity <- df_ch$TumorPurity * 10   # HR per 0.1 (i.e. per 10%) increase
  df_ch$Age_group <- ifelse(df_ch$Age >= 65, ">=65", "<65")
  
  # clean the pathologic T / N / margin variables (column names as given in the original report)
  df_ch$pT     <- clean_pT(df_ch$`Pathology T Stage`)
  df_ch$pN     <- clean_pN(df_ch$`Pathology N stage`)
  df_ch$Margin <- clean_margin_changhai(df_ch$`Surgical Margin`)
  
  for (a in analytes) {
    for (c in cutoffs) {
      res <- run_multivariable_cox(df_ch, "Changhai", "BCR_time", "BCR_status", 
                                   a, c, covars, perform_diagnostic = TRUE)
      if (!is.null(res)) {
        all_results[[length(all_results) + 1]] <- res
        fit <- attr(res, "fit")
        if (!is.null(fit)) {
          fit_list[[paste("Changhai", a, c, sep = "_")]] <- fit
        }
      }
    }
  }
  
  # ------------------ 7.2 TCGA cohort pre-processing ------------------
  df_tcga <- cohort_data$TCGA
  df_tcga <- df_tcga[!is.na(df_tcga$age_at_initial_pathologic_diagnosis), ]
  df_tcga$BCR_status <- df_tcga$PFI
  df_tcga$BCR_time <- as.numeric(df_tcga$PFI.time) / 30.4   # convert days to months
  df_tcga$Age <- df_tcga$age_at_initial_pathologic_diagnosis
  df_tcga$GleasonScore <- as.numeric(as.character(df_tcga$gleason_score))
  df_tcga$PSA <- log2(as.numeric(as.character(df_tcga$psa_value)) + 1)
  df_tcga$TumorPurity <- df_tcga$TumorPurity * 10
  df_tcga$Age_group <- ifelse(df_tcga$Age >= 65, ">=65", "<65")
  
  # clean the pathologic T / N / margin variables (columns pathologic_T / pathologic_N / residual_tumor)
  df_tcga$pT     <- clean_pT(df_tcga$pathologic_T)
  df_tcga$pN     <- clean_pN(df_tcga$pathologic_N)
  df_tcga$Margin <- clean_margin_tcga(df_tcga$residual_tumor)
  
  for (a in analytes) {
    for (c in cutoffs) {
      res <- run_multivariable_cox(df_tcga, "TCGA", "BCR_time", "BCR_status", 
                                   a, c, covars, perform_diagnostic = TRUE)
      if (!is.null(res)) {
        all_results[[length(all_results) + 1]] <- res
        fit <- attr(res, "fit")
        if (!is.null(fit)) {
          fit_list[[paste("TCGA", a, c, sep = "_")]] <- fit
        }
      }
    }
  }
  
  # ------------------ 7.3 Save the results ------------------
  final_df <- bind_rows(all_results)
  dir.create(file.path(result_dir, expr_type), showWarnings = FALSE, recursive = TRUE)
  write.csv(final_df, file.path(result_dir, expr_type, "multivariable_cox_results_stratified_age65.csv"), 
            row.names = FALSE)
  
  if (nrow(final_df) > 0) {
    diag_summary <- final_df %>%
      select(Cohort, Analyte, Cutoff, PH_global_p, PH_variable_p, linearity_p)
    write.csv(diag_summary, 
              file.path(result_dir, expr_type, "model_diagnostics_stratified_age65.csv"), 
              row.names = FALSE)
  }
  
  # ------------------ 7.4 Draw the forest plots ------------------
  if (length(fit_list) > 0) {
    dir.create(file.path(figure_dir, expr_type), showWarnings = FALSE, recursive = TRUE)
    # cairo_pdf() embeds the font subset so the exported PDF stays editable
    # (text remains text); base pdf() does not embed fonts and the file cannot
    # be edited in Illustrator/Acrobat.
    cairo_pdf(file.path(figure_dir, expr_type, "Multivariable_Cox_Forest_StratifiedAge65.pdf"),
        width = 10, height = 6, family = "sans", onefile = TRUE)
    
    for (name in names(fit_list)) {
      fit <- fit_list[[name]]
      tryCatch({
        p <- plot_forest_style(fit, title = name)
        if (!is.null(p)) {
          grid.newpage()
          grid.draw(p)
        }
      }, error = function(e) {
        cat("Forest plot failed:", name, " - ", e$message, "\n")
      })
    }
    dev.off()
    cat("Forest plots saved to:", file.path(figure_dir, expr_type, "Multivariable_Cox_Forest_StratifiedAge65.pdf"), "\n")
  }
  
  cat("Finished the ", expr_type, " analysis!\n")
}

cat("\n========== All done! ==========\n")
