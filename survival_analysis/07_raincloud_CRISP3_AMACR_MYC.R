# Nature-style raincloud: CRISP3 / AMACR / MYC expression by 3 Gleason strata
# Group 1: GS 6 + GS 3+4
# Group 2: GS 4+3 + GS 8
# Group 3: GS 9 + GS 10
# Pairwise two-sided Wilcoxon rank-sum tests are shown above each pair.
# Expected location: <project_root>/survival_analysis/07_raincloud_CRISP3_AMACR_MYC.R

library(ggplot2)
library(ggdist)
library(dplyr)
library(patchwork)
library(ggpubr)

if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}

find_base_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg_idx <- grep("^--file=", args)
  if (length(file_arg_idx) > 0) {
    script_path <- sub("^--file=", "", args[file_arg_idx[1]])
    return(dirname(dirname(normalizePath(script_path, winslash = "/"))))
  }
  ok <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
  if (ok) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(dirname(normalizePath(ctx$path, winslash = "/"))))
    }
  }
  stop("Could not locate the project root. Run this file from <root>/survival_analysis/.")
}

env_project <- Sys.getenv("RAINCLOUD_PROJECT_ROOT", unset = "")
base_dir <- if (nzchar(env_project)) env_project else find_base_dir()
result_dir <- file.path(base_dir, "result")
figure_dir <- file.path(base_dir, "figure", "CRISP3_AMACR_MYC")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

genes_of_interest <- c("CRISP3", "AMACR", "MYC")
width_mm <- 183
height_mm <- 92
export_dpi <- 600

save_pub_r <- function(fig, stem, width_mm = 183, height_mm = 92,
                       dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h,
                       family = "Arial", onefile = TRUE)
  print(fig)
  dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(fig)
  dev.off()

  raster_path <- paste0(stem, "_600dpi.tiff")
  suppressWarnings({
    ragg::agg_tiff(raster_path, width = w, height = h, units = "in",
                   res = dpi, compression = "lzw")
    print(fig)
    dev.off()
  })
  if (!file.exists(raster_path)) {
    # Fallback: base R tiff device to an ASCII temp file, then copy.
    tmp <- tempfile(fileext = ".tiff")
    grDevices::tiff(tmp, width = round(w * dpi), height = round(h * dpi),
                    res = dpi, compression = "lzw")
    print(fig)
    dev.off()
    ok <- file.copy(tmp, raster_path, overwrite = TRUE)
    file.remove(tmp)
    if (!ok) warning("Could not write raster export to: ", raster_path)
  }
}

# ---------- Nature-style theme (R track) ----------
palette_contract <- c(
  neutral_dark = "#272727",
  neutral_mid = "#767676",
  neutral_light = "#D8D8D8",
  signal_blue = "#3182BD",
  signal_teal = "#33B5A5",
  accent_red = "#D24B40",
  accent_orange = "#E28E2C"
)

theme_nature_contract <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      axis.text.x = element_text(size = base_size - 1.2, colour = "black"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.7),
      strip.text = element_text(size = base_size - 0.3, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold",
                                hjust = 0.5, colour = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(t = 4, r = 4, b = 2, l = 4)
    )
}

# ---------- Gleason classification ----------
gleason_six <- function(primary, secondary, gs_score) {
  p <- suppressWarnings(as.numeric(as.character(primary)))
  s <- suppressWarnings(as.numeric(as.character(secondary)))
  gs <- suppressWarnings(as.numeric(as.character(gs_score)))
  total <- p + s
  total[is.na(total)] <- gs[is.na(total)]

  out <- rep(NA_character_, length(total))
  out[!is.na(total) & total <= 6] <- "GS 6"
  out[!is.na(total) & total == 7 & p == 3 & s == 4] <- "GS 3+4"
  out[!is.na(total) & total == 7 & p == 4 & s == 3] <- "GS 4+3"
  out[!is.na(total) & total == 8] <- "GS 8"
  out[!is.na(total) & total == 9] <- "GS 9"
  out[!is.na(total) & total == 10] <- "GS 10"
  factor(out, levels = c("GS 6", "GS 3+4", "GS 4+3",
                         "GS 8", "GS 9", "GS 10"))
}

map_three_groups <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x %in% c("GS 6", "GS 3+4")] <- "g1"
  out[!is.na(x) & x %in% c("GS 4+3", "GS 8")] <- "g2"
  out[!is.na(x) & x %in% c("GS 9", "GS 10")] <- "g3"
  factor(out, levels = c("g1", "g2", "g3"))
}

# ---------- pairwise stats ----------
sig_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

pairwise_wilcox <- function(dat) {
  levs <- c("g1", "g2", "g3")
  pairs <- list(c("g1", "g2"), c("g2", "g3"), c("g1", "g3"))
  out <- lapply(pairs, function(pr) {
    a <- dat[["Expr"]][dat[["Group"]] == pr[1]]
    b <- dat[["Expr"]][dat[["Group"]] == pr[2]]
    pv <- tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA_real_)
    data.frame(g1 = pr[1], g2 = pr[2], p = pv, label = sig_star(pv))
  })
  bind_rows(out)
}

# ---------- panel ----------
raincloud_panel <- function(dat, gene, y_label, ymax_room = 0.18) {
  dat <- dat[!is.na(dat[["Group"]]), , drop = FALSE]
  if (nrow(dat) < 2) return(NULL)

  cnt <- table(factor(dat[["Group"]], levels = c("g1", "g2", "g3")))
  base_labs <- c("GS6 & 3+4", "GS4+3 & 8", "GS9 & 10")
  x_labs <- sprintf("%s\n(n=%d)", base_labs, as.integer(cnt))

  stats <- pairwise_wilcox(dat)
  ymin <- min(dat[["Expr"]], na.rm = TRUE)
  ymax <- max(dat[["Expr"]], na.rm = TRUE)
  yrng <- diff(range(c(ymin, ymax)))
  if (!is.finite(yrng) || yrng == 0) yrng <- 1

  y_lo <- ymin - yrng * 0.04
  y_hi <- ymax + yrng * ymax_room
  bracket_pos <- c(ymax + yrng * 0.03,
                   ymax + yrng * 0.09,
                   ymax + yrng * 0.16)
  bracket_df <- data.frame(
    group1 = stats[["g1"]],
    group2 = stats[["g2"]],
    y.position = bracket_pos,
    label = stats[["label"]]
  )

  p <- ggplot(dat, aes(x = Group, y = Expr)) +
    stat_halfeye(adjust = 0.55, width = 0.55,
                 justification = -0.18, .width = 0,
                 fill = palette_contract[["neutral_light"]],
                 color = NA, alpha = 0.8) +
    geom_boxplot(width = 0.10, outlier.shape = NA,
                 fill = "white", colour = palette_contract[["neutral_dark"]],
                 linewidth = 0.3) +
    geom_jitter(width = 0.09, height = 0, size = 0.22,
                alpha = 0.30, colour = palette_contract[["signal_blue"]]) +
    geom_bracket(data = bracket_df,
                 aes(xmin = group1, xmax = group2,
                     y.position = y.position, label = label),
                 colour = palette_contract[["neutral_dark"]],
                 size = 0.25, label.size = 2.2,
                 tip.length = 0.015) +
    labs(title = gene, y = y_label, x = NULL) +
    scale_x_discrete(drop = FALSE, labels = x_labs) +
    coord_cartesian(ylim = c(y_lo, y_hi), clip = "off") +
    theme_nature_contract()

  attr(p, "pairwise_stats") <- stats
  p
}

# ---------- main loop ----------
for (expr_type in c("FPKM", "TPM")) {
  load(file.path(result_dir, expr_type, "cohort_clinical_with_scores.RData"))

  for (cohort in c("Changhai", "TCGA")) {
    clin <- cohort_data[[cohort]]

    if (cohort == "Changhai") {
      expr_mat <- expr_ch
      g1_col <- "G1"
      g2_col <- "G2"
      gs_col <- "GS"
      common <- intersect(colnames(expr_mat), rownames(clin))
      clin <- clin[common, , drop = FALSE]
      expr_sub <- expr_mat[genes_of_interest, common, drop = FALSE]
    } else {
      expr_mat <- expr_tcga
      g1_col <- "primary_pattern"
      g2_col <- "secondary_pattern"
      gs_col <- "gleason_score"
      pref <- substr(colnames(expr_mat), 1, 15)
      common <- intersect(rownames(clin), pref)
      idx <- match(common, pref)
      clin <- clin[common, , drop = FALSE]
      expr_sub <- expr_mat[genes_of_interest, idx, drop = FALSE]
    }

    grade6 <- gleason_six(clin[[g1_col]], clin[[g2_col]], clin[[gs_col]])
    group3 <- map_three_groups(grade6)

    long_list <- vector("list", length(genes_of_interest))
    for (i in seq_along(genes_of_interest)) {
      g <- genes_of_interest[i]
      long_list[[i]] <- data.frame(
        Sample = rownames(clin),
        Cohort = cohort,
        ExprType = expr_type,
        Gene = g,
        Expr = as.numeric(expr_sub[g, ]),
        Group = group3,
        stringsAsFactors = FALSE
      )
    }
    long <- bind_rows(long_list)

    y_label <- if (expr_type == "FPKM") {
      "Expression (log2 FPKM+1)"
    } else {
      "Expression (log2 TPM+1)"
    }

    cat("\n[", cohort, "/", expr_type, "] samples:",
        length(unique(long[["Sample"]])), "\n", sep = "")
    cnt <- table(factor(long[["Group"]][long[["Gene"]] == "CRISP3"],
                        levels = c("g1", "g2", "g3")), useNA = "ifany")
    cat("3-group counts: "); print(cnt)

    panels <- list()
    for (g in genes_of_interest) {
      dat <- long[long[["Gene"]] == g, , drop = FALSE]
      p <- raincloud_panel(dat, g, y_label)
      if (!is.null(p)) panels[[g]] <- p
    }

    fig <- wrap_plots(panels, ncol = 3) +
      plot_annotation(tag_levels = "a") &
      theme(plot.tag = element_text(size = 8, face = "bold"))

    stem <- file.path(figure_dir, paste0(cohort, "_", expr_type,
                                         "_CRISP3_AMACR_MYC_raincloud"))
    save_pub_r(fig, stem, width_mm, height_mm, export_dpi)

    # console QA: pairwise tests
    for (g in genes_of_interest) {
      st <- attr(panels[[g]], "pairwise_stats")
      cat(g, ": ")
      cat(paste(sprintf("%s vs %s p=%.3g (%s)", st[["g1"]], st[["g2"]],
                        st[["p"]], st[["label"]]), collapse = "; "), "\n")
    }

    cat("Saved:", paste0(stem, ".pdf / .svg / _600dpi.tiff"), "\n")
  }
}

cat("All Nature-style 3-group raincloud PDFs completed.\n")
