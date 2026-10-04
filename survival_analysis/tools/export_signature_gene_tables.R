# Export the two signature gene lists as supplementary tables for the response letter
#   (Reviewer 2, comment 4: "the complete gene lists used for the Up and Down
#    signatures should be provided in a supplementary table").
#
# Direction mapping - verified, do NOT swap:
#   * The RData variable names follow the RAW pseudotime direction.
#   * The manuscript / revised survival pipeline defines the signatures by the
#     malignant-dedifferentiation direction, where they are inverted:
#         Up signature   = down_genes (n = 95)  = "伪时间去分化上调基因集" (pseudotime dedifferentiation up-regulated gene set)
#                                                 (pseudotime state-end up-regulated genes)
#                                                 -> used for Up_score in 01_calculate_scores.R
#         Down signature = up_genes   (n = 246) = down-regulated with dedifferentiation
#                                                 -> used for Down_score in 01_calculate_scores.R
#   Cross-checks: the column "伪时间去分化上调基因集" (pseudotime dedifferentiation up-regulated gene set) of the workbook "4个代表恶性特征基因集汇总.xlsx" (summary of 4 representative malignant-feature gene sets) starts
#   AMACR, CRISP3, AIDA, SPON2, ... == down_genes; and
#   down_genes == P1_down ∩ P3_down ∩ P4_down, up_genes == P1_up ∩ P3_up ∩ P4_up.

base_dir <- "/path/to/survival_project"
out_dir  <- "/path/to/response_letter_tables"

load(file.path(base_dir, "P1_P3_P4_overlap_pesudotime_genes_new.RData"))
stopifnot(length(down_genes) == 95, length(up_genes) == 246)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_table <- function(genes) {
  data.frame(No. = seq_along(genes), `Gene symbol` = genes,
             check.names = FALSE, stringsAsFactors = FALSE)
}

up_sig   <- make_table(down_genes)   # Up signature   (n = 95)
down_sig <- make_table(up_genes)     # Down signature (n = 246)

writexl::write_xlsx(list(Up_signature = up_sig),
                    file.path(out_dir, "Table_S_Up_signature_gene_list_n95.xlsx"))
writexl::write_xlsx(list(Down_signature = down_sig),
                    file.path(out_dir, "Table_S_Down_signature_gene_list_n246.xlsx"))

write.csv(up_sig,   file.path(out_dir, "Table_S_Up_signature_gene_list_n95.csv"), row.names = FALSE)
write.csv(down_sig, file.path(out_dir, "Table_S_Down_signature_gene_list_n246.csv"), row.names = FALSE)

cat("Up signature   (pseudotime state-end up-regulated, used for Up_score)   n =", nrow(up_sig),
    "\n  first:", paste(head(down_genes, 8), collapse = ", "),
    "\n  last :", paste(tail(down_genes, 5), collapse = ", "), "\n")
cat("Down signature (down with dedifferentiation, used for Down_score)      n =", nrow(down_sig),
    "\n  first:", paste(head(up_genes, 8), collapse = ", "),
    "\n  last :", paste(tail(up_genes, 5), collapse = ", "), "\n\n")
cat("duplicated entries:", sum(duplicated(down_genes)), "/", sum(duplicated(up_genes)), "\n")
cat("written to:", out_dir, "\n")
