suppressMessages({library(Seurat); library(msigdbr); library(dplyr); library(data.table); library(GSVA); library(limma); library(stringr); library(ggplot2); library(pheatmap); library(RColorBrewer)})

# ---------- 1. load & subset (no purity filter per user decision; exclude PIN cluster 5) ----------
load("/path/to/PR/integrated/prad_merged_integrated.rda")
# purity filter intentionally SKIPPED (user decision 2026-09-09): P2(PR_P2)/P7(P616) have all-NA Tumour_purity
pure_Tumor <- subset(prad, subset = integrated_snn_res.0.4 %in% c(0,1,2,9,12))  # no cluster5/PIN
anno <- rep("cell", ncol(pure_Tumor)); names(anno) <- colnames(pure_Tumor)
anno[pure_Tumor$integrated_snn_res.0.4 %in% c(1)]   <- "subtype_low"
anno[pure_Tumor$integrated_snn_res.0.4 %in% c(0,2)] <- "subtype_middle"
anno[pure_Tumor$integrated_snn_res.0.4 %in% c(9,12)]<- "subtype_high"
pure_Tumor <- AddMetaData(pure_Tumor, anno, col.name = "tumor_subtype")
pure_Tumor$tumor_subtype <- factor(pure_Tumor$tumor_subtype, levels = c("subtype_low","subtype_middle","subtype_high"))

patient_map <- c("PR_P1"="P1","PR_P2"="P2","PR_P3"="P3","PR_P4"="P4","PR_P5"="P5","P60"="P6","P616"="P7")
pt <- patient_map[as.character(pure_Tumor$orig.ident)]
pure_Tumor <- AddMetaData(pure_Tumor, pt, col.name = "patient")
pure_Tumor$patient <- factor(pure_Tumor$patient, levels = paste0("P",1:7))

cat("=== spots per subtype x patient ===\n")
print(table(pure_Tumor$tumor_subtype, pure_Tumor$patient))
flush.console()

# ---------- 2. gene sets (same construction as original scripts) ----------
# hallmark (category H)
h  <- msigdbr(species="Homo sapiens", category="H")
H  <- select(h, gs_name, gene_symbol) %>% as.data.frame %>% split(., .$gs_name) %>% lapply(., function(x) x$gene_symbol)
gs_h <- lapply(H, unique)
keep_h <- names(which(table(unlist(gs_h)) < 2)); gs_h <- lapply(gs_h, function(x) intersect(keep_h, x)); gs_h <- gs_h[lapply(gs_h, length) > 0]

# cell-junction (C5 GO + C2 KEGG)
h5 <- msigdbr(species="Homo sapiens", category="C5")
H5 <- select(h5, gs_name, gene_symbol) %>% as.data.frame %>% split(., .$gs_name) %>% lapply(., function(x) x$gene_symbol)
H5 <- H5[grep("GOBP|GO", names(H5))]; H5_GO <- H5[grep("EXTRACELLULAR|DESMOSOME", names(H5))]
h2 <- msigdbr(species="Homo sapiens", category="C2")
H2 <- select(h2, gs_name, gene_symbol) %>% as.data.frame %>% split(., .$gs_name) %>% lapply(., function(x) x$gene_symbol)
H2 <- H2[grep("KEGG", names(H2))]; H2_KEGG <- H2[grep("ADHERENS|JUNCTION|ADHESION|ECM", names(H2))]
Hj <- c(H2_KEGG, H5_GO)
gs_j <- lapply(Hj, unique)
keep_j <- names(which(table(unlist(gs_j)) < 3)); gs_j <- lapply(gs_j, function(x) intersect(keep_j, x)); gs_j <- gs_j[lapply(gs_j, length) > 0]

gs_all <- c(gs_h, gs_j)
cat("total gene sets for GSVA:", length(gs_all), "\n"); flush.console()

# ---------- 3. expression matrix & GSVA (once) ----------
DefaultAssay(pure_Tumor) <- "Spatial"
expr <- as.matrix(pure_Tumor@assays$Spatial@counts)
expr <- log2(expr + 1)
uni_matrix <- expr[which(rowSums(expr) > 0),]; rm(expr); gc()

gsva_es <- gsva(uni_matrix, gs_all, kcdf="Gaussian", verbose=T, parallel.sz = parallel::detectCores())

# ---------- 4. select pathways ----------
hall_rows <- c("APICAL_JUNCTION","EPITHELIAL_MESENCHYMAL_TRANSITION","PI3K_AKT_MTOR_SIGNALING","TGF_BETA_SIGNALING","WNT_BETA_CATENIN_SIGNALING","NOTCH_SIGNALING","MYC_TARGETS_V1","MYC_TARGETS_V2","HEDGEHOG_SIGNALING","G2M_CHECKPOINT","MITOTIC_SPINDLE","E2F_TARGETS","P53_PATHWAY","APOPTOSIS","DNA_REPAIR","UV_RESPONSE_UP","UV_RESPONSE_DN","HYPOXIA","ANGIOGENESIS","INFLAMMATORY_RESPONSE","KRAS_SIGNALING_UP","KRAS_SIGNALING_DN","CHOLESTEROL_HOMEOSTASIS","REACTIVE_OXYGEN_SPECIES_PATHWAY","GLYCOLYSIS","FATTY_ACID_METABOLISM")
junc_rows <- c("KEGG_ADHERENS_JUNCTION","KEGG_CELL_ADHESION_MOLECULES_CAMS","KEGG_FOCAL_ADHESION","KEGG_GAP_JUNCTION","KEGG_TIGHT_JUNCTION","GOBP_DESMOSOME_ORGANIZATION","GOCC_DESMOSOME","GOCC_HEMIDESMOSOME")
rownames(gsva_es) <- str_replace(rownames(gsva_es), "HALLMARK_", "")
m_h <- gsva_es[hall_rows, , drop=FALSE]
m_j <- gsva_es[junc_rows, , drop=FALSE]
if (any(!hall_rows %in% rownames(gsva_es))) stop("missing hallmark rows")
if (any(!junc_rows %in% rownames(gsva_es))) stop("missing junction rows")

# ---------- 5. column means: subtype-major, patient P1..P7 inside ----------
subtype_levels <- c("subtype_low","subtype_middle","subtype_high")
bc <- colnames(pure_Tumor)
allc <- character()
for (s in subtype_levels) for (p in paste0("P",1:7)) {
  sel <- bc[pure_Tumor$tumor_subtype==s & as.character(pure_Tumor$patient)==p]
  if (length(sel) > 0) allc <- c(allc, paste0(s,"_",p))
}
cat("columns:", paste(allc, collapse=" | "), "\n"); flush.console()
mean_h <- sapply(allc, function(cn) {ss <- strsplit(cn,"_")[[1]]; s <- paste(ss[1:(length(ss)-1)], collapse="_"); p <- ss[length(ss)]
  sel <- bc[pure_Tumor$tumor_subtype==s & as.character(pure_Tumor$patient)==p]; rowMeans(m_h[, sel, drop=FALSE])})
mean_j <- sapply(allc, function(cn) {ss <- strsplit(cn,"_")[[1]]; s <- paste(ss[1:(length(ss)-1)], collapse="_"); p <- ss[length(ss)]
  sel <- bc[pure_Tumor$tumor_subtype==s & as.character(pure_Tumor$patient)==p]; rowMeans(m_j[, sel, drop=FALSE])})
rownames(mean_h) <- hall_rows; rownames(mean_j) <- junc_rows
subtype_of <- sapply(strsplit(allc,"_"), function(x) paste(x[-length(x)], collapse="_"))
patient_of <- sapply(strsplit(allc,"_"), function(x) x[length(x)])
gaps <- cumsum(rle(subtype_of)$lengths); gaps <- gaps[-length(gaps)]
cat("gap cols:", paste(gaps, collapse=","), " ncol:", ncol(mean_h), "\n"); flush.console()

anno_col_df <- data.frame(Subtype = factor(subtype_of, levels=subtype_levels),
                          Patient = factor(patient_of, levels=paste0("P",1:7)),
                          row.names = allc, check.names=FALSE)
sub_cols <- c("subtype_low"="#4393c3","subtype_middle"="#bdbdbd","subtype_high"="#d6604d")
pat_cols <- setNames(brewer.pal(7,"Set2"), paste0("P",1:7))
ann_colors <- list(Subtype=sub_cols, Patient=pat_cols)

# row annotation for hallmark (fixed 26-row category order)
anno_row <- data.frame(Pathway = c(rep("Migration and invasion pathways",2), rep("PI3K_AKT_mTOR pathway",1), rep("Stemness pathways",6), rep("Cell cycle pathways",5), rep("DNA repair pathways",3), rep("Other important pathways",9)))
rownames(anno_row) <- hall_rows
pathwaycolor <- c("Migration and invasion pathways"='#8dd3c7',"PI3K_AKT_mTOR pathway"='#ffffb3',"Stemness pathways"='#bebada',"Cell cycle pathways"='#fb8072',"DNA repair pathways"='#80b1d3',"Other important pathways"='#fdb462')
ann_colors_row <- list(Pathway=pathwaycolor)

heatcol <- colorRampPalette(c("#0571b0","#f7f7f7","#ca0020"))(100)
ann_col_all <- c(list(Subtype=sub_cols, Patient=pat_cols), list(Pathway=pathwaycolor))

# ---------- 6. hallmark heatmap (no PIN, by patient) ----------
p1 <- pheatmap(mean_h, cluster_rows=T, cluster_cols=F,
      color=heatcol, border_color="grey", scale="row",
      show_colnames=T, show_rownames=T, fontsize=12, fontsize_row=10,
      angle_col=90, gaps_col=gaps,
      annotation_col=anno_col_df, annotation_row=anno_row,
      annotation_colors=ann_col_all,
      labels_col=patient_of,
      cellwidth=20, cellheight=15,
      filename="tumor_GSVA_by_patient.pdf", width=13, height=10)
cat("FIG1 hallmark by-patient SAVED\n"); flush.console()

# ---------- 7. cell junction heatmap (by patient) ----------
p2 <- pheatmap(mean_j, cluster_rows=T, cluster_cols=F, cutree_rows=2,
      color=heatcol, border_color="grey", scale="row",
      show_colnames=T, show_rownames=T, fontsize=12, fontsize_row=10,
      angle_col=90, gaps_col=gaps,
      annotation_col=anno_col_df, annotation_colors=ann_colors,
      labels_col=patient_of,
      cellwidth=25, cellheight=22,
      filename="cell_junction_GSVA_by_patient.pdf", width=13, height=7)
cat("FIG2 junction by-patient SAVED\n")
cat("ALL DONE\n")
