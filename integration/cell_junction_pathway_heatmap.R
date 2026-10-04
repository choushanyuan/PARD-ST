library(Seurat)
library(dplyr)
library(ggplot2)

load("/path/to/PR/integrated/prad_merged_integrated.rda")
Idents(prad)<-prad$integrated_snn_res.0.4

pure_Tumor<-subset(prad,subset=Tumour_purity > 0.75)
pure_Tumor<-subset(pure_Tumor,subset=integrated_snn_res.0.4 %in% c(0,1,2,9,12))
anno<-rep("cell", ncol(pure_Tumor))
names(anno)<-colnames(pure_Tumor)

anno[pure_Tumor$integrated_snn_res.0.4 %in% c(1)] <- "subtype_low"
anno[pure_Tumor$integrated_snn_res.0.4 %in% c(0,2)] <- "subtype_middle"
anno[pure_Tumor$integrated_snn_res.0.4 %in% c(9,12)] <- "subtype_high"
pure_Tumor<- AddMetaData(pure_Tumor, metadata = anno, col.name = "tumor_subtype")
pure_Tumor$tumor_subtype<-factor(pure_Tumor$tumor_subtype,levels=c("subtype_low","subtype_middle","subtype_high"))
cat("after purity+cluster filter:", ncol(pure_Tumor), "spots\n")
print(table(pure_Tumor$tumor_subtype))

# GSVA of the cell-junction-related pathways
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

h <- msigdbr(species = "Homo sapiens", category = "C5")
H <- select(h, gs_name, gene_symbol) %>%
  as.data.frame %>%
  split(., .$gs_name) %>%
  lapply(., function(x)(x$gene_symbol))
H<-H[grep("GOBP|GO", names(H))]
H_GO<-H[grep("EXTRACELLULAR|DESMOSOME", names(H))]

h <- msigdbr(species = "Homo sapiens", category = "C2")
H <- select(h, gs_name, gene_symbol) %>%
  as.data.frame %>%
  split(., .$gs_name) %>%
  lapply(., function(x)(x$gene_symbol))
H<-H[grep("KEGG", names(H))]
H_KEGG<-H[grep("ADHERENS|JUNCTION|ADHESION|ECM", names(H))]
H<-c(H_KEGG,H_GO)

gs <- lapply(H, unique)
count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 3))
gs <- lapply(gs, function(x) intersect(keep, x))
gs <- gs[lapply(gs, length) > 0]
sel <- c("KEGG_ADHERENS_JUNCTION","KEGG_CELL_ADHESION_MOLECULES_CAMS","KEGG_FOCAL_ADHESION","KEGG_GAP_JUNCTION","KEGG_TIGHT_JUNCTION","GOBP_DESMOSOME_ORGANIZATION","GOCC_DESMOSOME","GOCC_HEMIDESMOSOME")
cat("selected pathways present in filtered gs:", sum(sel %in% names(gs)), "/", length(sel), "\n")
if (any(!sel %in% names(gs))) stop("MISSING PATHWAYS: ", paste(sel[!sel %in% names(gs)], collapse=","))
for (s in sel) cat(s, "genes:", length(gs[[s]]), "\n")

ST_obj_tumor<-pure_Tumor
DefaultAssay(ST_obj_tumor)<-"Spatial"
expr=as.matrix(ST_obj_tumor@assays$Spatial@counts)
expr <-log2(expr+1)
uni_matrix <- expr[which(rowSums(expr) > 0),]
rm(expr); gc()

gsva_es <- gsva(uni_matrix,gs,
                kcdf="Gaussian",
                verbose=T,
                parallel.sz = parallel::detectCores()
)

gsva_es_selected_pathway<-gsva_es[c("KEGG_ADHERENS_JUNCTION",
"KEGG_CELL_ADHESION_MOLECULES_CAMS",
"KEGG_FOCAL_ADHESION",
"KEGG_GAP_JUNCTION",
"KEGG_TIGHT_JUNCTION",
"GOBP_DESMOSOME_ORGANIZATION",
"GOCC_DESMOSOME",
"GOCC_HEMIDESMOSOME"),]

anno_col<-data.frame(ST_obj_tumor$tumor_subtype)
colnames(anno_col)<-"anno"
anno_col$barcode<-rownames(anno_col)
anno_col<-anno_col[order(anno_col$anno),]
gsva_es_selected_pathway<-gsva_es_selected_pathway[,anno_col$barcode]
rownames(anno_col)<-colnames(gsva_es_selected_pathway)

gsva_es_mean<-data.frame(pathway=rownames(gsva_es_selected_pathway))
for (i in names(table(anno_col$anno))) {
  gsva_es_mean<-cbind(gsva_es_mean,data.frame(apply(gsva_es_selected_pathway[,anno_col[which(anno_col$anno==i),]$barcode],1,mean)))
}
gsva_es_mean<-gsva_es_mean[,-1]
colnames(gsva_es_mean)<-names(table(anno_col$anno))
cat("mean matrix col order:", colnames(gsva_es_mean), "\n")

library(pheatmap)
p<-pheatmap(gsva_es_mean,fontsize=12,fontsize_row=10,cluster_rows = T,cluster_cols = F,
            color=colorRampPalette(c("#0571b0","#f7f7f7","#ca0020"))(100),
            show_colnames = T,border_color = "grey",scale = "row",show_rownames =T,
            angle_col = 90,
            gaps_col = c(2),
            cutree_rows=2,
            cellwidth = 15, cellheight = 10,
            filename="cell_junction_GSVA.pdf", width=8, height=6)
cat("CELL JUNCTION HEATMAP SAVED OK\n")
