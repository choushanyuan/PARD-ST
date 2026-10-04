# Clustering after de-noising
```R
# Load the R packages
library(Seurat)
library(ggplot2)
library(dplyr)

# Read the spatial transcriptomics data
PR_P2<- Load10X_Spatial(
  data.dir="/path/to/Spa-PR/PR/P2/0_PR_02-loupe-Spaceranger-results/PR_P2_Spaceranger_for_sprod/outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)
Idents(PR_P2)<-"orig.ident"
PR_P2 <- RenameIdents(PR_P2, 'SeuratProject' = "PR_P2")
PR_P2$orig.ident<-Idents(PR_P2)

# Read the de-noised matrix
tmp <- read.table("/path/to/Spa-PR/PR/P2/2_sprod_denoising/sprod_Denoised_matrix.txt"
                  , header=T, row.names = 1, sep="\t")
sprod_count <- as(as.matrix(t(tmp)), "dgCMatrix")
commongene <- intersect(rownames(sprod_count), rownames(PR_P2))
PR_P2 <- PR_P2[commongene,]
sprod_count <- sprod_count[commongene,]
identical(rownames(PR_P2), rownames(sprod_count))
PR_P2[['Sprod']] = CreateAssayObject(counts = sprod_count)
# Set the default assay
DefaultAssay(PR_P2)<-"Sprod"

# Filtering
mt.genes <- grep(pattern = "^MT-", x = rownames(PR_P2), value = TRUE)
PR_P2$percent.mito <- (Matrix::colSums(PR_P2@assays$Spatial@counts[mt.genes, ])/Matrix::colSums(PR_P2@assays$Spatial@counts))*100
#remove mt genes
objgene.2 <- Matrix::rowSums(PR_P2@assays$Spatial@counts != 0)
objgene.2 <- objgene.2[which(objgene.2 >= 5)]
genes_to_keep <- setdiff(names(objgene.2),mt.genes)
PR_P2 <- subset(PR_P2,features =genes_to_keep, subset = nFeature_Spatial > 200 & percent.mito < 20)

# Fix the data type of the coordinate columns
PR_P2@images$slice1@coordinates$row<-as.numeric(PR_P2@images$slice1@coordinates$row)
PR_P2@images$slice1@coordinates$col<-as.numeric(PR_P2@images$slice1@coordinates$col)
PR_P2@images$slice1@coordinates$imagerow<-as.numeric(PR_P2@images$slice1@coordinates$imagerow)
PR_P2@images$slice1@coordinates$imagecol<-as.numeric(PR_P2@images$slice1@coordinates$imagecol)
PR_P2@images$slice1@coordinates$tissue<-as.numeric(PR_P2@images$slice1@coordinates$tissue)


# Dimensionality reduction and clustering
# umap
PR_P2<- FindVariableFeatures(PR_P2,selection.method = "vst", nfeatures = 3000)
PR_P2<- ScaleData(PR_P2)

PR_P2 <- RunPCA(PR_P2, assay = "Sprod", verbose = FALSE, dims = 1:30)
PR_P2 <- FindNeighbors(PR_P2, reduction = "pca", dims = 1:30)
PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.4)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P2,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.4.pdf",width=14)

PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.6)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P2,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.6.pdf",width=14)

PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.8)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P2,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.8.pdf",width=14)

DefaultAssay(PR_P2)<-"Sprod"
p3<-SpatialFeaturePlot(PR_P2, features = c("EPCAM", "PTPRC","PECAM1","ACTA2","CD8A","CD8B","CD3D","CD3G","CD4","DNMT1","TET1","TET2","TET3","PDCD1LG2","VIM","CDH1","CDH2"),ncol=2)
ggsave(p3,file="marker_genes.pdf",width=18,height=45)
```
# Annotate with the scRNA-seq data (without de-noising)
```R
load("/path/to/Spa-PR/PCa_Single_Cell.rda")
allen_reference<-subset(PRADObj,subset=CellType!="Undefined") # drop the "Undefined" cells
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
PR_P2<- SCTransform(PR_P2, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
p1<-DimPlot(allen_reference, group.by = "CellType", label = TRUE)
ggsave(p1,file="ref.pdf")
anchors <- FindTransferAnchors(reference = allen_reference, query = PR_P2, normalization.method = "SCT",
                               reduction="cca" # "pcaproject": recommended when reference and query both come from scRNA-seq
                               # "lsiproject": recommended when reference and query both come from scATAC-seq
                               # "rpca": unsupervised; finds the optimal reconstruction subspace so that its principal components capture most of the variance of the sample
                               # "cca": unsupervised; reduces both datasets and finds an optimally correlated subspace, suitable for cross-modality learning
)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$CellType, prediction.assay = TRUE, 
                                  weight.reduction = PR_P2[["pca"]],dims = 1:30)
PR_P2[["predictions"]] <- predictions.assay
DefaultAssay(PR_P2) <- "predictions"

p2<-SpatialFeaturePlot(PR_P2, features = rownames(predictions.assay)[1:9], pt.size.factor = 1.6, ncol = 3, crop = TRUE)
ggsave(p2,file="anno.pdf",width=20,height=20)
```
# Clustering without de-noising
```R
# SCTransform normalisation
PR_P2 <- SCTransform(PR_P2, assay = "Spatial", verbose = FALSE)
DefaultAssay(PR_P2)<-"SCT"
# Dimensionality reduction and clustering
# umap

PR_P2 <- RunPCA(PR_P2, assay = "SCT", verbose = FALSE, dims = 1:30)
PR_P2 <- FindNeighbors(PR_P2, reduction = "pca", dims = 1:30)
PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.4)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P2,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.4.pdf",width=14)

PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.6)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P2,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.6.pdf",width=14)

PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.8)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P2,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.8.pdf",width=14)

PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.5)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P2,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.5.pdf",width=14)

DefaultAssay(PR_P2)<-"Sprod"
PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.6)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
plot1<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))

DefaultAssay(PR_P2)<-"SCT"
PR_P2 <- FindClusters(PR_P2, verbose = FALSE, resolution = 0.8)
PR_P2 <- RunUMAP(PR_P2, reduction = "pca", dims = 1:30)
plot2<-SpatialDimPlot(PR_P2,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))

ggsave(plot1+plot2,file="sprod_0.6_vs_unsprod_0.8.pdf",width=14)
```
# Annotation provided by the pathologist
```R
library(Seurat)
load("/path/to/Spa-PR/PR/P2/1_Seurat-cluster-and-anno/PR_P2_cluster_and_sprod.RData")
anno<-rep("cell", ncol(PR_P2))
names(anno)<-colnames(PR_P2)
anno[PR_P2$Sprod_snn_res.0.6 %in% c(9)] <- "GS4_forming_ductal_carcinoma"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(20)] <- "GS4_1"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(10)] <- "GS4_2"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(11)] <- "Tumor_1"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(5,12)] <- "Normal_with_Tumor"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(4)] <- "Tumor_2"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(19)] <- "Tumor_3"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(1)] <- "Tumor_4"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(16)] <- "Tumor_5"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(18)] <- "Tumor_6"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(3,17,8,0,2,22)] <- "undefined"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(21)] <- "inflammation"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(13,14,15,6)] <- "Normal_with_inflammation"
anno[PR_P2$Sprod_snn_res.0.6 %in% c(7)] <- "Normal"



# Extract the barcodes of between_PIN_and_ductal_carcinoma
between_PIN_and_ductal_carcinoma<- Load10X_Spatial(
  data.dir="/path/to/Spa-PR/PR/P2/0_PR_02-loupe-Spaceranger-results/PR_P2_Spaceranger_between_PIN_and_ductal_carcinoma/outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)

anno[colnames(PR_P2) %in% colnames(between_PIN_and_ductal_carcinoma)] <- "between_PIN_and_ductal_carcinoma"

# Extract the barcodes of ductal_carcinoma
ductal_carcinoma<- Load10X_Spatial(
  data.dir="/path/to/Spa-PR/PR/P2/0_PR_02-loupe-Spaceranger-results/PR_P2_Spaceranger_ductal_carcinoma/outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)

anno[colnames(PR_P2) %in% colnames(ductal_carcinoma)] <- "ductal_carcinoma"

# Extract the barcodes of PIN
PIN<- Load10X_Spatial(
  data.dir="/path/to/Spa-PR/PR/P2/0_PR_02-loupe-Spaceranger-results/PR_P2_Spaceranger_PIN/outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)

anno[colnames(PR_P2) %in% colnames(PIN)] <- "PIN"

# Extract the barcodes of vascular_wall
vascular_wall<- Load10X_Spatial(
  data.dir="/path/to/Spa-PR/PR/P2/0_PR_02-loupe-Spaceranger-results/PR_P2_Spaceranger_vascular_wall/outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)

anno[colnames(PR_P2) %in% colnames(vascular_wall)] <- "vascular_wall"

PR_P2<- AddMetaData(PR_P2, metadata = anno, col.name = "anno_cluster")
library(ggplot2)
plot<-SpatialDimPlot(PR_P2,group.by="anno_cluster",label=T,label.size=3,repel=T,label.color = "black")+scale_fill_manual(values=c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA'))
library(ggplot2)
ggsave(plot,file="anno_cluster_results.pdf")
save(PR_P2,file="PR_P2_cluster_and_anno.RData")

```
