# Clustering after de-noising
```R
# Load the R packages
library(Seurat)
library(ggplot2)
library(dplyr)

# Read the spatial transcriptomics data
PR_P4<- Load10X_Spatial(
  data.dir="/path/to/Spa-PR/PR/P4/0_PR_04-loupe-Spaceranger-results/PR_P4_Spaceranger_for_sprod/outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)
Idents(PR_P4)<-"orig.ident"
PR_P4 <- RenameIdents(PR_P4, 'SeuratProject' = "PR_P4")
PR_P4$orig.ident<-Idents(PR_P4)

# Read the de-noised matrix
tmp <- read.table("/path/to/Spa-PR/PR/P4/2_sprod_denoising/sprod_Denoised_matrix.txt"
                  , header=T, row.names = 1, sep="\t")
sprod_count <- as(as.matrix(t(tmp)), "dgCMatrix")
commongene <- intersect(rownames(sprod_count), rownames(PR_P4))
PR_P4 <- PR_P4[commongene,]
sprod_count <- sprod_count[commongene,]
identical(rownames(PR_P4), rownames(sprod_count))
PR_P4[['Sprod']] = CreateAssayObject(counts = sprod_count)
# Set the default assay
DefaultAssay(PR_P4)<-"Sprod"

# Filtering
mt.genes <- grep(pattern = "^MT-", x = rownames(PR_P4), value = TRUE)
PR_P4$percent.mito <- (Matrix::colSums(PR_P4@assays$Spatial@counts[mt.genes, ])/Matrix::colSums(PR_P4@assays$Spatial@counts))*100
#remove mt genes
objgene.2 <- Matrix::rowSums(PR_P4@assays$Spatial@counts != 0)
objgene.2 <- objgene.2[which(objgene.2 >= 5)]
genes_to_keep <- setdiff(names(objgene.2),mt.genes)
PR_P4 <- subset(PR_P4,features =genes_to_keep, subset = nFeature_Spatial > 200 & percent.mito < 20)

# Fix the data type of the coordinate columns
PR_P4@images$slice1@coordinates$row<-as.numeric(PR_P4@images$slice1@coordinates$row)
PR_P4@images$slice1@coordinates$col<-as.numeric(PR_P4@images$slice1@coordinates$col)
PR_P4@images$slice1@coordinates$imagerow<-as.numeric(PR_P4@images$slice1@coordinates$imagerow)
PR_P4@images$slice1@coordinates$imagecol<-as.numeric(PR_P4@images$slice1@coordinates$imagecol)
PR_P4@images$slice1@coordinates$tissue<-as.numeric(PR_P4@images$slice1@coordinates$tissue)


# Dimensionality reduction and clustering
# umap
PR_P4<- FindVariableFeatures(PR_P4,selection.method = "vst", nfeatures = 3000)
PR_P4<- ScaleData(PR_P4)

PR_P4 <- RunPCA(PR_P4, assay = "Sprod", verbose = FALSE, dims = 1:30)
PR_P4 <- FindNeighbors(PR_P4, reduction = "pca", dims = 1:30)
PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.4)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P4,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.4.pdf",width=14)

PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.6)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P4,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.6.pdf",width=14)

PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.8)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P4,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.8.pdf",width=14)

DefaultAssay(PR_P4)<-"Sprod"
p3<-SpatialFeaturePlot(PR_P4, features = c("EPCAM", "PTPRC","PECAM1","ACTA2","CD8A","CD8B","CD3D","CD3G","CD4","DNMT1","TET1","TET2","TET3","PDCD1LG2","VIM","CDH1","CDH2"),ncol=2)
ggsave(p3,file="marker_genes.pdf",width=18,height=45)
```
# Annotate with the scRNA-seq data (without de-noising)
```R
load("/path/to/Spa-PR/PCa_Single_Cell.rda")
allen_reference<-subset(PRADObj,subset=CellType!="Undefined") # drop the "Undefined" cells
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
PR_P4<- SCTransform(PR_P4, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
p1<-DimPlot(allen_reference, group.by = "CellType", label = TRUE)
ggsave(p1,file="ref.pdf")
anchors <- FindTransferAnchors(reference = allen_reference, query = PR_P4, normalization.method = "SCT",
                               reduction="cca" # "pcaproject": recommended when reference and query both come from scRNA-seq
                               # "lsiproject": recommended when reference and query both come from scATAC-seq
                               # "rpca": unsupervised; finds the optimal reconstruction subspace so that its principal components capture most of the variance of the sample
                               # "cca": unsupervised; reduces both datasets and finds an optimally correlated subspace, suitable for cross-modality learning
)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$CellType, prediction.assay = TRUE, 
                                  weight.reduction = PR_P4[["pca"]],dims = 1:30)
PR_P4[["predictions"]] <- predictions.assay
DefaultAssay(PR_P4) <- "predictions"

p2<-SpatialFeaturePlot(PR_P4, features = rownames(predictions.assay)[1:9], pt.size.factor = 1.6, ncol = 3, crop = TRUE)
ggsave(p2,file="anno.pdf",width=20,height=20)
```
# Clustering without de-noising
```R
# SCTransform normalisation
PR_P4 <- SCTransform(PR_P4, assay = "Spatial", verbose = FALSE)
DefaultAssay(PR_P4)<-"SCT"
# Dimensionality reduction and clustering
# umap

PR_P4 <- RunPCA(PR_P4, assay = "SCT", verbose = FALSE, dims = 1:30)
PR_P4 <- FindNeighbors(PR_P4, reduction = "pca", dims = 1:30)
PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.4)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P4,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.4.pdf",width=14)

PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.6)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P4,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.6.pdf",width=14)

PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.8)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P4,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.8.pdf",width=14)

PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.5)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P4,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.5.pdf",width=14)

DefaultAssay(PR_P4)<-"Sprod"
PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.6)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
plot1<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))

DefaultAssay(PR_P4)<-"SCT"
PR_P4 <- FindClusters(PR_P4, verbose = FALSE, resolution = 0.8)
PR_P4 <- RunUMAP(PR_P4, reduction = "pca", dims = 1:30)
plot2<-SpatialDimPlot(PR_P4,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))

ggsave(plot1+plot2,file="sprod_0.6_vs_unsprod_0.8.pdf",width=14)

save(PR_P4,file="PR_P4_cluster_and_sprod.RData")
```
# Annotation provided by the pathologist
```R
library(Seurat)
load("/path/to/Spa-PR/PR/P4/1_Seurat-cluster-and-anno/PR_P4_cluster_and_sprod.RData")
anno<-rep("cell", ncol(PR_P4))
names(anno)<-colnames(PR_P4)
anno[PR_P4$Sprod_snn_res.0.6 %in% c(5,13)] <- "GS4_forming_ductal_carcinoma"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(6,24,8,20,22)] <- "Normal_with_Tumor"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(16)] <- "Normal_with_immune"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(1)] <- "Normal_with_undefined_Tumor"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(19)] <- "GS5_1"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(21)] <- "GS4_1"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(7,9,15)] <- "Normal"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(23)] <- "GS5_2"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(10)] <- "GS4_and_GS5"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(12)] <- "GS5_3"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(18)] <- "GS4_2"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(0,17)] <- "GS4_with_GS5"
anno[PR_P4$Sprod_snn_res.0.6 %in% c(2,3,4,11,14)] <- "GS3"


PR_P4<- AddMetaData(PR_P4, metadata = anno, col.name = "anno_cluster")

library(ggplot2)
plot<-SpatialDimPlot(PR_P4,group.by="anno_cluster",label=T,label.size=3,repel=T,label.color = "black")+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#5C4033"))

ggsave(plot,file="anno_cluster_results.pdf")
save(PR_P4,file="PR_P4_cluster_and_anno.RData")
```
