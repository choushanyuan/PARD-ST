```R
# Install the R packages
library(Seurat)
library(SpatialCPie)
library(Spaniel)
library(SingleR)
library(infercnv)
library(clustree)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(devtools)
library(Matrix)
library(cowplot)
library(SeuratData)
library(patchwork)
library(dplyr)
library(hdf5r)

# Read the spatial transcriptomics data
P616<- Load10X_Spatial(
  data.dir="/path/to/Spa-PR/1_space/P616_l_much/outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)
Idents(P616)<-"orig.ident"
P616<- RenameIdents(P616, 'SeuratProject' = "P616")
P616$orig.ident<-Idents(P616)

# Normalisation
P616<- SCTransform(P616, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction and clustering
# umap
P616 <- RunPCA(P616, assay = "SCT", verbose = FALSE, dims = 1:30)
P616 <- FindNeighbors(P616, reduction = "pca", dims = 1:30)
P616 <- FindClusters(P616, verbose = FALSE, resolution = 0.8)
P616 <- RunUMAP(P616, reduction = "pca", dims = 1:30)
p1<-DimPlot(P616,reduction="umap",label=T)
p2<-SpatialDimPlot(P616,label=T,label.size=3)
ggsave(p1+p2,file="Dimension.pdf",width=14)

# Identify spatially variable features
P616<- FindSpatiallyVariableFeatures(P616, assay = "SCT", features = VariableFeatures(P616)[1:1000], 
                                    selection.method = "markvariogram")
# Visualise the top 20 features found by this method
top.features <- head(SpatiallyVariableFeatures(P616, selection.method = "markvariogram"), 20)
p3<-SpatialFeaturePlot(P616, features = top.features, ncol = 4, alpha = c(0.1, 1))
ggsave(p3,file="SpatialFeaturePlot.pdf")

#integrate scRNA-seq
P616_ref<-readRDS("/path/to/Spa-PR/scRNA_20210828.rds")

allen_reference <- P616_ref
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
P616<- SCTransform(P616, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
p1<-DimPlot(allen_reference, group.by = "celltype2", label = TRUE)
ggsave(p1,file="ref.pdf")
anchors <- FindTransferAnchors(reference = allen_reference, query = P616, normalization.method = "SCT",
                               reduction="cca" # "pcaproject": recommended when reference and query both come from scRNA-seq
                               # "lsiproject": recommended when reference and query both come from scATAC-seq
                               # "rpca": unsupervised; finds the optimal reconstruction subspace so that its principal components capture most of the variance of the sample
                               # "cca": unsupervised; reduces both datasets and finds an optimally correlated subspace, suitable for cross-modality learning
)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$celltype2, prediction.assay = TRUE, 
                                  weight.reduction = P616[["pca"]],dims = 1:30)
P616[["predictions"]] <- predictions.assay
DefaultAssay(P616) <- "predictions"

group <- rep("cell", ncol(P616[["predictions"]][]))
names(group) <- colnames(P616[["predictions"]][])
for(i in 1:ncol(P616[["predictions"]][])){
  group[i]<-rownames(P616[["predictions"]][])[which.max(P616[["predictions"]][,names(group)[i]])]
}
P616<- AddMetaData(P616, metadata = group, col.name = "celltype_max_score")
# Colour palette adapted from Liang Yuan (senior labmate), adjusted
my_cols <- c('B-cell'='#4B4BF7', 'Endothelial-cells'='#F68282', 'Epithelial-cells'='#B95FBB', 'Fiberblasts'='#1FA195', 'Monocyte'='#AC8F14', 'Neutrophils'='#A4DFF2', 'NK-cell'='#25aff5', 'T-cells'='#faf4cf', 'Pre-B-cell-CD34-'='#CCB1F1', 'subtype-0'='#aeadb3', 'subtype-1'='#ff9a36', 'subtype-2'='#31C53F', 'subtype-3'='#E6C122')
p4<-SpatialDimPlot(P616,group.by="celltype_max_score",label=T,label.size=3)
p5<-DimPlot(P616,reduction="umap",group.by="celltype_max_score",label=T)
ggsave(p4+p5,file="celltype_max_score.pdf",width=15,height=10)

# Run FindAllMarkers (NB: the Idents line below assigns to sRCC_P1 - a copy-paste leftover - so P616 keeps its own cluster Idents)
Idents(sRCC_P1)<-P616$celltype_max_score
DefaultAssay(P616) <- "SCT"
all.markers <- FindAllMarkers(object = P616,only.pos = T,test.use = "MAST")
top20<-all.markers %>% group_by(cluster) %>% top_n(n=20,wt=avg_log2FC)
```
