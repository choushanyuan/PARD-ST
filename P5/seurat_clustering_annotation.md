# Clustering after de-noising
```R
# Load the R packages
library(Seurat)
library(ggplot2)
library(dplyr)

# Read the spatial transcriptomics data
PR_P5<- Load10X_Spatial(
  data.dir="/path/to/Spa-PR/PR/P5/0_PR_05-loupe-Spaceranger-results/PR_P5_Spaceranger_tif/outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)
Idents(PR_P5)<-"orig.ident"
PR_P5 <- RenameIdents(PR_P5, 'SeuratProject' = "PR_P5")
PR_P5$orig.ident<-Idents(PR_P5)

# Read the de-noised matrix
tmp <- read.table("/path/to/Spa-PR/PR/P5/2_sprod_denoising/output/sprod_Denoised_matrix.txt"
                  , header=T, row.names = 1, sep="\t")
sprod_count <- as(as.matrix(t(tmp)), "dgCMatrix")
commongene <- intersect(rownames(sprod_count), rownames(PR_P5))
PR_P5 <- PR_P5[commongene,]
sprod_count <- sprod_count[commongene,]
identical(rownames(PR_P5), rownames(sprod_count))
PR_P5[['Sprod']] = CreateAssayObject(counts = sprod_count)
# Set the default assay
DefaultAssay(PR_P5)<-"Sprod"

# Filtering
mt.genes <- grep(pattern = "^MT-", x = rownames(PR_P5), value = TRUE)
PR_P5$percent.mito <- (Matrix::colSums(PR_P5@assays$Spatial@counts[mt.genes, ])/Matrix::colSums(PR_P5@assays$Spatial@counts))*100
#remove mt genes
objgene.2 <- Matrix::rowSums(PR_P5@assays$Spatial@counts != 0)
objgene.2 <- objgene.2[which(objgene.2 >= 5)]
genes_to_keep <- setdiff(names(objgene.2),mt.genes)
PR_P5 <- subset(PR_P5,features =genes_to_keep, subset = nFeature_Spatial > 200 & percent.mito < 20)

# Fix the data type of the coordinate columns
PR_P5@images$slice1@coordinates$row<-as.numeric(PR_P5@images$slice1@coordinates$row)
PR_P5@images$slice1@coordinates$col<-as.numeric(PR_P5@images$slice1@coordinates$col)
PR_P5@images$slice1@coordinates$imagerow<-as.numeric(PR_P5@images$slice1@coordinates$imagerow)
PR_P5@images$slice1@coordinates$imagecol<-as.numeric(PR_P5@images$slice1@coordinates$imagecol)
PR_P5@images$slice1@coordinates$tissue<-as.numeric(PR_P5@images$slice1@coordinates$tissue)


# Dimensionality reduction and clustering
# umap
PR_P5<- FindVariableFeatures(PR_P5,selection.method = "vst", nfeatures = 3000)
PR_P5<- ScaleData(PR_P5)

PR_P5 <- RunPCA(PR_P5, assay = "Sprod", verbose = FALSE, dims = 1:30)
PR_P5 <- FindNeighbors(PR_P5, reduction = "pca", dims = 1:30)
PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.4)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P5,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.4.pdf",width=14)

PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.6)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P5,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.6.pdf",width=14)

PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.8)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P5,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_sprod_0.8.pdf",width=14)

DefaultAssay(PR_P5)<-"Sprod"
p3<-SpatialFeaturePlot(PR_P5, features = c("EPCAM", "PTPRC","PECAM1","ACTA2","CD8A","CD8B","CD3D","CD3G","CD4","DNMT1","TET1","TET2","TET3","PDCD1LG2","VIM","CDH1","CDH2"),ncol=2)
ggsave(p3,file="marker_genes.pdf",width=18,height=45)
```
# Annotate with the scRNA-seq data (without de-noising)
```R
load("/path/to/Spa-PR/PCa_Single_Cell.rda")
allen_reference<-subset(PRADObj,subset=CellType!="Undefined") # drop the "Undefined" cells
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
PR_P5<- SCTransform(PR_P5, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
p1<-DimPlot(allen_reference, group.by = "CellType", label = TRUE)
ggsave(p1,file="ref.pdf")
anchors <- FindTransferAnchors(reference = allen_reference, query = PR_P5, normalization.method = "SCT",
                               reduction="cca" # "pcaproject": recommended when reference and query both come from scRNA-seq
                               # "lsiproject": recommended when reference and query both come from scATAC-seq
                               # "rpca": unsupervised; finds the optimal reconstruction subspace so that its principal components capture most of the variance of the sample
                               # "cca": unsupervised; reduces both datasets and finds an optimally correlated subspace, suitable for cross-modality learning
)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$CellType, prediction.assay = TRUE, 
                                  weight.reduction = PR_P5[["pca"]],dims = 1:30)
PR_P5[["predictions"]] <- predictions.assay
DefaultAssay(PR_P5) <- "predictions"

p2<-SpatialFeaturePlot(PR_P5, features = rownames(predictions.assay)[1:9], pt.size.factor = 1.6, ncol = 3, crop = TRUE)
ggsave(p2,file="anno.pdf",width=20,height=20)
```
# Clustering without de-noising
```R
# SCTransform normalisation
PR_P5 <- SCTransform(PR_P5, assay = "Spatial", verbose = FALSE)
DefaultAssay(PR_P5)<-"SCT"
# Dimensionality reduction and clustering
# umap

PR_P5 <- RunPCA(PR_P5, assay = "SCT", verbose = FALSE, dims = 1:30)
PR_P5 <- FindNeighbors(PR_P5, reduction = "pca", dims = 1:30)
PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.4)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P5,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.4.pdf",width=14)

PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.6)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P5,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.6.pdf",width=14)

PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.8)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P5,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.8.pdf",width=14)

PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.5)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
p1<-DimPlot(PR_P5,reduction="umap",label=T)
p2<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))
ggsave(p1+p2,file="Dimension_0.5.pdf",width=14)

DefaultAssay(PR_P5)<-"Sprod"
PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.6)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
plot1<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))

DefaultAssay(PR_P5)<-"SCT"
PR_P5 <- FindClusters(PR_P5, verbose = FALSE, resolution = 0.8)
PR_P5 <- RunUMAP(PR_P5, reduction = "pca", dims = 1:30)
plot2<-SpatialDimPlot(PR_P5,label=T,label.size=3,pt.size.factor = 1.2,alpha=c(0.7,0.7))

ggsave(plot1+plot2,file="sprod_0.6_vs_unsprod_0.8.pdf",width=14)
```
# Annotation provided by the pathologist
```R
library(Seurat)
anno<-rep("GS5", ncol(PR_P5))
names(anno)<-colnames(PR_P5)
anno[PR_P5$Sprod_snn_res.0.6 %in% c(8,23,15,13,20,16,12)] <- "Normal"
anno[PR_P5$Sprod_snn_res.0.6 %in% c(4)] <- "GS4_with_some_GS5"
anno[PR_P5$Sprod_snn_res.0.6 %in% c(5)] <- "GS4_adn_GS5"

PR_P5<- AddMetaData(PR_P5, metadata = anno, col.name = "anno_cluster")
library(ggplot2)
plot<-SpatialDimPlot(PR_P5,group.by="anno_cluster",label=T,label.size=3,repel=T,label.color = "black")+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a"))

ggsave(plot,file="anno_cluster_results.pdf")
save(PR_P5,file="PR_P5_cluster_and_anno.RData")
```
