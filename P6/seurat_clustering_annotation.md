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
P60<- Load10X_Spatial(
     data.dir="/path/to/Spa-PR/1_loupe/P60-loupe/outs/",
     filename = "filtered_feature_bc_matrix.h5",
     assay = "Spatial",
     slice = "slice1",
     filter.matrix = TRUE,
     to.upper = FALSE,
)
Idents(P60)<-"orig.ident"
P60<- RenameIdents(P60, 'SeuratProject' = "P60")
P60$orig.ident<-Idents(P60)

# Normalisation
P60<- SCTransform(P60, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction and clustering
# umap
P60 <- RunPCA(P60, assay = "SCT", verbose = FALSE, dims = 1:30)
P60 <- FindNeighbors(P60, reduction = "pca", dims = 1:30)
P60 <- FindClusters(P60, verbose = FALSE, resolution = 0.6)
P60 <- RunUMAP(P60, reduction = "pca", dims = 1:30)
p1<-DimPlot(P60,reduction="umap",label=T)
p2<-SpatialDimPlot(P60,label=T,label.size=3)
ggsave(p1+p2,file="Dimension.pdf",width=14)

# Identify spatially variable features
P60<- FindSpatiallyVariableFeatures(P60, assay = "SCT", features = VariableFeatures(P60)[1:1000], 
    selection.method = "markvariogram")
# Visualise the top 20 features found by this method
top.features <- head(SpatiallyVariableFeatures(P60, selection.method = "markvariogram"), 20)
p3<-SpatialFeaturePlot(P60, features = top.features, ncol = 4, alpha = c(0.1, 1))
ggsave(p3,file="SpatialFeaturePlot.pdf")

#integrate scRNA-seq
P60_ref<-readRDS("/path/to/Spa-PR/scRNA_20210828.rds")

allen_reference <- P60_ref
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
P60<- SCTransform(P60, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
p1<-DimPlot(allen_reference, group.by = "celltype2", label = TRUE)
ggsave(p1,file="ref.pdf")
anchors <- FindTransferAnchors(reference = allen_reference, query = P60, normalization.method = "SCT",
                    reduction="cca" # "pcaproject": recommended when reference and query both come from scRNA-seq
                                   # "lsiproject": recommended when reference and query both come from scATAC-seq
                                   # "rpca": unsupervised; finds the optimal reconstruction subspace so that its principal components capture most of the variance of the sample
                                   # "cca": unsupervised; reduces both datasets and finds an optimally correlated subspace, suitable for cross-modality learning
)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$celltype2, prediction.assay = TRUE, 
                                  weight.reduction = P60[["pca"]],dims = 1:30)
P60[["predictions"]] <- predictions.assay
DefaultAssay(P60) <- "predictions"

group <- rep("cell", ncol(P60[["predictions"]][]))
names(group) <- colnames(P60[["predictions"]][])
for(i in 1:ncol(P60[["predictions"]][])){
group[i]<-rownames(P60[["predictions"]][])[which.max(P60[["predictions"]][,names(group)[i]])]
}
P60<- AddMetaData(P60, metadata = group, col.name = "celltype_max_score")
# Colour palette adapted from Liang Yuan (senior labmate), adjusted
my_cols <- c('B-cell'='#4B4BF7', 'Endothelial-cells'='#F68282', 'Epithelial-cells'='#B95FBB', 'Fiberblasts'='#1FA195', 'Monocyte'='#AC8F14', 'Neutrophils'='#A4DFF2', 'NK-cell'='#25aff5', 'T-cells'='#faf4cf', 'Pre-B-cell-CD34-'='#CCB1F1', 'subtype-0'='#aeadb3', 'subtype-1'='#ff9a36', 'subtype-2'='#31C53F', 'subtype-3'='#E6C122')

p4<-SpatialDimPlot(P60,group.by="celltype_max_score",label=T,label.size=3,cols=my_cols)
p5<-DimPlot(P60,reduction="umap",group.by="celltype_max_score",label=T,cols=my_cols)
ggsave(p4+p5,file="celltype_max_score.pdf",width=15,height=10)

# Run FindAllMarkers (NB: the Idents line below assigns to sRCC_P1 - a copy-paste leftover - so P60 keeps its own cluster Idents)
Idents(sRCC_P1)<-P60$celltype_max_score
 DefaultAssay(P60) <- "SCT"
all.markers <- FindAllMarkers(object = P60,only.pos = T,test.use = "MAST")
top20<-all.markers %>% group_by(cluster) %>% top_n(n=20,wt=avg_log2FC)
```
# Violin plots of the marker genes
```R

library(Seurat)
load("sRCC_P1_de_noise.RData")
spatial_count <- sRCC_P1@assays$Sprod@counts
spatial_count<-spatial_count[c("DNMT1","TET1","TET2","TET3","PDCD1LG2"),]
spatial_count<-as.matrix(spatial_count)
spatial_count<-t(spatial_count)

load("prop_data.RData")
spatial_count<-spatial_count[rownames(prop_data),]
spatial_count<-cbind(spatial_count,prop_data$s_r_related)
colnames(spatial_count)<-c("DNMT1","TET1","TET2","TET3","PDCD1LG2","s_r_related")

load("5_markergene_counts_matrix.RData")
spatial_count<-as.data.frame(spatial_count)
for (i in c(1:5)) {
  spatial_count[,i]<-as.numeric(spatial_count[,i])
}


# Load the packages
library(ggplot2)
library(tidyr)
library(ggpubr)
library(tidyverse)
library(hrbrthemes)
library(viridis)


cell_percent_boxplot <- function(i, tg_new, my_comparisons) {
  value<-i
  name<-'s_r_related' # change manually depending on what the x axis should mean
  tg_boxplot<-tg_new[,c(value,'s_r_related')]
  #boxplot
  tg_boxplot[,name]<-factor(tg_boxplot[,name],levels = c('normal_Epi','ccRCC','sRCC')) # change the levels to match the meaning of the x axis
  colnames(tg_boxplot)<-c("expression","s_r_related")
  p<-tg_boxplot %>%
    ggplot( aes(x=s_r_related, y=expression, fill=s_r_related)) +
    geom_violin() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(i) +
    xlab("")
  plot<-p+stat_compare_means(comparisons = my_comparisons,label = "p.signif")
  return(plot)
  
}
my_comparisons<-list(c('normal_Epi','ccRCC'),c('ccRCC','sRCC'),c('normal_Epi','sRCC'))
tg_colnames<-c("DNMT1","TET1","TET2","TET3","PDCD1LG2")

for (i in tg_colnames) {
  plot<-cell_percent_boxplot(i, tg_new=spatial_count, my_comparisons)
  assign(paste("p_", i, sep=""), plot)
}

library(showtext)
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto()

library(patchwork)
p<-p_DNMT1+p_TET1+p_TET2+p_TET3+p_PDCD1LG2+plot_layout(ncol=3,nrow=2) # adjust the panel layout
p

```
