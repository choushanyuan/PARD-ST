# Prepare the input files for inferCNV
```R
library(Seurat)
load("/path/to/Spa-PR/PR/P3/1_Seurat-cluster-and-anno/PR_P3_cluster_and_anno.RData")
DefaultAssay(PR_P3)<-"Spatial"

# Use the "Normal" cluster as the reference


Idents(PR_P3)<-"anno_cluster"
# Select the abnormal (malignant) epithelium
PR_P3.epi <- subset(PR_P3, idents = c("Normal","vascular_wall","Normal_with_inflammation"), invert=T) # drop the normal spots; the rest are the abnormal epithelium used for CNV
PR_P3.epi$celltype.case<-Idents(PR_P3.epi)
Idents(PR_P3.epi)<-"celltype.case"

PR_P3.epi.split.list <- SplitObject(PR_P3.epi, split.by = "orig.ident")
# Select the reference
PR_P3.fib <- subset(PR_P3, idents = c("Normal","Normal_with_inflammation"))
PR_P3.fib$celltype.case<-Idents(PR_P3.fib)
Idents(PR_P3.fib)<-"celltype.case"

PR_P3.fib.split.list <- SplitObject(PR_P3.fib, split.by = "orig.ident")

### PR_P3 ###
epi_PR_P3<-PR_P3.epi.split.list[["PR_P3"]]
DefaultAssay(epi_PR_P3) <- "Spatial"
cellAnnota_PR_P3_Epi <- subset(epi_PR_P3@meta.data, select='celltype.case')
exprMatrix_PR_P3_Epi <- as.matrix(GetAssayData(epi_PR_P3[["Spatial"]], slot='counts'))


fib_PR_P3<-PR_P3.fib.split.list[["PR_P3"]]
DefaultAssay(fib_PR_P3) <- "Spatial"
cellAnnota_PR_P3_fib <- subset(fib_PR_P3@meta.data, select='celltype.case')
exprMatrix_PR_P3_fib <- as.matrix(GetAssayData(fib_PR_P3[["Spatial"]], slot='counts'))

cellAnnota_PR_P3 <-rbind(cellAnnota_PR_P3_Epi,cellAnnota_PR_P3_fib)
exprMatrix_PR_P3<-cbind(exprMatrix_PR_P3_Epi, exprMatrix_PR_P3_fib)
exprMatrix_PR_P3 <- subset(exprMatrix_PR_P3,select=rownames(cellAnnota_PR_P3))


#for infercnv
write.table(exprMatrix_PR_P3, paste('PR_P3_PR_P3','_exprMatrix.txt',sep=""), col.names=NA, sep='\t',quote = F)
write.table(cellAnnota_PR_P3, paste('PR_P3_PR_P3','_cellAnnota.txt',sep=""), col.names=F, sep='\t',quote = F)
```
#Run inferCNV
```R
library(infercnv)
library(dplyr)


##Create inferCNV object
rm(list=ls())
options(stringsAsFactors = F)
expFile='PR_P3_PR_P3_exprMatrix.txt'
groupFiles='PR_P3_PR_P3_cellAnnota.txt'

options(scipen = 100)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= '/path/to/Spa-PR/PR/P60-loupe/3_inferCNV/hg38_gencode_v27.txt',
                                    ref_group_names=c("Normal","Normal_with_inflammation")
)

save(infercnv_obj,file='infercnv_PR_P3.rda')
#perform infercnv operations to reveal cnv signal
#load("infercnv_PR_P3.rda")
# Run per spot; used afterwards to build the dendrogram
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='./output_dir_PR_P3_for_tree',
                             cluster_by_groups=F,
                             denoise=TRUE,
                             HMM=TRUE,
                             analysis_mode='subclusters',
                             tumor_subcluster_partition_method = "random_trees")
# Run per cluster for a quick look
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='./output_dir_PR_P3',
                             cluster_by_groups=T,
                             denoise=TRUE,
                             HMM=F,
                             analysis_mode='sample',
                             tumor_subcluster_partition_method = "random_trees")
```

# Cut the tree to inspect the subclonal relationships

```R
rm(list=ls())
options(stringsAsFactors = F)
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Seurat)
library(infercnv)
library(miscTools)

#  Import inferCNV dendrogram
infercnv.dend <- read.dendrogram(file = "/path/to/Spa-PR/PR/P3/3_infercnv/output_dir_PR_P3_for_tree/infercnv.observations_dendrogram.txt")
# Cut tree 
infercnv.labels <- cutree(infercnv.dend, k = 5, order_clusters_as_data = FALSE)
#infercnv.labels <- cutree(infercnv.dend, h = 30, order_clusters_as_data = FALSE)
table(infercnv.labels)
# Color labels
the_bars <- as.data.frame(tableau_color_pal("Tableau 20")(20)[infercnv.labels])
colnames(the_bars) <- "inferCNV_tree"
the_bars$inferCNV_tree <- as.character(the_bars$inferCNV_tree)

pdf("inferCNV_dendrogram.pdf",height = 5,width = 10)
infercnv.dend %>% set("labels",rep("", nobs(infercnv.dend)) )  %>% plot(main="inferCNV dendrogram") %>%
  colored_bars(colors = as.data.frame(the_bars), dend = infercnv.dend, sort_by_labels_order = FALSE, add = T, y_scale=100 , y_shift = 0)
dev.off()
```
```R
load("/path/to/Spa-PR/PR/P3/1_Seurat-cluster-and-anno/PR_P3_cluster_and_anno.RData")
Idents(PR_P3)<-"anno_cluster"
# Select the abnormal (malignant) epithelium
PR_P3_no_Normal <- subset(PR_P3, idents = c("Normal","vascular_wall","Normal_with_inflammation"), invert=T)

library(ggplot2)
PR_P3_no_Normal <- AddMetaData(PR_P3_no_Normal, metadata = infercnv.labels, col.name = "CNV_cluster")
# Apply the colours of the dendrogram to the tissue section as well
rownames(the_bars)<-names(infercnv.labels)
PR_P3_no_Normal <- AddMetaData(PR_P3_no_Normal, metadata = the_bars, col.name = "CNV_color")
table(PR_P3_no_Normal$CNV_cluster,PR_P3_no_Normal$CNV_color)
p<-SpatialDimPlot(PR_P3_no_Normal,group.by="CNV_cluster",label=T,label.size=3,repel=T)+scale_fill_manual(values=c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F","#8CD17D","#B6992D","#F1CE63"))
ggsave(p,file="CNV_cluster_plot.pdf")
```
