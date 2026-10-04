# Prepare the input files for inferCNV
```R
library(Seurat)
load("/path/to/Spa-PR/PR/P616-loupe-much/1_Seurat-cluster-and-anno/P616_denoised_cluster.RData")
DefaultAssay(P616)<-"Spatial" # must be set, otherwise STEP 3 throws an error

# clusters (0,1,2,3,6,8,9,11,13) are tumour
# clusters (4,5,12) are normal


Idents(P616)<-"seurat_clusters"
P616.epi <- subset(P616, idents =c(0,1,2,3,6,8,9,11,13) )
P616.epi$celltype.case<-Idents(P616.epi)
Idents(P616.epi)<-"celltype.case"

P616.epi.split.list <- SplitObject(P616.epi, split.by = "orig.ident")

P616.fib <- subset(P616, idents = c(4,5,12))
P616.fib$celltype.case<-Idents(P616.fib)
Idents(P616.fib)<-"celltype.case"

P616.fib.split.list <- SplitObject(P616.fib, split.by = "orig.ident")

### P616 ###
epi_P616<-P616.epi.split.list[["P616"]]
DefaultAssay(epi_P616) <- "Spatial"
cellAnnota_P616_Epi <- subset(epi_P616@meta.data, select='celltype.case')
exprMatrix_P616_Epi <- as.matrix(GetAssayData(epi_P616[["Spatial"]], slot='counts'))


fib_P616<-P616.fib.split.list[["P616"]]
DefaultAssay(fib_P616) <- "Spatial"
cellAnnota_P616_fib <- subset(fib_P616@meta.data, select='celltype.case')
exprMatrix_P616_fib <- as.matrix(GetAssayData(fib_P616[["Spatial"]], slot='counts'))

cellAnnota_P616 <-rbind(cellAnnota_P616_Epi,cellAnnota_P616_fib)
exprMatrix_P616<-cbind(exprMatrix_P616_Epi, exprMatrix_P616_fib)
exprMatrix_P616 <- subset(exprMatrix_P616,select=rownames(cellAnnota_P616))


#for infercnv
write.table(exprMatrix_P616, paste('P616_P616','_exprMatrix.txt',sep=""), col.names=NA, sep='\t',quote = F)
write.table(cellAnnota_P616, paste('P616_P616','_cellAnnota.txt',sep=""), col.names=F, sep='\t',quote = F)
```
#Run inferCNV
```R
library(infercnv)
library(dplyr)


##Create inferCNV object
rm(list=ls())
options(stringsAsFactors = F)
expFile='P616_P616_exprMatrix.txt'
groupFiles='P616_P616_cellAnnota.txt'

options(scipen = 100)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= '/path/to/Spa-PR/PR/P60-loupe/3_inferCNV/hg38_gencode_v27.txt',
                                    ref_group_names=c("4","5","12")
)

save(infercnv_obj,file='infercnv_P616.rda')
#perform infercnv operations to reveal cnv signal
#load("infercnv_P616.rda")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='./output_dir_P616',
                             cluster_by_groups=F,
                             denoise=TRUE,
                             HMM=TRUE,
                             analysis_mode='subclusters',
                             tumor_subcluster_partition_method = "random_trees")



```
# Build the tree
```

```
# Annotate the tree
```R
library(Seurat)
load("/path/to/Spa-PR/PR/P60-loupe/1_Seurat-cluster-and-anno/P60_denoised_cluster.RData")
df<-read.table("/path/to/Spa-PR/PR/P60-loupe/3_inferCNV/output_dir_P60/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",header=T)
rownames(df)<-df$cell
df$cluster<-"Fib"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.1.1"),]$cluster<-"D"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.1.2"),]$cluster<-"E"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.2.1"),]$cluster<-"G"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.2.2"),]$cluster<-"H"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.1.1"),]$cluster<-"K"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.1.2"),]$cluster<-"L"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.2.1"),]$cluster<-"N"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.2.2"),]$cluster<-"O"


P60 <- AddMetaData(P60, metadata = df$cluster, col.name = "CNV_cluster")
table(P60$Sprod_snn_res.0.4,P60$CNV_cluster)
tb<-table(P60$Sprod_snn_res.0.4,P60$CNV_cluster)
tb<-as.data.frame(tb)
tb$M<-tb$all_observations.all_observations.1.2.2.2
p<-SpatialDimPlot(P60,group.by="CNV_cluster",label=T,label.size=3,repel=T)+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00"))
```
#Cut a Tree (Dendrogram/hclust/phylo) into Groups of Data
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
infercnv.dend <- read.dendrogram(file = "infercnv.observations_dendrogram.txt")
# Cut tree 
infercnv.labels <- cutree(infercnv.dend, k = 4, order_clusters_as_data = FALSE)
#infercnv.labels <- cutree(infercnv.dend, h = 60, order_clusters_as_data = FALSE)
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
load("/path/to/Spa-PR/PR/P60-loupe/1_Seurat-cluster-and-anno/P60_denoised_cluster.RData")
P60_no_Fib<-subset(P60,subset=Sprod_snn_res.0.4!=c(7))
library(ggplot2)
P60_no_Fib <- AddMetaData(P60_no_Fib, metadata = infercnv.labels, col.name = "CNV_cluster")
p<-SpatialDimPlot(P60_no_Fib,group.by="CNV_cluster",label=T,label.size=3,repel=T)
ggsave(p,file="CNV_cluster_plot.pdf")

df<-read.table("/path/to/Spa-PR/PR/P60-loupe/3_inferCNV/output_dir_P60/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",header=T)
rownames(df)<-df$cell
df$cluster<-"Fib"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.1.1"),]$cluster<-"D"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.1.2"),]$cluster<-"E"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.2.1"),]$cluster<-"G"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.2.2"),]$cluster<-"H"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.1.1"),]$cluster<-"K"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.1.2"),]$cluster<-"L"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.2.1"),]$cluster<-"N"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.2.2"),]$cluster<-"O"
df_no_Fib<-subset(df,subset=cluster!="Fib")
P60_no_Fib <- AddMetaData(P60_no_Fib, metadata = df_no_Fib$cluster, col.name = "CNV_cluster2")
table(P60_no_Fib$CNV_cluster,P60_no_Fib$CNV_cluster2)
```
