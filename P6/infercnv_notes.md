# Prepare the input files for inferCNV
```R
library(Seurat)
load("/path/to/Spa-PR/PR/P60-loupe/1_Seurat-cluster-and-anno/P60_denoised_cluster.RData")
DefaultAssay(P60)<-"Spatial"
#SpatialDimPlot(P60 , group.by = "seurat_clusters",label = TRUE, label.size = 1)
#cluster7 is fibroblast


Idents(P60)<-"seurat_clusters"
P60.epi <- subset(P60, idents = c("7"), invert=T) # drop cluster 7
P60.epi$celltype.case<-Idents(P60.epi)
Idents(P60.epi)<-"celltype.case"

P60.epi.split.list <- SplitObject(P60.epi, split.by = "orig.ident")

P60.fib <- subset(P60, idents = c("7"))
P60.fib$celltype.case<-Idents(P60.fib)
Idents(P60.fib)<-"celltype.case"

P60.fib.split.list <- SplitObject(P60.fib, split.by = "orig.ident")

### P60 ###
epi_P60<-P60.epi.split.list[["P60"]]
DefaultAssay(epi_P60) <- "Spatial"
cellAnnota_P60_Epi <- subset(epi_P60@meta.data, select='celltype.case')
exprMatrix_P60_Epi <- as.matrix(GetAssayData(epi_P60[["Spatial"]], slot='counts'))


fib_P60<-P60.fib.split.list[["P60"]]
DefaultAssay(fib_P60) <- "Spatial"
cellAnnota_P60_fib <- subset(fib_P60@meta.data, select='celltype.case')
exprMatrix_P60_fib <- as.matrix(GetAssayData(fib_P60[["Spatial"]], slot='counts'))

cellAnnota_P60 <-rbind(cellAnnota_P60_Epi,cellAnnota_P60_fib)
exprMatrix_P60<-cbind(exprMatrix_P60_Epi, exprMatrix_P60_fib)
exprMatrix_P60 <- subset(exprMatrix_P60,select=rownames(cellAnnota_P60))


#for infercnv
write.table(exprMatrix_P60, paste('P60_P60','_exprMatrix.txt',sep=""), col.names=NA, sep='\t',quote = F)
write.table(cellAnnota_P60, paste('P60_P60','_cellAnnota.txt',sep=""), col.names=F, sep='\t',quote = F)
```
#Run inferCNV
```R
library(infercnv)
library(dplyr)


##Create inferCNV object
rm(list=ls())
options(stringsAsFactors = F)
expFile='P60_P60_exprMatrix.txt'
groupFiles='P60_P60_cellAnnota.txt'

options(scipen = 100)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= '/path/to/Spa-PR/PR/P60-loupe/3_inferCNV/hg38_gencode_v27.txt',
				                    ref_group_names=c("7")
					)

save(infercnv_obj,file='infercnv_P60.rda')
#perform infercnv operations to reveal cnv signal
#load("infercnv_P60.rda")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='./output_dir_P60_new',
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
df<-read.table("/path/to/Spa-PR/PR/P60-loupe/3_inferCNV/output_dir_P60_new/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",header=T)
rownames(df)<-df$cell
df$cluster<-"Fib"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.1.1"),]$cluster<-"D"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.1.2"),]$cluster<-"E"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.2.1"),]$cluster<-"G"
df[which(df$cell_group_name=="all_observations.all_observations.1.1.2.2"),]$cluster<-"H"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.1.1"),]$cluster<-"K"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.1.2"),]$cluster<-"L"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.2.1"),]$cluster<-"O"
df[which(df$cell_group_name=="all_observations.all_observations.1.2.2.2"),]$cluster<-"P"


P60 <- AddMetaData(P60, metadata = df$cluster, col.name = "CNV_cluster")
table(P60$Sprod_snn_res.0.4,P60$CNV_cluster)
tb<-table(P60$Sprod_snn_res.0.4,P60$CNV_cluster)
tb<-as.data.frame(tb)
tb$M<-tb$all_observations.all_observations.1.2.2.2
p<-SpatialDimPlot(P60,group.by="CNV_cluster",label=T,label.size=3,repel=T)+scale_fill_manual(values=c("'#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4'))

# Annotate genes on the heatmap
dat<-read.table("/path/to/Spa-PR/PR/P60-loupe/3_inferCNV/output_dir_P60_new/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.pred_cnv_genes.dat",header=T)
rownames(df)<-df$cell
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
infercnv.dend <- read.dendrogram(file = "/path/to/Spa-PR/PR/P60-loupe/3_inferCNV/output_dir_P60/infercnv.observations_dendrogram.txt")
# Cut tree 
infercnv.labels <- cutree(infercnv.dend, k = 5, order_clusters_as_data = FALSE)
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
# Apply the colours of the dendrogram to the tissue section as well
rownames(the_bars)<-names(infercnv.labels)
P60_no_Fib <- AddMetaData(P60_no_Fib, metadata = the_bars, col.name = "CNV_color")
table(P60_no_Fib$CNV_cluster,P60_no_Fib$CNV_color)
p<-SpatialDimPlot(P60_no_Fib,group.by="CNV_cluster",label=T,label.size=3,repel=T)+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#5C4033"))
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
