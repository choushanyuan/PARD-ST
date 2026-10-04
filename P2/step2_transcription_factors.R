###################################################
###################step2###########################
###################################################
#### Transcription factors ####
library(dorothea)
library(Seurat)
library(pheatmap)
library(tibble)
library(dplyr)
library(viper)
library(tidyr)
library(ggplot2)

dorothea_regulon_human<-get(data("dorothea_hs",package="dorothea"))
regulon<-dorothea_regulon_human %>% dplyr::filter(confidence %in% c("A","B","C"))

load("/path/to/Spa-PR/PR/P2/8_part1_analysis/ST_obj_new.RData")
ST_obj_tumor<-subset(ST_obj,subset=purity_cluster %in% c("Tumor_cured","Tumor_GS4"))
load("/path/to/Spa-PR/PR/P2/8_part1_analysis/infercnv.labels.RData")
infercnv.labels[which(infercnv.labels==1)]<-"subclone_B"
infercnv.labels[which(infercnv.labels==2)]<-"subclone_A"
ST_obj_tumor <- AddMetaData(ST_obj_tumor, metadata = infercnv.labels, col.name = "CNV_cluster")
ST_obj_tumor <-subset(ST_obj_tumor,subset=CNV_cluster!="NA")
ST_obj<-ST_obj_tumor

DefaultAssay(ST_obj)<-"Spatial"
ST_obj<-NormalizeData(ST_obj)
all.genes<-rownames(ST_obj)
ST_obj<-ScaleData(ST_obj,features=all.genes)

load("cell_state.RData")
ST_obj[['State']]<-0
ST_obj[['State']][rownames(cell_state),] <- cell_state
Idents(ST_obj)<-"State"
ST_obj@assays$RNA<-ST_obj@assays$Spatial
ST_obj<-run_viper(ST_obj,regulon,options=list(method="scale",minsize=4,eset.filter=F,cores=1,verbose=F))

DefaultAssay(ST_obj)<-"dorothea"
ST_obj<-ScaleData(ST_obj)

viper_scores_df<-GetAssayData(ST_obj,slot="scale.data",assay="dorothea") %>% data.frame(check.names=F) %>% t()
CellsClusters<-data.frame(cell=names(ST_obj$State),cell_type=as.character(ST_obj$State),check.names=F)
viper_scores_clusters<-viper_scores_df %>% data.frame() %>% rownames_to_column("cell") %>% gather(tf,activity,-cell) %>% inner_join(CellsClusters)
summarized_viper_scores<-viper_scores_clusters %>% group_by(tf,cell_type) %>% summarise(avg=mean(activity),std=sd(activity))
save(viper_scores_df,viper_scores_clusters,summarized_viper_scores,file="TF_results.RData")

# Run locally
library(tidyverse)
load("TF_results.RData")
highly_variable_tfs<-summarized_viper_scores %>% 
  group_by(tf) %>% 
  mutate(var=var(avg)) %>% 
  ungroup() %>% 
  top_n(350,var) %>% 
  distinct(tf)#50x7
summarized_viper_scores_df<-summarized_viper_scores %>% 
  semi_join(highly_variable_tfs,by="tf") %>% 
  dplyr::select(-std) %>% 
  spread(tf,avg) %>% 
  data.frame(row.names=1,check.names=F)


library(pheatmap)
df<-t(summarized_viper_scores_df)[,c("1","2","3","4","5","6","7")]
p<-pheatmap(df,fontsize=12,fontsize_row=10,cluster_rows = T,cluster_cols = T,
            color=colorRampPalette(c("#0571b0","#f7f7f7","#ca0020"))(100),
            show_colnames = T,border_color = "grey",scale = "row",show_rownames =T,
            angle_col = 90,
            cutree_col=2,
            cutree_rows=2,
            cellwidth = 15, cellheight = 10)
p
ggsave(p,file="TF.pdf",width = 5,height = 9)
