```R
library(Seurat)
load("/path/to/Spa-PR/PR/P4/1_Seurat-cluster-and-anno/PR_P4_cluster_and_anno.RData")
exprSet<-PR_P4@assays$SCT@data

```
```R
dat=exprSet
library(estimate)

pro='PR_P4'
input.f=paste0(pro,'_estimate_input.txt')
output.f=paste0(pro,'_estimate_gene.gct')
output.ds=paste0(pro,'_estimate_score.gct')
write.table(dat,file = input.f,sep = '\t',quote = F)

filterCommonGenes(input.f=input.f,
                  output.f=output.f ,
                  id="GeneSymbol")
estimateScore(input.ds = output.f,
              output.ds=output.ds,
              platform="illumina")   ## note the 'platform' argument
scores=read.table(output.ds,skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores<-as.data.frame(scores)
scores$Tumour_purity<-cos(0.6049872018+0.0001467884*scores$ESTIMATEScore)
```
# Map tumour purity, stromal score and immune score onto the tissue section
```R
library(ggplot2)
# Tumour purity
PR_P4[['Tumour_purity']] <- scores$Tumour_purity
PR_P4[['Stromal_Score']] <- scores$StromalScore
PR_P4[['Immune_Score']] <- scores$ImmuneScore


p1<-SpatialFeaturePlot(PR_P4, features = c("Tumour_purity"))
p2<-SpatialFeaturePlot(PR_P4, features = c("Stromal_Score"))
p3<-SpatialFeaturePlot(PR_P4, features = c("Immune_Score"))
plot<-p1+p2+p3

ggsave(plot,file="PR_P4_tumor_zone_estimate_spatialplot.pdf",width=15,height=5)
save(PR_P4,file="/path/to/Spa-PR/PR/P4/1_Seurat-cluster-and-anno/PR_P4_cluster_and_anno.RData")
```
# Tumour purity of each Sprod_cluster (violin plot)
```R
# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(showtext)
font_add('Arial','/Library/Fonts/Arial.ttf') # load the font; on macOS fonts live in /Library/Fonts
showtext_auto() # showtext must be enabled or the font is lost in ggsave(), which opens and closes the graphics device

data<-as.data.frame(PR_P4$Sprod_snn_res.0.6)
data$Tumor_purity<-PR_P4$Tumour_purity
colnames(data)<-c("Cluster","Score")

cluster_anno <- data.frame(
  name=c( 0:24  ),
  value=c("GS4_with_GS5","Normal_with_undefined_Tumor","GS3","GS3","GS3",
          "GS4_forming_ductal_carcinoma","Normal_with_Tumor","Normal",
          "GS4_forming_ductal_carcinoma","Normal","GS4_and_GS5",
          "GS3","GS5_3","GS4_forming_ductal_carcinoma","GS3","Normal",
          "Normal_with_immune","GS4_with_GS5","GS4_2","GS5_1","Normal_with_Tumor",
          "GS4_1","Normal_with_Tumor","GS5_2","Normal_with_Tumor"
          )
)
data$Cluster<-as.character(data$Cluster)
cluster_anno$name<-as.character(cluster_anno$name)
colnames(cluster_anno)<-c("Cluster","Anno")
p<-data %>%
  left_join(cluster_anno) %>%
  mutate(myaxis = paste0(Cluster, "\n", Anno)) %>%
  ggplot( aes(x=myaxis, y=Score, fill=Cluster)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("PR_P4_Sprod_cluster_tumor_purity_violinplot") +
  xlab("")
p<-p + theme(axis.text.x = element_text(size = 10, family = "myFont", color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 75))
ggsave(p,file="PR_P4_tumor_purity_groupby_sprod_cluster_violinplot.pdf",width=15)
```
# Tumour purity of each anno_cluster (violin plot)
```R
data<-as.data.frame(PR_P4$anno_cluster)
data$Tumor_purity<-PR_P4$Tumour_purity
colnames(data)<-c("Cluster","Score")
sample_size = data %>% group_by(Cluster) %>% summarize(num=n())

p<-data %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(Cluster,"\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=Score, fill=Cluster))+
  geom_violin(width=1.4) +
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#5C4033"))+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  #scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("PR_P4_Anno_cluster_tumor_purity_violinplot") +
  xlab("")
p<-p + theme(axis.text.x = element_text(size = 10, family = "myFont", color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 45))
ggsave(p,file="PR_P4_tumor_purity_groupby_anno_cluster_violinplot.pdf",width=15)
```
