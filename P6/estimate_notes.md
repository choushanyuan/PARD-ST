```R
library(Seurat)
load("/path/to/Spa-PR/PR/P60-loupe/1_Seurat-cluster-and-anno/P60_CNV_cluster.RData")
exprSet<-P60_no_Fib@assays$SCT@data

```
```R
dat=exprSet
library(estimate)

pro='P60'
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
# Tumour purity
P60_no_Fib[['Tumour_purity']] <- scores$Tumour_purity
P60_no_Fib[['Stromal_Score']] <- scores$StromalScore
P60_no_Fib[['Immune_Score']] <- scores$ImmuneScore

Idents(P60_no_Fib)<-"CNV_cluster"

P60_Tumor<- subset(P60_no_Fib,subset=CNV_cluster %in% c("5","3","2"))


p1<-SpatialFeaturePlot(P60_Tumor, features = c("Tumour_purity"))
p2<-SpatialFeaturePlot(P60_Tumor, features = c("Stromal_Score"))
p3<-SpatialFeaturePlot(P60_Tumor, features = c("Immune_Score"))
plot<-p1+p2+p3

ggsave(plot,file="P60_tumor_zone_estimate_spatialplot.pdf",width=15,height=5)
save(P60_no_Fib,file="/path/to/Spa-PR/PR/P60-loupe/1_Seurat-cluster-and-anno/P60_CNV_cluster.RData")
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

data<-as.data.frame(P60$Sprod_snn_res.0.4)
data$Tumor_purity<-P60$Tumour_purity
colnames(data)<-c("Cluster","Score")

cluster_anno <- data.frame(
  name=c( 0:13  ),
  value=c("GS3_2","GS3_2","GS3_2","GS3_and_GS4_1","GS4","GS3_2","GS3_1","Fibroblasts",
          "GS3_and_GS4_2","GS3_and_GS4_2","GS4","GS3_and_GS4_1","GS4","GS3_and_GS4_2"
          
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
  ggtitle("P60_Sprod_cluster_tumor_purity_violinplot") +
  xlab("")
p<-p + theme(axis.text.x = element_text(size = 10, family = "myFont", color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 75))
ggsave(p,file="P60_tumor_purity_groupby_sprod_cluster_violinplot.pdf",width=15)
```
# Tumour purity of each anno_cluster (violin plot)
```R
data<-as.data.frame(P60$anno_cluster)
data$Tumor_purity<-P60$Tumour_purity
colnames(data)<-c("Cluster","Score")
sample_size = data %>% group_by(Cluster) %>% summarize(num=n())

p<-data %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(Cluster,"\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=Score, fill=Cluster))+
  geom_violin(width=1.4) +
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c"))+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  #scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("P60_Anno_cluster_tumor_purity_violinplot") +
  xlab("")
p<-p + theme(axis.text.x = element_text(size = 10, family = "myFont", color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 45))
ggsave(p,file="P60_tumor_purity_groupby_anno_cluster_violinplot.pdf",width=15)
```
