# Load the integrated clustering data
```R
library(Seurat)
library(dplyr)
library(ggplot2)

load("/path/to/Spa-PR/PR/intergretd_analysis/intergrated_sample/sc_anno_by_Spa/integrated_prad_estimate_results.RData")
Idents(prad)<-prad$integrated_snn_res.0.4
```
# Compute transcriptome diversity
```R
library(matrixStats)
tumor_region_spots<-prad@assays$integrated@data
tumor_region_spots<-as.data.frame(tumor_region_spots)
# Compute the Pearson correlation coefficients
correlation_matrix <- cor(tumor_region_spots[VariableFeatures(prad),])

# Compute the median absolute deviation (MAD) of the Pearson correlations
mad_correlation<-c()
for (i in colnames(correlation_matrix)) {
  cor<-correlation_matrix[,i]
  cor<-cor[which(names(cor) != i)]
  cor_median<-median(cor)
  cor<-abs(cor-cor_median)
  mad_correlation<-c(mad_correlation,median(cor))
}

# Compute the transcriptome diversity degree
transcriptome_diversity_degree <- 1.4826 * mad_correlation
names(transcriptome_diversity_degree)<-colnames(correlation_matrix)
save(transcriptome_diversity_degree,file="transcriptome_diversity_degree.RData")
```
# Compute spatial continuity
```R
library(dplyr)
# Assume 'spots' is a data frame holding the coordinates and other information of every spot
# Assume 'x' and 'y' are the names of the coordinate columns
find_neighboring_spots <- function(spots, target_spot, distance_threshold) {
  neighbors <- spots %>%
    filter(sqrt((x - target_spot$x)^2 + (y - target_spot$y)^2) <= distance_threshold) %>%
    filter((y - target_spot$y)^2 <4)
  return(neighbors)
}
distance_threshold <- 2  # distance threshold that defines the neighbourhood radius
spatial_continuity<-c()
spot_name<-c()
for (i in names(table(prad$orig.ident))) {
  PR_epi<-subset(prad,subset=orig.ident==i)
  df_position<-data.frame(id=colnames(PR_epi),row=PR_epi@images$slice1@coordinates$row,col=PR_epi@images$slice1@coordinates$col)
  rownames(df_position)<-df_position$id
  df_position<-df_position[,c(2,3)]
  colnames(df_position)<-c("x","y")
  df_label<-PR_epi$integrated_snn_res.0.4
  for (j in rownames(df_position)) {
    target_spot <- df_position[j, ]
    neighboring_spots <- find_neighboring_spots(df_position, target_spot, distance_threshold)
    neighboring_spots_label<-df_label[rownames(neighboring_spots)]
    count <- (length(neighboring_spots_label[neighboring_spots_label == df_label[j]])-1)/(length(neighboring_spots_label)-1)
    spatial_continuity<-c(spatial_continuity,count)
  }
  spot_name<-c(spot_name,rownames(df_position))
}
names(spatial_continuity)<-spot_name
save(spatial_continuity,file="spatial_continuity.RData")
```
# Scatter plot
```R
prad$transcriptome_diversity_degree<-1
prad$transcriptome_diversity_degree[names(transcriptome_diversity_degree)]<-transcriptome_diversity_degree

prad$spatial_continuity<-1
prad$spatial_continuity[names(spatial_continuity)]<-spatial_continuity

df<-data.frame(id=colnames(prad),patient=prad$orig.ident,cluster=prad$integrated_snn_res.0.4,transcriptome_diversity_degree=prad$transcriptome_diversity_degree,spatial_continuity=prad$spatial_continuity)

df[which(df$patient=="P60"),]$patient<-"P6"
df[which(df$patient=="P616"),]$patient<-"P7"
df[which(df$patient=="PR_P1"),]$patient<-"P1"
df[which(df$patient=="PR_P2"),]$patient<-"P2"
df[which(df$patient=="PR_P3"),]$patient<-"P3"
df[which(df$patient=="PR_P4"),]$patient<-"P4"
df[which(df$patient=="PR_P5"),]$patient<-"P5"
df$new_cluster<-paste(df$patient,df$cluster,sep="_")

library(tidyverse)

df_2 <- df %>%
  group_by(cluster) %>%
  summarize(transcriptome_diversity_degree=mean(transcriptome_diversity_degree),spatial_continuity=mean(spatial_continuity)) 

df_2<-as.data.frame(df_2)

library(ggplot2)
library(hrbrthemes)
allcolour=c("0"="#dd3497",
            "1"="#f768a1",
            "2"="#ae017e",
            "3"="#78c679",
            "4"="#d9f0a3",
            "5"="#fb9a99",
            "6"="#ff7f00",
            "7"="#1f78b4",
            "8"="#b15928",
            "9"="#7a0177",
            "10"="#e31a1c",
            "11"="#006837",
            "12"="#49006a")
p<-ggplot(df_2, aes(x=spatial_continuity, y=transcriptome_diversity_degree, color=cluster)) + 
    geom_point(size=6) +
    scale_color_manual(values=allcolour)+
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5)) + 
    theme(panel.grid = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) 



ggsave(p,file="test.pdf")


```
```R
df_2 <- df %>%
  group_by(new_cluster) %>%
  summarize(transcriptome_diversity_degree=mean(transcriptome_diversity_degree),spatial_continuity=mean(spatial_continuity)) 
df_2$patient<-sapply(strsplit(df_2$new_cluster, "_"), function(x) x[1])
df_2$cluster<-sapply(strsplit(df_2$new_cluster, "_"), function(x) x[2])
p <- ggplot(df, aes(x = A, y = B, shape = C, color = D)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c(1, 2, 3)) +
  scale_color_manual(values = c("red", "blue", "green"))
```
