#/path/to/Spa-PR/sRCC/1_Seurat-cluster-and-anno/sRCC_P1_cluster_and_anno.RData
# sRCC_02: /path/to/Spa-PR/sRCC_02/1_Seurat-cluster-and-anno/sRCC_02-cluster-results.RData
.libPaths(c(.libPaths(),"/path/to/R/x86_64-conda_cos6-linux-gnu-library/4.0","/path/to/.conda/envs/env/lib/R/library"))

library(Seurat)
library(tidyverse)
library(scales)
library(pals)
library(cowplot)
library(fgsea)
library(msigdbr)
library(data.table)


load("/path/to/Spa-PR/PR/P60-loupe/1_Seurat-cluster-and-anno/P60_denoised_cluster.RData")
load("/path/to/Spa-PR/PR/P616-loupe-much/1_Seurat-cluster-and-anno/P616_denoised_cluster.RData")
load("/path/to/Spa-PR/PR/P1/1_Seurat-cluster-and-anno/PR_P1_cluster_and_anno.RData")
load("/path/to/Spa-PR/PR/P2/1_Seurat-cluster-and-anno/PR_P2_cluster_and_anno.RData")
load("/path/to/Spa-PR/PR/P3/1_Seurat-cluster-and-anno/PR_P3_cluster_and_anno.RData")
load("/path/to/Spa-PR/PR/P4/1_Seurat-cluster-and-anno/PR_P4_cluster_and_anno.RData")
load("/path/to/Spa-PR/PR/P5/1_Seurat-cluster-and-anno/PR_P5_cluster_and_anno.RData")
options(future.globals.maxSize = 20000 * 1024^2)

# Run SCTransform integration workflow: https://satijalab.org/seurat/archive/v3.0/integration.html

prad.list <- list("P60" = P60, "P616" = P616, "PR_P1" = PR_P1, "PR_P2" = PR_P2, "PR_P3" = PR_P3, "PR_P4" = PR_P4, "PR_P5" = PR_P5)


genes.common <- Reduce(intersect, list(rownames(P60), rownames(P616), rownames(PR_P1), rownames(PR_P2), rownames(PR_P3), rownames(PR_P4), rownames(PR_P5)))
list.features <- SelectIntegrationFeatures(object.list = prad.list,
                                           nfeatures = 3000,
                                           assay = c("SCT", "SCT", "SCT", "SCT", "SCT", "SCT", "SCT"))

prad.list <- PrepSCTIntegration(object.list = prad.list,
                               anchor.features = list.features,
                               assay = "SCT",
                               verbose = F)
prad.anchors <- FindIntegrationAnchors(object.list = prad.list,
                                      normalization.method = "SCT",
                                      anchor.features = list.features,
                                      verbose = F)                       
prad <- IntegrateData(anchorset = prad.anchors,
                     features.to.integrate = genes.common,
                     normalization.method = "SCT", 
                     verbose = F)


save(prad, file = "prad_merged_nonPCA.rda")
# Rerun dimensionality reduction and clustering on integrated object.
prad <- RunPCA(prad, npcs = 30, verbose = F)
prad <- FindNeighbors(prad, reduction = "pca", dims = 1:30)
prad <- FindClusters(prad, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6), verbose = F)
prad <- RunUMAP(prad, reduction = "pca", dims = 1:30)

save(prad, file="prad_merged_integrated.rda")

#
Idents(prad)<-"integrated_snn_res.0.2"
prad@misc$markers <- FindAllMarkers(object = prad, assay = 'Sprod',only.pos = TRUE, test.use = 'wilcox')
write.table(prad@misc$markers,file='prad_resolution0.2.integrated.txt',row.names = FALSE,quote = FALSE,sep = '\t')
top30 <- prad@misc$markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(top30,"prad_resolution0.2.integrated_top30_markers.csv")



p0 <- DimPlot(prad, reduction = "pca", group.by = c("integrated_snn_res.0.2", "orig.ident"))
ggsave(p0, filename = "merge.cluster.show.pdf",height = 5, width = 12)

p1 <- SpatialDimPlot(prad,group.by = c("integrated_snn_res.0.2"),ncol=4)
ggsave(p1, filename = "merge.Spatial.dimplot.pdf",height = 8, width = 20)

p2 <- SpatialFeaturePlot(prad, features = c("KRT5", "COL4A1")) 
ggsave(p2, filename = "merge.SpatialFeature.pdf",height = 8, width = 16)

p2 <- SpatialFeaturePlot(prad, features = c("KRT5", "COL4A1")) 
# Change the colours
for(i in 1:14)
{
    #p2[[i]] <- p2[[i]]+scale_fill_gradientn(colours = c("lightgray", "mistyrose", "red", "darkred", "black"))
    p2[[i]] <- p2[[i]]+scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11,'PiYG'))
} 
ggsave(p2, filename = "merge.SpatialFeature.col.pdf",height = 8, width = 26)



save(prad, file = "PRAD_merged.rda")


