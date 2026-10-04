library(CellChat)
library(Seurat)
library(tidyverse)
library(viridis)
library(RColorBrewer)


#P2
load("integrated_label.RData")
load("/path/to/Spa-PR/PR/P2/8_part1_analysis/ST_obj_new.RData")
ST_obj[['barcode']] <- paste(colnames(ST_obj), "_4", sep = "")

df_label_2<-df_label[which(df_label$barcode %in% ST_obj$barcode),]
rownames(df_label_2)<-df_label_2$barcode
df_label_2<-df_label_2[ST_obj$barcode,]
ST_obj[['label']] <-"temp"
df_label_3<-data.frame(label_1=colnames(ST_obj),label_2=ST_obj$barcode)
rownames(df_label_3)<-df_label_3$label_2
ST_obj[['label']][df_label_3$label_1,] <- df_label_2$lable

DefaultAssay(ST_obj)<-"Spatial"
ST_obj<-NormalizeData(ST_obj)
all.genes<-rownames(ST_obj)
ST_obj<-ScaleData(ST_obj,features=all.genes)
ST_obj$label<-factor(ST_obj$label,levels=c("PIN-like","Subtype-low","Subtype-medium","Subtype-high",
                                           "Cluster_3","Cluster_4","Cluster_6","Cluster_7","Cluster_8")) # differs per patient
Idents(ST_obj)<-ST_obj$label


# Expression matrix
data.input = Seurat::GetAssayData(ST_obj, slot = "data", assay = "Spatial") 
# Meta information
meta = data.frame(labels = Idents(ST_obj), # custom name
                  row.names = names(Idents(ST_obj))) # manually create a dataframe consisting of the cell labels
unique(meta$labels)

# Spatial image information
spatial.locs = Seurat::GetTissueCoordinates(ST_obj, scale = NULL, 
                                            cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameter information
scale.factors = jsonlite::fromJSON(txt = 
                                     file.path("/path/to/Spa-PR/PR/P2/0_PR_02-loupe-Spaceranger-results/PR_P2_Spaceranger_for_sprod/outs/spatial", 'scalefactors_json.json'))
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)
# Build the CellChat object
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "labels", # the name defined in the meta data above is "labels"
                           datatype = "spatial", ###
                           coordinates = spatial.locs, 
                           scale.factors = scale.factors)

# Set the reference database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# Pre-processing
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) # use 1 on a laptop
## Identify over-expressed genes
cellchat <- identifyOverExpressedGenes(cellchat)
# Identify over-expressed ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat)

# Infer the cell-cell communication network
cellchat <- computeCommunProb(cellchat, 
                              #type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, 
                              scale.distance = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Note 1: 'type' defaults to triMean, producing fewer but stronger interactions; 'type = "truncatedMean"' needs the 'trim' argument and produces more interactions.
# Note 2: distance.use = FALSE filters out interactions between spots that are far apart in space.
# Note 3: the spatial default is the truncatedMean + trim combination; tune it empirically.
# Compute the communication of all ligand-receptor pairs belonging to each signalling pathway
cellchat <- computeCommunProbPathway(cellchat)
# Aggregate the communication across cell types
cellchat <- aggregateNet(cellchat)
save(cellchat,file="P2_cellchat_results.RData")
#### Visualisation
my_color<-c("PIN-like"="#e31a1c",
            "Subtype-low"="#f768a1",
            "Subtype-medium"="#ae017e",
            "Subtype-high"="#49006a",
            "Cluster_3"="#78c679",
            "Cluster_4"="#d9f0a3",
            "Cluster_6"="#ff7f00",
            "Cluster_7"="#1f78b4",
            "Cluster_8"="#b15928"
            #"Cluster_11"="#006837",
            )
# Communication between cell types (circle plot)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("P2_aggregated_cell_cell_communication_network.pdf")
netVisual_circle(cellchat@net$count, color.use=my_color,vertex.weight = rowSums(cellchat@net$count), 
                 weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,color.use=my_color,vertex.weight = rowSums(cellchat@net$weight), 
                 weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
dev.off()
# Heatmap
pdf("P2_heatmap.pdf")
netVisual_heatmap(cellchat, color.use=my_color, measure = "count", color.heatmap = c("Blues"))
netVisual_heatmap(cellchat, color.use=my_color, measure = "weight", color.heatmap = c("Blues"))
dev.off()

# Find the pathways with the strongest communication
cellchat<-computeCommunProbPathway(cellchat)
df.netp<-subsetCommunication(cellchat,slot.name="netP")
# Define a helper function
get_top_10 <- function(df.netp, filename) {
  pathway_name<-names(table(df.netp$pathway_name))
  tmp<-data.frame(pathway_name=pathway_name,mean=NA)
  rownames(tmp)<-tmp$pathway_name
  for (i in pathway_name) {
    pathway_sum<-df.netp[which(df.netp$pathway_name==i),]
    tmp[i,2]<-mean(pathway_sum$prob)
  }
  tmp_order<-tmp[order(-tmp$mean),]
  top_10_pathway<-df.netp[which(df.netp$pathway_name %in% rownames(head(tmp_order,10))),]
  top_10_pathway$pathway_name<-factor(top_10_pathway$pathway_name,levels=rownames(head(tmp_order,10)))
  
  # Bottom Right
  p<-ggplot(top_10_pathway, aes(x=pathway_name, y=prob, fill=pathway_name)) + 
    geom_boxplot() +
    theme(legend.position="none") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5)) + # black panel border
    theme(panel.background = element_rect(fill = "white")) + # pure white background
    scale_fill_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd"))
  
  ggsave(p,file=filename)
}
# Top 10 pathways between tumour and non-tumour
df.netp_with_non_tumor<-df.netp[which((df.netp$source %in% c("PIN-like","Subtype-low","Subtype-medium","Subtype-high") & df.netp$target %in% c("Cluster_3","Cluster_4","Cluster_6","Cluster_7","Cluster_8","Cluster_11")) | (df.netp$target %in% c("PIN-like","Subtype-low","Subtype-medium","Subtype-high") & df.netp$source %in% c("Cluster_3","Cluster_4","Cluster_6","Cluster_7","Cluster_8","Cluster_11"))),]
get_top_10(df.netp_with_non_tumor,"P2_Top_10_with_non_tumor_pathway.pdf")
# Top 10 pathways within the tumour compartment only
df.netp_only_tumor<-df.netp[which(df.netp$source %in% c("PIN-like","Subtype-low","Subtype-medium","Subtype-high") & df.netp$target %in% c("PIN-like","Subtype-low","Subtype-medium","Subtype-high")),]
get_top_10(df.netp_only_tumor,"P2_Top_10_only_tumor_pathway.pdf")
# Show individual signalling pathways
cellchat@netP$pathways
pathways.show1 <- c("COLLAGEN")
pathways.show2<- c("LAMININ")
pathways.show3 <- c("FN1")
pathways.show4 <- c("MK")
pathways.show5 <- c("MIF")
pathways.show6 <- c("APP")
pathways.show7 <- c("RA")


levels(cellchat@idents)   
vertex.receiver = c(1:4)  # indices into levels(cellchat@idents)
pdf("P2_pathway.pdf",width=14)
netVisual_aggregate(cellchat, signaling = pathways.show1, color.use=my_color,                     
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
netVisual_aggregate(cellchat, signaling = pathways.show2, color.use=my_color,                     
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
netVisual_aggregate(cellchat, signaling = pathways.show3, color.use=my_color,                     
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
netVisual_aggregate(cellchat, signaling = pathways.show4, color.use=my_color,                     
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
netVisual_aggregate(cellchat, signaling = pathways.show5, color.use=my_color,                     
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
netVisual_aggregate(cellchat, signaling = pathways.show6, color.use=my_color,                     
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
netVisual_aggregate(cellchat, signaling = pathways.show7, color.use=my_color,                     
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
dev.off()
# Circle plot in spatial coordinates
# Spatial plot
par(mfrow=c(1,1))
pdf("P2_pathway_ST.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show1, layout = "spatial",sources.use = c(5:9), targets.use = c(1:4), color.use=my_color,
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = pathways.show2, layout = "spatial", sources.use = c(5:9), targets.use = c(1:4),color.use=my_color,
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = pathways.show3, layout = "spatial", sources.use = c(5:9), targets.use = c(1:4),color.use=my_color,
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = pathways.show4, layout = "spatial", sources.use = c(5:9), targets.use = c(1:4),color.use=my_color,
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = pathways.show5, layout = "spatial", sources.use = c(5:9), targets.use = c(1:4),color.use=my_color,
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = pathways.show6, layout = "spatial", sources.use = c(5:9), targets.use = c(1:4),color.use=my_color,
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = pathways.show7, layout = "spatial", sources.use = c(5:9), targets.use = c(1:4),color.use=my_color,
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)
dev.off()
# Heatmap of the pathways
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
pdf("P2_pathway_heatmap.pdf")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show1, color.use=my_color,
                                  width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show2, color.use=my_color,
                                  width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show3, color.use=my_color,
                                  width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show4, color.use=my_color,
                                  width = 8, height = 2.5, font.size = 10)
dev.off()
# Ligand-receptor bubble plot
levels(cellchat@idents)
pdf("P2_bubble_source.pdf",height=10)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(5:9),                  
                 signaling = pathways.show1,title.name=pathways.show1,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(5:9),                  
                 signaling = pathways.show2,title.name=pathways.show2,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(5:9),                  
                 signaling = pathways.show3,title.name=pathways.show3,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(5:9),                  
                 signaling = pathways.show4,title.name=pathways.show4,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(5:9),                  
                 signaling = pathways.show5,title.name=pathways.show5,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(5:9),                  
                 signaling = pathways.show6,title.name=pathways.show6,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(5:9),                  
                 signaling = pathways.show7,title.name=pathways.show7,
                 remove.isolate = FALSE)
dev.off()
pdf("P2_bubble_target.pdf",height=10)
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4),                  
                 signaling = pathways.show1,title.name=pathways.show1,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4),                  
                 signaling = pathways.show2,title.name=pathways.show2,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4),                  
                 signaling = pathways.show3,title.name=pathways.show3,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4),                  
                 signaling = pathways.show4,title.name=pathways.show4,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4),                  
                 signaling = pathways.show5,title.name=pathways.show5,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4),                  
                 signaling = pathways.show6,title.name=pathways.show6,
                 remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4),                  
                 signaling = pathways.show7,title.name=pathways.show7,
                 remove.isolate = FALSE)
dev.off()
