# Add a helper function
get_fraction <- function(HSMM_myo) {
  data_df <- t(reducedDimS(HSMM_myo)) %>%
    as.data.frame() %>%
    select_('Component 1' = 1, 'Component 2' = 2) %>%
    rownames_to_column("Cells") %>%
    mutate(pData(HSMM_myo)$State,
           
           pData(HSMM_myo)$Pseudotime,
           
           pData(HSMM_myo)$orig.ident,
           
           pData(HSMM_myo)$clusters)
  
  
  colnames(data_df) <- c("cells","Component_1","Component_2","State",
                         "Pseudotime","orig.ident","celltype")
  
  Cellratio <- prop.table(table(data_df$State, data_df$celltype), margin = 2) # proportion of each cell group per sample
  Cellratio <- as.data.frame(Cellratio)
  colnames(Cellratio) <- c('State',"celltype","Freq")
  return(Cellratio)
}
###################################################
###################step1###########################
###################################################
library(Seurat)
library(monocle)
load("/path/to/Spa-PR/PR/P3/12_part1_analysis_2/ST_obj_new.RData")
my_color<-c("PIN"="#fbb4b9","GS3"="#f768a1","GS4"="#c51b8a","GS4+5"="#7a0177")

PR_epi<-subset(ST_obj,subset=anno_all %in% c("PIN","GS4","GS4+5","GS3"))

PR_epi<- SCTransform(PR_epi, assay = "Spatial", verbose = FALSE)
PR_epi<- RunPCA(PR_epi, assay = "SCT", verbose = FALSE, dims = 1:30)
expression_matrix = PR_epi@assays$Spatial@counts
cell_metadata  <- data.frame(group = PR_epi[['orig.ident']],clusters = PR_epi@meta.data$anno_all)
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), stringsAsFactors = F) #gene_id=rownames(expression_matrix)
rownames(gene_annotation) <- rownames(expression_matrix)
pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(expression_matrix,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily=negbinomial.size())
# Estimate size factors
HSMM_myo <- estimateSizeFactors(HSMM)
# Estimate dispersions
HSMM_myo <- estimateDispersions(HSMM_myo)
# Take the highly variable genes selected by Seurat
expressed_genes<-VariableFeatures(PR_epi)
# Differential test
diff_test_res1 <- differentialGeneTest(HSMM_myo[expressed_genes,],fullModelFormulaStr = '~clusters', cores = 4)
# Select the differentially expressed genes
ordering_genes <- subset(diff_test_res1, qval < 0.01)
ordering_genes <-ordering_genes[order(ordering_genes$qval,decreasing=F),]
write.table(ordering_genes,file="train_monocles_DEG_gene.txt",col.names=T,row.names=F,sep="\t",quote=F)
ordergene<-rownames(ordering_genes)
# Filter genes
HSMM_myo <- setOrderingFilter(HSMM_myo, ordergene)
# Reduce the dimension
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, method = 'DDRTree')
# Order the cells
HSMM_myo <- orderCells(HSMM_myo)
#HSMM_myo <- orderCells(HSMM_myo,root_state=3)

# Inspect the differentiation relationship between the subpopulations
HSMM_myo$clusters<-as.character(HSMM_myo$clusters)

pdf("P3_monocle_clusters.pdf")
plot_cell_trajectory(HSMM_myo, color_by = "clusters",show_branch_points = T,cell_size =1.5)+scale_color_manual(values=my_color)
dev.off()
pdf("P3_monocle_State.pdf")
plot_cell_trajectory(HSMM_myo, color_by = "State",show_branch_points = T,cell_size =1.5)
dev.off()
#HSMM_myo_2<-HSMM_myo
#HSMM_myo_2$Pseudotime<-max(HSMM_myo_2$Pseudotime)-HSMM_myo_2$Pseudotime
pdf("P3_monocle_pseudotime_plot.pdf")
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime",show_branch_points = T,cell_size =1.5)+scale_colour_distiller(palette = "Spectral")
dev.off()
# Split the plot by subpopulation
pdf("P3_monocle_clusters_split.pdf")
plot_cell_trajectory(HSMM_myo, color_by = "clusters",show_branch_points = T,cell_size =1.5) + facet_wrap(~clusters, nrow = 5)+scale_color_manual(values=my_color)
dev.off()
# Map the (inferred) differentiation time back onto the tissue section
cell_Pseudotime <- data.frame(pData(HSMM_myo)$Pseudotime)
rownames(cell_Pseudotime) <- rownames(cell_metadata)

ST_obj[['Pseudotime']] <- NA
ST_obj[['Pseudotime']][rownames(cell_Pseudotime),] <- cell_Pseudotime
pdf("P3_monocle_pseudotime_spatialfeatureplot.pdf")
SpatialFeaturePlot(ST_obj, features = c("Pseudotime"),pt.size.factor = 1.5)
dev.off()
# Map the monocle state back onto the tissue section
cell_state <- data.frame(pData(HSMM_myo)$State)
rownames(cell_state) <- rownames(cell_metadata)
ST_obj[['State']] <- 0
ST_obj[['State']][rownames(cell_state),] <- cell_state
save(cell_state,file="cell_state.RData")
Cellratio<-get_fraction(HSMM_myo)
save(Cellratio,file = "Cellratio.RData")



