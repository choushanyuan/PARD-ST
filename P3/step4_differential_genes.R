###################################################
###################step4###########################
###################################################
DefaultAssay(ST_obj)<-"Spatial"
ST_obj<-NormalizeData(ST_obj)
all.genes<-rownames(ST_obj)
ST_obj<-ScaleData(ST_obj,features=all.genes)
Idents(ST_obj)<-"State"

markers <- FindMarkers(ST_obj, ident.1 = "3", ident.2 = "1", only.pos =F) 
markers<-markers[which(markers$p_val_adj<0.05),]
markers<-markers[order(-markers$avg_log2FC),]
save(markers,file="P3_3_vs_1_DEG.RData")

markers <- FindMarkers(ST_obj, ident.1 = "2", ident.2 = "1", only.pos =F) 
markers<-markers[which(markers$p_val_adj<0.05),]
markers<-markers[order(-markers$avg_log2FC),]
save(markers,file="P3_2_vs_1_DEG.RData")

markers <- FindMarkers(ST_obj, ident.1 = "3", ident.2 = "2", only.pos =F) 
markers<-markers[which(markers$p_val_adj<0.05),]
markers<-markers[order(-markers$avg_log2FC),]
save(markers,file="P3_3_vs_2_DEG.RData")
###################################################
###################step5###########################
###################################################
# Volcano plot
# Read DEG data
load("P3_3_vs_1_DEG.RData")
x<-markers
x$gene<-rownames(x)

# Set thresholds
logFCcut <- 1#1 2
pvalCut <- 0.05 #0.05 0.001
# Set colours for up-/down-regulated genes
x[,7] <- ifelse((x$p_val < pvalCut & x$avg_log2FC > logFCcut), "red", ifelse((x$p_val < pvalCut & x$avg_log2FC < -logFCcut), "blue","grey30"))

# Set point size
size <- ifelse((x$p_val < pvalCut & abs(x$avg_log2FC) > logFCcut), 4, 2)

# Set the min/max of the x and y axes
head(x)
x_right<-3
tail(x)
x_left<--2.5
xmin <- (range(x$avg_log2FC)[1]-(range(x$avg_log2FC)[1]+5))
xmax <- (range(x$avg_log2FC)[1]+(5-range(x$avg_log2FC)[1]))
ymin <- 0
ymax <- range(-log(x$p_val_adj))[2]+1
ymax
ymax <- 150



# Construct the plot object
p1 <- ggplot(data=x, aes(x=avg_log2FC, y=-log10(p_val_adj), label = gene)) +
  geom_point(alpha = 0.6, size=size, colour=x[,7]) +
  scale_color_manual(values = c("lightgrey", "navy", "red")) + 
  
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value-adj")), title="P3_State-3_vs_State-1") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(x_left, -logFCcut, 0, logFCcut,x_right), # tick positions
    labels = c(x_left, -logFCcut, 0, logFCcut,x_right),
    limits = c(x_left, x_right) # x-axis range; symmetric looks better
  ) +
  # Draw the threshold lines
  geom_vline(xintercept = logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_vline(xintercept = -logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_hline(yintercept = -log10(pvalCut), color="grey40", linetype="longdash", size=0.5) +
  
  guides(colour = guide_legend(override.aes = list(shape=16)))+
  
  theme_bw(base_size = 12) +
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black", size=8),
        plot.title = element_text(hjust = 0.8),
        axis.text.x = element_text(face="bold", color="black", size=15),
        axis.text.y = element_text(face="bold",  color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=16),
        axis.title.y = element_text(face="bold",color="black", size=16))
p1

# Label the genes with |logFC| > 1 by type
p2<-p1 + geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val), 
                             label = ifelse((avg_log2FC > 1 & p_val < pvalCut)|(avg_log2FC < -1 & p_val < pvalCut), rownames(x),"")),
                         colour="black", size = 5, box.padding = unit(0.35, "lines"), 
                         point.padding = unit(0.3, "lines"))
p2

ggsave(p2,file="P3_3_vs_1_volcano.pdf", width=6,height=6)
########################################################
# Read DEG data
load("P3_2_vs_1_DEG.RData")
x<-markers
x$gene<-rownames(x)

# Set thresholds
logFCcut <- 1#1 2
pvalCut <- 0.05 #0.05 0.001
# Set colours for up-/down-regulated genes
x[,7] <- ifelse((x$p_val < pvalCut & x$avg_log2FC > logFCcut), "red", ifelse((x$p_val < pvalCut & x$avg_log2FC < -logFCcut), "blue","grey30"))

# Set point size
size <- ifelse((x$p_val < pvalCut & abs(x$avg_log2FC) > logFCcut), 4, 2)

# Set the min/max of the x and y axes
head(x)
x_right<-2
tail(x)
x_left<--2
xmin <- (range(x$avg_log2FC)[1]-(range(x$avg_log2FC)[1]+5))
xmax <- (range(x$avg_log2FC)[1]+(5-range(x$avg_log2FC)[1]))
ymin <- 0
ymax <- range(-log(x$p_val_adj))[2]+1
ymax
ymax <- 90



# Construct the plot object
p1 <- ggplot(data=x, aes(x=avg_log2FC, y=-log10(p_val_adj), label = gene)) +
  geom_point(alpha = 0.6, size=size, colour=x[,7]) +
  scale_color_manual(values = c("lightgrey", "navy", "red")) + 
  
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value-adj")), title="P3_State-2_vs_State-1") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(x_left, -logFCcut, 0, logFCcut,x_right), # tick positions
    labels = c(x_left, -logFCcut, 0, logFCcut,x_right),
    limits = c(x_left, x_right) # x-axis range; symmetric looks better
  ) +
  # Draw the threshold lines
  geom_vline(xintercept = logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_vline(xintercept = -logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_hline(yintercept = -log10(pvalCut), color="grey40", linetype="longdash", size=0.5) +
  
  guides(colour = guide_legend(override.aes = list(shape=16)))+
  
  theme_bw(base_size = 12) +
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black", size=8),
        plot.title = element_text(hjust = 0.8),
        axis.text.x = element_text(face="bold", color="black", size=15),
        axis.text.y = element_text(face="bold",  color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=16),
        axis.title.y = element_text(face="bold",color="black", size=16))
p1

# Label the genes with |logFC| > 1 by type
p2<-p1 + geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val), 
                             label = ifelse((avg_log2FC > 1 & p_val < pvalCut)|(avg_log2FC < -1 & p_val < pvalCut), rownames(x),"")),
                         colour="black", size = 5, box.padding = unit(0.35, "lines"), 
                         point.padding = unit(0.3, "lines"))
p2

ggsave(p2,file="P3_2_vs_1_volcano.pdf", width=6,height=6)
########################################################
# Read DEG data
load("P3_3_vs_2_DEG.RData")
x<-markers
x$gene<-rownames(x)

# Set thresholds
logFCcut <- 1#1 2
pvalCut <- 0.05 #0.05 0.001
# Set colours for up-/down-regulated genes
x[,7] <- ifelse((x$p_val < pvalCut & x$avg_log2FC > logFCcut), "red", ifelse((x$p_val < pvalCut & x$avg_log2FC < -logFCcut), "blue","grey30"))

# Set point size
size <- ifelse((x$p_val < pvalCut & abs(x$avg_log2FC) > logFCcut), 4, 2)

# Set the min/max of the x and y axes
head(x)
x_right<-2
tail(x)
x_left<--2
xmin <- (range(x$avg_log2FC)[1]-(range(x$avg_log2FC)[1]+5))
xmax <- (range(x$avg_log2FC)[1]+(5-range(x$avg_log2FC)[1]))
ymin <- 0
ymax <- range(-log(x$p_val_adj))[2]+1
ymax
ymax <- 150



# Construct the plot object
p1 <- ggplot(data=x, aes(x=avg_log2FC, y=-log10(p_val_adj), label = gene)) +
  geom_point(alpha = 0.6, size=size, colour=x[,7]) +
  scale_color_manual(values = c("lightgrey", "navy", "red")) + 
  
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value-adj")), title="P3_State-3_vs_State-2") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(x_left, -logFCcut, 0, logFCcut,x_right), # tick positions
    labels = c(x_left, -logFCcut, 0, logFCcut,x_right),
    limits = c(x_left, x_right) # x-axis range; symmetric looks better
  ) +
  # Draw the threshold lines
  geom_vline(xintercept = logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_vline(xintercept = -logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_hline(yintercept = -log10(pvalCut), color="grey40", linetype="longdash", size=0.5) +
  
  guides(colour = guide_legend(override.aes = list(shape=16)))+
  
  theme_bw(base_size = 12) +
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black", size=8),
        plot.title = element_text(hjust = 0.8),
        axis.text.x = element_text(face="bold", color="black", size=15),
        axis.text.y = element_text(face="bold",  color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=16),
        axis.title.y = element_text(face="bold",color="black", size=16))
p1

# Label the genes with |logFC| > 1 by type
p2<-p1 + geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val), 
                             label = ifelse((avg_log2FC > 1 & p_val < pvalCut)|(avg_log2FC < -1 & p_val < pvalCut), rownames(x),"")),
                         colour="black", size = 5, box.padding = unit(0.35, "lines"), 
                         point.padding = unit(0.3, "lines"))
p2

ggsave(p2,file="P3_3_vs_2_volcano.pdf", width=6,height=6)
