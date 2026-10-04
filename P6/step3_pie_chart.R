setwd("/path/to/pseudotime_analysis")
load("Cellratio.RData")
my_color<-c("GS3_1"="#fa9fb5","GS3+4_1"="#c51b8a","GS4"="#7a0177","GS3_2"="#fcc5c0","GS3+4_2"="#f768a1")
data=subset(Cellratio,State=='1')
rownames(data)<-data$celltype
data<-data[c("GS3_1","GS3+4_1","GS4","GS3_2","GS3+4_2"),]
pie(data$Freq, labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='2')
rownames(data)<-data$celltype
data<-data[c("GS3_1","GS3+4_1","GS4","GS3_2","GS3+4_2"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='3')
rownames(data)<-data$celltype
data<-data[c("GS3_1","GS3+4_1","GS4","GS3_2","GS3+4_2"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)
