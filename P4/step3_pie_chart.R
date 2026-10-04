setwd("/path/to/pseudotime_analysis")
load("Cellratio.RData")
my_color<-c("GS3"="#fcc5c0","GS4"="#fa9fb5","GS4+5_2"="#f768a1","GS4+5_1"="#c51b8a","GS5"="#7a0177")
data=subset(Cellratio,State=='1')
rownames(data)<-data$celltype
data<-data[c("GS3","GS4","GS4+5_1","GS4+5_2","GS5"),]
pie(data$Freq, labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='2')
rownames(data)<-data$celltype
data<-data[c("GS3","GS4","GS4+5_1","GS4+5_2","GS5"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='3')
rownames(data)<-data$celltype
data<-data[c("GS3","GS4","GS4+5_1","GS4+5_2","GS5"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)
