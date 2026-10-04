setwd("/path/to/pseudotime_analysis")
load("Cellratio.RData")
my_color<-c("GS4-1"="#ae017e","GS4-2"="#f768a1","GS3"="#fbb4b9")
data=subset(Cellratio,State=='1')
rownames(data)<-data$celltype
data<-data[c("GS3","GS4-1","GS4-2"),]
pie(data$Freq, labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='3')
rownames(data)<-data$celltype
data<-data[c("GS3","GS4-1","GS4-2"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='4')
rownames(data)<-data$celltype
data<-data[c("GS3","GS4-1","GS4-2"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)
