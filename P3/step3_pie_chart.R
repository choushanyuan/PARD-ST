setwd("/path/to/pseudotime_analysis")
load("Cellratio.RData")
my_color<-c("PIN"="#fbb4b9","GS3"="#f768a1","GS4"="#c51b8a","GS4+5"="#7a0177")
data=subset(Cellratio,State=='1')
rownames(data)<-data$celltype
data<-data[c("PIN","GS3","GS4","GS4+5"),]
pie(data$Freq, labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='2')
rownames(data)<-data$celltype
data<-data[c("PIN","GS3","GS4","GS4+5"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='3')
rownames(data)<-data$celltype
data<-data[c("PIN","GS3","GS4","GS4+5"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)
