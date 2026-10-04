setwd("/path/to/pseudotime_analysis")
load("Cellratio.RData")
my_color<-c("PIN-Low"="#fcc5c0","PIN-High"="#fa9fb5","GS4"="#f768a1","GS4+5"="#c51b8a","GS5"="#7a0177")
data=subset(Cellratio,State=='1')
rownames(data)<-data$celltype
data<-data[c("PIN-Low","PIN-High","GS4","GS4+5","GS5"),]
pie(data$Freq, labels = c("PIN-Low","PIN-High","GS4","GS4+5","GS5"), border="white", col=my_color)

data=subset(Cellratio,State=='2')
rownames(data)<-data$celltype
data<-data[c("PIN-Low","PIN-High","GS4","GS4+5","GS5"),]
pie(data$Freq , labels = c("PIN-Low","PIN-High","GS4","GS4+5","GS5"), border="white", col=my_color)

data=subset(Cellratio,State=='3')
rownames(data)<-data$celltype
data<-data[c("PIN-Low","PIN-High","GS4","GS4+5","GS5"),]
pie(data$Freq , labels = c("PIN-Low","PIN-High","GS4","GS4+5","GS5"), border="white", col=my_color)
