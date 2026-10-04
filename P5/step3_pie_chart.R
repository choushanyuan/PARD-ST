setwd("/path/to/pseudotime_analysis")
load("Cellratio.RData")
my_color<-c("GS4"="#fcc5c0","GS5-1"="#7a0177","GS5-2"="#ae017e","GS5-3"="#dd3497","GS5-4"="#f768a1","GS5-5"="#fa9fb5")
data=subset(Cellratio,State=='1')
rownames(data)<-data$celltype
data<-data[c("GS4","GS5-1","GS5-2","GS5-3","GS5-4","GS5-5"),]
pie(data$Freq, labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='4')
rownames(data)<-data$celltype
data<-data[c("GS4","GS5-1","GS5-2","GS5-3","GS5-4","GS5-5"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='5')
rownames(data)<-data$celltype
data<-data[c("GS4","GS5-1","GS5-2","GS5-3","GS5-4","GS5-5"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='7')
rownames(data)<-data$celltype
data<-data[c("GS4","GS5-1","GS5-2","GS5-3","GS5-4","GS5-5"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)
