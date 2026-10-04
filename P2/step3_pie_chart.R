setwd("/path/to/pseudotime_analysis")
load("Cellratio.RData")
my_color<-c("subclone_A"="#feb24c","subclone_B"="#f03b20")
data=subset(Cellratio,State=='1')
rownames(data)<-data$celltype
data<-data[c("subclone_A","subclone_B"),]
pie(data$Freq, labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='3')
rownames(data)<-data$celltype
data<-data[c("subclone_A","subclone_B"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)

data=subset(Cellratio,State=='5')
rownames(data)<-data$celltype
data<-data[c("subclone_A","subclone_B"),]
pie(data$Freq , labels = c(""), border="white", col=my_color)
