# 
# Then annotation was dowloaded from the hta20sttranscriptcluster.db which contatins all the information
# for affymetrix HTA 2.0 array (transcriptome arry), 
# ACCUM=acession number
# SYMBOL= genesymbol
# DESC= description
# GO= GO annotations
# 
# The normalised data + annotation was saved to a R workspace to save memory,
# code used to normalize and download annotation below:

setwd("../hta_Uppsala/")
source("http://bioconductor.org/biocLite.R")
library(oligo)
library(hta20sttranscriptcluster.db)
library(limma)
library(genefilter)


Annot <- data.frame(ACCNUM=sapply(contents(hta20sttranscriptclusterACCNUM), paste, collapse=", "), 
                    SYMBOL=sapply(contents(hta20sttranscriptclusterSYMBOL), paste, collapse=", "), 
                    DESC=sapply(contents(hta20sttranscriptclusterGENENAME), paste, collapse=", "),
                    GO=sapply(contents(hta20sttranscriptclusterGO), paste, collapse=", "))


list<-list.celfiles() #list the cel files in current workspace
celfiles<-read.celfiles(list) #reads the files into R
eset <- rma(celfiles, normalize=TRUE) #normalizes the celfiles
data<-as.data.frame(exprs(eset)) # subset the normalized expression values
annotation(eset)<-"hta20sttranscriptcluster.db"
sampleNames(eset)<-sub("\\_\\(HTA-2_0).CEL","",sampleNames(eset))
# pm<-pm(celfiles)
# 
# 
# pdf(file = "MAplot.pdf")
# MAplot(eset)
# dev.off()

eset.filterd <- nsFilter(eset, remove.dupEntrez = FALSE,var.cutoff = 0.5)$eset


all<-merge(data,Annot,by.x=0,by.y=0,) #merge by rownames
symbol<-all$SYMBOL  #subset symbol as only needed in this task
names(symbol)<-all$Row.names # add probset id as rownames to id

biogroup<-read.csv2("Annot.csv", header = T,)
rownames(biogroup)<-biogroup$probe.number.in.HTA

data<-exprs(eset)


aaa<-data[,colnames(data)%in%rownames(biogroup[biogroup$group.annotation ==0,])]
cont<-data[,colnames(data)%in%rownames(biogroup[biogroup$group.annotation ==1,])]

mean.aaa<-apply(aaa,1,mean)
mean.cont<-apply(cont,1,mean)

pval<-rep(NA, length(mean.aaa))
fc<-rep(NA,length(mean.aaa))

for(i in 1:length(mean.aaa)){
  pval[i]<-t.test(aaa[i,],cont[i,])$p.val
  fc[i]<-mean.aaa[i]-mean.cont[i]
  
}

groups<-c(rep("aaa",ncol(aaa)),rep("cont",ncol(cont)))
all.data<-c(aaa,cont)
pval.adj<-p.adjust(pval,method = "fdr", n=length(pval))
pca1 = prcomp(cbind(aaa,cont))
df<-as.data.frame(pca1$rotation[,1])
df$x<-pca1$rotation[,1]
df$y <- pca1$rotation[,2]
df$z<-pca1$rotation[,3]
df$groups1 <- groups

ggplot(df, aes(x,y)) + geom_point(aes(color=groups1), size=2) + scale_color_manual(values = c("orange", "purple"))
plot3d(df$x, df$y, df$z,col=c(rep("red", ncol(aaa)), rep("purple", ncol(cont))), cex=20)
