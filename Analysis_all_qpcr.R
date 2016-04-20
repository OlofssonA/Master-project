setwd("..//Master-project/")
load("AllDataQpcrRemovedData.RData")
library(reshape)
library(ggplot2)
library(rgl)


#'Takes the output from thermofisher cloud and formats the table by taking out samples as columnnames and 
#'targets as rownames with the ddcq as row values for each sample
#'table is the input table and row is the number of targets(number of genes/rows)
dataform<-function(table,targets, dat, sample, rowname){
  samp<-nrow(table)/targets #number of samples in data
  data<-matrix(0,nrow = targets, ncol = samp) #creates the new table, rows=genes and cols=samples
  rownames(data)<-as.character(table[1:targets,rowname]) #set the rownames to the genenames
  colnames(data)<-colnames(data, do.NULL = F) #creates colnames
  k=targets #counter for the end of the sample
  m=1 #counter for the begining of the sample
  for (i in seq(from = 1,to = samp)){ #sets the ddcq values for each sample
    colnames(data)[i]<-table[(k*i),sample] #names the column according to sample
    data[,i]<-table[m:(i*k),dat] #Adds data according to sample
    m=(i*k+1) #increases beginging of sample
    
  }
  data
}


#Format the data on the same form as micro array
dataAll<-dataform(all.sum,758,3,1,2)
#dataAll<-dataform(all.sumAmpFilter,758,3,1,2)



#Function to normalize a vector by subtract it by the mean 
norm<-function(vec){
  vec[vec>38]<-NA #Any value above 35 not included
  meanV<-mean(na.omit(vec)) #Removes NA from the mean calculation
  norm<-rep(0, length(vec))
  for(i in 1:length(vec)){
    vec[i]<-vec[i]-meanV
    
  }
  vec
}


#Applys the normalizing function for each column/sample
normDat<-apply(dataAll,2,norm)


boxplot(dataAll, main="Raw data", xlab="Sample", ylab="Ct", na.rm=T)
boxplot(normDat, main="Normlized data, global normalization", na.rm=T)

#Phenodata
biogroup<-read.csv2(file = "..//Data//biogroup_OpenArray.csv")
aaa.samp<-biogroup[biogroup[,2]==1,] # Sample id collected from AAA patients
cont.samp<-biogroup[biogroup[,2]==0,] # Samples id collected from controls

#Divides the data into the biogroups
aaa<-normDat[,colnames(normDat)%in%aaa.samp[,1]]
cont<-normDat[,colnames(normDat)%in%cont.samp[,1]]

#Plot the distribution of the non normalized and normalized data
dens.normDat<-density(na.omit(melt(normDat)$value))
plot(dens.normDat, main="Distribution of Ct normalized data")
plot(density(na.omit(all.sum$Ct)), main="Distribution of raw data")

#Filter the data, filter out targets which has less then half valid ct value of the samples for each biogrup
filterData<-function(case,control){
  tf.filt<-rep(FALSE, nrow(case))
  for (i in 1:nrow(case)){
    tf.filt[i]<-(sum(is.na(case[i,])/ncol(case)) < 0.5 & sum(is.na(control[i,]))/ncol(control)<0.5) 
  }
  tf.filt
}

#Boolean vector with targets to keep
tf.filter<-filterData(aaa, cont)

#The filtered data for each biogroup
aaa.filter<-aaa[tf.filter,]
cont.filter<-cont[tf.filter,]

#Binds the filter data for PCA
all.filter<-as.matrix(cbind(aaa.filter,cont.filter))
groups<-c(rep("aaa",ncol(aaa.filter)),rep("control",ncol(cont.filter)))

#Function to replace the na with the mean, used with pca
f1 <- function(vec) {
  m <- mean(vec, na.rm = TRUE)
  vec[is.na(vec)] <- m
  return(vec)
}

## Another way of dealing with NAs, by setting the missing values to the mean
## and then center the pca it is equal to setting the NAs to zero
all.filter.mean <- apply(all.filter,1,f1)
all.filter.mean<-t(all.filter.mean)

#PCA, center=T when settting NAs to the mean of each target
pca1.dcq = prcomp(all.filter.mean, center = T)


#create a dataframe with the different loading vector from the pca
df<-as.data.frame(pca1.dcq$rotation[,1])
df$x<-pca1.dcq$rotation[,1]
df$y <- pca1.dcq$rotation[,2]
df$z<-pca1.dcq$rotation[,3]
df$groups1 <- groups # Add the group beloing for each sample
str(df)

ggplot(df, aes(x,y)) + geom_point(aes(color=groups1), size=2) + scale_color_manual(values = c("orange", "purple"))
plot3d(df$x, df$y, df$z,col=c(rep("red", ncol(aaa)), rep("purple", ncol(cont))), cex=20)
#rgl.snapshot("plots///pca3D_All.png", fmt = "png") #Takes a snapshot of the current rgl view


#Empty vector to store pvalues and foldchange
pval<-rep(1, nrow(aaa.filter))
fc.l<-rep(1, nrow(aaa.filter))
fc.l2<-rep(1, nrow(aaa.filter))
#Calculate pvalue from ttest and foldchange
for(i in 1:nrow(aaa.filter)){
  pval[i]<-t.test(aaa.filter[i,], cont.filter[i, ], na.action = na.omit)$p.val #Store only pvalue
  fc.l2[i]<-mean(aaa.filter[i,], na.rm=T)-(mean(cont.filter[i,],na.rm=T)) #Foldchange
  fc.l[i]<-2^(-mean(aaa.filter[i,], na.rm=T)-(-mean(cont.filter[i,],na.rm=T))) #Foldchange
}

names(pval)<-rownames(aaa.filter) #Add the targets name to the pvalue
names(fc.l2)<-rownames(aaa.filter) #Add the targets name to foldchange

#Adjust the pvalue with false discovery rate
pval.adj<-p.adjust(pval,method = "fdr")
k<-which(pval.adj < 0.05) #Gives the targets which are significant

#View(cbind(pval.adj[k],fc.l2[k]))



#Dataframe for plotting a volcano plot, store pvalue and foldchange
volcplot.df<-as.data.frame(pval.adj)
volcplot.df$fc.l2<-as.vector(fc.l2)

#Volcano plot, the significant targets higlighted in red
ggplot(volcplot.df, aes(fc.l2,-log10(pval)))+geom_point(aes(color=pval<0.05)) + 
  scale_color_manual(values = c("black", "red"))+geom_vline(xintercept = 0)+
  geom_hline(yintercept = -log10(0.05))+ 
  annotate("text", x = 2, y =5, label = "Fold change > 0 Uppregulated in AAA ")+
  xlab("Log2 Fold Change") +
  ylab("Adjusted log10 Pvalues") 




heatData<-all.filter[k,]
# 
# #If wanting to remove the endengous controls
# commonTargets<-c("RNU44_001094","RNU48_001006", "U6 rRNA_001973","ath-miR159a_000338")
# tf.commonTarget<-rownames(normDat)%in%commonTargets

# 3 ####### Heat Map #######
require(pheatmap)
library(gplots)
library(RColorBrewer)

pheatmap(heatData, cluster_cols = F)












