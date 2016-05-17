setwd("..//Master-project/")
load("AllData20160512.RData")
library(reshape)
library(ggplot2)
library(rgl)
library(RColorBrewer)

### Analysis of the data which has removed samples ##########################

#' Formate the output from parser_of_eds to a matrix with samples as columns
#' and miRNAs as rows.
#' 
#' @param table A table from parser_of_eds.R.
#' @param ntargets The number of unique miRNAs for each sample.
#' @param Ctval The column number which contains the Ct values.
#' @param SampleID The column which contains the Sample IDs.
#' @param miRID The column which contains the miRNAs Ids. 
#' @return matrix A matrix with Samples column wise and miRNAs row wise.

edsToMatrix<-function(table, ntargets, Ctval, SampleID, miRID){
  samp<-nrow(table)/ntargets #number of samples in data
  data<-matrix(0,nrow = ntargets, ncol = samp) #creates the new table, rows=genes and cols=samples
  rownames(data)<-as.character(table[1:ntargets,miRID]) #set the rownames to the genenames
  colnames(data)<-colnames(data, do.NULL = F) #creates colnames
  k=ntargets #counter for the end of the sample
  m=1 #counter for the begining of the sample
  for (i in seq(from = 1,to = samp)){ #sets the ddcq values for each sample
    colnames(data)[i]<-table[(k*i),SampleID] #names the column according to sample
    data[,i]<-table[m:(i*k),Ctval] #Adds data according to sample
    m=(i*k+1) #increases beginging of sample
    
  }
  data
}


#Format the data on the same form as microarray
dataAll<-edsToMatrix(all.sum,758,4,1,2)




#' Normalizes a vector by subracting the mean of the vector from each obejct.
#' 
#' @param vec A vector to be normalized
#' @return vec A normalized vector . 
norm<-function(vec){
  vec[vec>38]<-NA #Any value above 38 not included
  vec[vec<10]<-NA #Any value below 10 is not included
  vec-mean(na.omit(vec))
}


#Applys the normalizing function for each column/sample
normDat<-apply(dataAll,2,norm)
normDat<-normDat[-1,] #Remove the first row as it is a endogenous control
normDat<-normDat[-755:-757,] #Remove the the two last rows as it is endogenous control


# Boxplot to check the distribution of the raw data and normalized data
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] #Importes colors
colDens<-colorRampPalette(brewer.pal(11, "Spectral"))(379) #Colour the data in classic densisty colours #takes the 74 most unique colours

boxplot(dataAll, main="Raw data",na.rm=T, xaxt='n',xlab="Samples", ylab="Ct", col=colDens)
boxplot(normDat, main="Normlized data", na.rm=T, xaxt="n", xlab="Samples", ylab="Ct",  col=colDens)

#Plot the distribution of the non normalized and normalized data, one line for each sample
naindexNorm <- apply(normDat, 2, function(x) !all(is.na(x))) #index of all the data which is not NA
densNorm <- apply(normDat[, naindexNorm], 2, density, na.rm = TRUE) #Calculates the density column wise
xNorm <- do.call(cbind, lapply(densNorm, function(d) d$x)) #extracts the xvalues 
yNorm <- do.call(cbind, lapply(densNorm, function(d) d$y)) #extracts the y values
matplot(xNorm, yNorm, col=colDens, type="l", xlab="Ct values", ylab="Density", main="Density, normalized data") #plots the data

indexRaw <- apply(dataAll, 2, function(x) !all(is.na(x)))
densRaw <- apply(dataAll[, indexRaw], 2, density, na.rm = TRUE) #Calculates the density column wise
xRaw <- do.call(cbind, lapply(densRaw, function(d) d$x)) #extracts the xvalues
yRaw <- do.call(cbind, lapply(densRaw, function(d) d$y)) #extracts the y values
matplot(xRaw, yRaw, col=colDens, type="l", ylab="Density", xlab="Ct values", main="Density, raw data") #plots the data

# Load in the phenotype data
biogroup<-read.csv2(file = "..//Data//biogroup_OpenArray.csv")
aaa.samp<-biogroup[biogroup[,2]==1,] # Sample id collected from AAA patients
cont.samp<-biogroup[biogroup[,2]==0,] # Samples id collected from controls

#Divides the data into the biogroups
aaa<-normDat[,colnames(normDat)%in%aaa.samp[,1]]
cont<-normDat[,colnames(normDat)%in%cont.samp[,1]]



#' A function to get a boolean vector of which miRNAs which is expressed in more
#' than 50 % of each biogroup.
#' 
#' @param case A matrix with the case Ct values
#' @param cont A matrix with the control Ct values
#' @return tf.filt A boolean vector with TRUE for the miRNAs passing the filtering
filterData<-function(case,control){
  funx <-function(row){
    sum(is.na(row))/length(row) < 0.5
  }
  tf.mat.case <- apply(case,1,funx)
  tf.mat.cont <- apply(control,1,funx)
  tf.filt <- tf.mat.case & tf.mat.cont 
  
}

#Boolean vector with targets to keep
tf.filter<-filterData(aaa, cont)

#The filtered data for each biogroup
aaa.filter<-aaa[tf.filter,]
cont.filter<-cont[tf.filter,]



#Binds the filter data for PCA
all.filter<-as.matrix(cbind(aaa.filter,cont.filter))
groups<-c(rep("aaa",ncol(aaa.filter)),rep("control",ncol(cont.filter)))


#' A function to replace NAs in a vector with the mean
#' @param vec A vector which to replace the NAs
#' @return vec A vector with the NAs replaced with the mean of the vector
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
pca = prcomp(all.filter.mean, center = T)


#create a dataframe with the different loading vector from the pca
df<-as.data.frame(pca$rotation[,1])
df$x<-pca$rotation[,1]
df$y <- pca$rotation[,2]
df$z<-pca$rotation[,3]
df$groups <- groups # Add the group beloing for each sample
str(df)

#Plots the PCA
ggplot(df, aes(x,y)) + geom_point(aes(color=groups), size=3) + scale_color_manual(values = c('#e34a33','#2c7fb8'))+
  labs( x="PC 1", y="PC 2", colour = "Biogroup", size=3)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=12),
        legend.text=element_text(size=12))

#plot3d(df$x, df$y, df$z,col=c(rep("red", ncol(aaa)), rep("purple", ncol(cont))), cex=20)
#rgl.snapshot("plots///pca3D_All.png", fmt = "png") #Takes a snapshot of the current rgl view



#' A function to perform a unpaired t-test on a vector
#' with known group belogning, to be used with apply(matrix, 2, case,na.rm=T).
#' @param vec A vector with all the data orderd according to group
#' @param case A boolean vector with which columns belong to one of the biogroups 
#' @return pval Returns a vector with pvalues from t-test
matTtest<-function(vec, case){
  pval<-t.test(vec[case], vec[!case], na.action=na.omit)$p.val
  fc<-(-mean(vec[case], na.rm=T))-(-mean(vec[!case], na.rm=T))
  return(pval)
}

#' A function to perform a unpaired t-test on a vector
#' with known group belogning, to be used with apply(matrix, 2, case,na.rm=T).
#' @param matrix A matrix with all the data orderd according to group
#' @param case A boolean vector with which columns belong to one of the biogroups 
#' @return fc Returns a vector with the log2 foldchange with respect to the case group
matFC<-function(vec, case){
  fc<-(-mean(vec[case], na.rm=T))-(-mean(vec[!case], na.rm=T))
  return(fc)
}

#Pvalues and fold change
pval<-apply(all.filter,1,matTtest, case=c(rep(T, ncol(aaa.filter)),rep(F,ncol(cont.filter))))
fc<-apply(all.filter,1,matFC, case=c(rep(T, ncol(aaa.filter)),rep(F,ncol(cont.filter))))


#Adjust the pvalue with false discovery rate
pval.adj<-p.adjust(pval,method = "BH")
k<-which(pval.adj < 0.05) #Gives the targets which are significant

View(cbind(pval[k],fc[k]))



#Dataframe for plotting a volcano plot, store pvalue and foldchange
volcplot.df<-as.data.frame(pval.adj)
volcplot.df$fc<-as.vector(fc)

#Volcano plot, the significant targets higlighted in red
ggplot(volcplot.df, aes(fc,-log10(pval.adj)))+geom_point(aes(color=pval.adj<0.05),size=3) + 
  scale_color_manual(values = c("black", "red"))+ 
  annotate("text", x = 0.7, y =4, label = "Fold change > 0 Uppregulated in AAA ")+
  xlab("Log2 Fold Change") +
  ylab("log10 adjusted Pvalues")+
  ggtitle("Data with removed samples")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        title=element_text(size=12),
        legend.text=element_text(size=12))
 

View(cbind(pval.adj,fc.l)




heatData<-all.filter[k,]















