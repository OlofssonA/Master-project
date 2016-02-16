#source("https://bioconductor.org/biocLite.R")
#library(HTqPCR)
library(ggplot2)
library(calibrate)
library(reshape)

setwd(".")

l<-read.csv2("..//export//SampleResults1.csv",skip=12, na.string="-", stringsAsFactors = FALSE)


#'Takes the output from thermofisher cloud and formats the table by taking out samples as columnnames and 
#'targets as rownames with the ddcq as row values for each sample
#'table is the input table and row is the number of targets(number of genes/rows)
#'@param table A table from thermo fisher cloud export
#'@param targets The number of targets each sample has
#'@param dat The column number which has the values for each sample
#'@param sample The column number which holds the sample id
dataform<-function(table,targets, dat, sample, rowname){
  samp<-nrow(table)/targets #number of samples in data
  data<-matrix(0,nrow = targets, ncol = samp) #creates the new table, rows=genes and cols=samples
  rownames(data)<-as.character(l[1:targets,rowname]) #set the rownames to the genenames
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

#Check the number of targets,
sum(l$Sample.Name%in%"Sample 44")

#Reformates the data with the mean Cq values
expr<-dataform(l,758,3,1,2)
expr.dcq<-dataform(l,758,6,1,2)

#Sets the NAs to 0
expr[is.na(expr)]<-0



#To plot the data
#Creates a datafram
d <- data.frame(x = 1:171, y = expr[3,], names = colnames(expr))
#Plots the dataframe with annotation to make it easier to read
ggplot(d, aes(x,y)) + geom_point(color=ifelse(expr[3,]<1,"red", "black"),size = 3) + 
  geom_text(aes(label=ifelse(y<1,colnames(expr), ""),hjust=-0.1,just=-2)) + xlab("Sample")+ylab("Mean Cq Endogenous Control") +
  theme(axis.title.x = element_text(face="bold", colour="black", size=15),
        axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y  = element_text(size=10,face="bold"))+scale_y_continuous(breaks=c(0, 10, 20,30),
                                                                             labels=c("N/A", "10", "20","30"))

#Boxplot of rawdata
ggplot(na.omit(melt(expr)), aes(as.factor(X2), value)) + geom_boxplot(color="black")+ylab("Raw mean Cq value")+xlab("Sample")
ggplot(na.omit(melt(expr.dcq)), aes(as.factor(X2), value)) + geom_boxplot(color="black")+ylab("Delta-Cq value (Normalized)")+xlab("Sample")


#Loads the info to be able to distinguise between the gorups
sampinfo<-read.csv2("..//export//WellResults1.csv",skip=12, na.string="-", stringsAsFactors = FALSE)

#Takes the rows which contains unique sample names
sample.group<- with(sampinfo, match(unique(Sample.Name), Sample.Name))

#Subsets the sample name and group identifier
sampgroup<-sampinfo[sample.group,4:5]


#seperate the sample names according to biogroups
Contsamp<-sampgroup[sampgroup$Biological.Group.Name%in%"Control",]
Diseasamp<-sampgroup[sampgroup$Biological.Group.Name%in%"Disease",]

#Seperate the samples according to biogorup
aaa<-expr[,colnames(expr)%in%Diseasamp$Sample.Name]
cont<-expr[,colnames(expr)%in%Contsamp$Sample.Name]

#Filter the data
tf.filt<-rep(FALSE, nrow(case))
for (i in 1:nrow(case)){
  tf.filt[i]<-(sum(is.na(case[i,])/ncol(case)) < 0.5 & sum(is.na(cont[i,]))/ncol(cont)<0.5)
}

#Take only the genes which pass the test
case.filter<-expr[tf.filt,]
cont.filter<-expr[tf.filt,]


#T-test
pval<-rep(0,nrow(case.filter))
for( i in 1:nrow(case.filter)){
  pval[i]<-t.test(case.filter[i,],cont.filter[i,])$p.value
}


#Pvalue adjustment
min(p.adjust(pval, method ="fdr", n = length(pval)))



