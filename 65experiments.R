source("https://bioconductor.org/biocLite.R")
library(HTqPCR)

setwd(".")

l<-read.csv2("..//export//SampleResults1.csv",skip=12, na.string="-", stringsAsFactors = FALSE)


#Takes the output from thermofisher cloud and formats the table by taking out samples as columnnames and 
#targets as rownames with the ddcq as row values for each sample
#table is the input table and row is the number of targets(number of genes/rows)
dataform<-function(table,targets){
  samp<-nrow(table)/targets #number of samples in data
  data<-matrix(0,nrow = row, ncol = samp) #creates the new table, rows=genes and cols=samples
  rownames(data)<-as.character(l[1:targets,2]) #set the rownames to the genenames
  colnames(data)<-colnames(data, do.NULL = F) #creates colnames
  k=targets #counter for the end of the sample
  m=1 #counter for the begining of the sample
  for (i in seq(from = 1,to = samps)){ #sets the ddcq values for each sample
    colnames(data)[i]<-l[(k*i),1] #names the column according to sample
    data[,i]<-l[m:(i*k),17] #Adds data according to sample
    m=(i*k+1) #increases beginging of sample
    
  }
  
}

