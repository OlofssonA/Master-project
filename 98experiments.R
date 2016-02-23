#source("https://bioconductor.org/biocLite.R")
#library(HTqPCR)
library(ggplot2)
library(calibrate)
library(reshape)
# loading the package
library(dendextend)
library(rgl)

setwd(".")

data<-read.csv2("..//export//96samp//dcq_96samp.csv",skip=12, na.string="-", stringsAsFactors = FALSE)
data.glob<-read.csv2("..//export//96samp//glob_96samp.csv",skip=12, na.string="-", stringsAsFactors = FALSE)


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

#source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

sum(data$Sample.Name%in%"1")

#Get the data on the form which i liek
expr.raw<-dataform(data,758,3,1,2)
expr.dcq<-dataform(data,758,6,1,2)
expr.dcqGlob<-dataform(data.glob,758,6,1,2)
expr.globRaw<-dataform(data.glob,758,3,1,2)

#To plot the data
#Creates a datafram
d <- data.frame(x = 1:208, y = expr.raw[3,], sampnames = colnames(expr.raw))
d[is.na(d)]<-0
#Plots the dataframe with annotation to make it easier to read
ggplot(d, aes(x,y)) + geom_point(color=ifelse(d$y <1,"red", "black"),size = 3) + 
  geom_text(aes(label=ifelse(y<1,colnames(expr.raw), ""),hjust=-0.1,just=-2)) + xlab("Sample")+ylab("Mean Cq Endogenous Control") +
  theme(axis.title.x = element_text(face="bold", colour="black", size=15),
        axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y  = element_text(size=10,face="bold"))+scale_y_continuous(breaks=c(0, 10, 20,30),
                                                                             labels=c("N/A", "10", "20","30"))

#Remove all the columns which are NA for all mirnas
expr.dcq<-expr.dcq[, !apply(is.na(expr.dcq), 2, all)]
expr.dcqGlob<-expr.dcqGlob[, !apply(is.na(expr.dcqGlob), 2, all)]

#Boxplot of rawdata
p1<-ggplot(na.omit(melt(expr.raw)), aes(as.factor(X2), value)) + geom_boxplot(color="black")+ylab("Raw mean Cq value")+xlab("Sample")+ggtitle("Raw data")
p2<-ggplot(na.omit(melt(expr.dcq)), aes(as.factor(X2), value)) + geom_boxplot(color="black")+ylab("Delta-Cq value (Normalized)")+xlab("Sample")+ ggtitle("Endogenous control")
p3<-ggplot(na.omit(melt(expr.dcqGlob)), aes(as.factor(X2), value)) + geom_boxplot(color="black")+ylab("Delta-Cq value (Normalized)")+xlab("Sample")+ggtitle("Global normalization")
multiplot(p1, p2, p3, cols=1)

#To get biogroups
sampinfo<-read.csv2("..//export//96samp//glob_96well.csv",skip=12, na.string="-", stringsAsFactors = FALSE)
sample.group<- with(sampinfo, match(unique(Sample.Name), Sample.Name))
sampgroup<-sampinfo[sample.group,4:5]


#seperate the sample names according to biogroups
Contsamp<-sampgroup[sampgroup$Biological.Group.Name%in%"Control",]
Diseasamp<-sampgroup[sampgroup$Biological.Group.Name%in%"Disease",]

#Seperate the samples according to biogorup
aaa.dcq<-expr.dcq[,colnames(expr.dcq)%in%Diseasamp$Sample.Name]
cont.dcq<-expr.dcq[,colnames(expr.dcq)%in%Contsamp$Sample.Name]
aaa.glob<-expr.dcqGlob[,colnames(expr.dcqGlob)%in%Diseasamp$Sample.Name]
cont.glob<-expr.dcqGlob[,colnames(expr.dcqGlob)%in%Contsamp$Sample.Name]

#Filter the data
filterData<-function(case,control){
  tf.filt<-rep(FALSE, nrow(case))
  for (i in 1:nrow(case)){
    tf.filt[i]<-(sum(is.na(case[i,])/ncol(case)) < 0.3 & sum(is.na(control[i,]))/ncol(control)<0.3)
  }
  tf.filt
}

tf.filter.dcq<-filterData(aaa.dcq,cont.dcq)
tf.filter.glob<-filterData(aaa.glob,cont.glob)

#Take only the genes which pass the test
aaa.dcq.filter<-aaa.dcq[tf.filter.dcq,]
cont.dcq.filter<-cont.dcq[tf.filter.dcq,]
aaa.glob.filter<-aaa.glob[tf.filter.glob,]
cont.glob.filter<-cont.glob[tf.filter.glob,]

all.filter.dcq<-as.matrix(cbind(aaa.dcq.filter,cont.dcq.filter))
all.filter.glob<-as.matrix(cbind(aaa.glob.filter,cont.glob.filter))
groups.dcq<-c(rep("aaa",ncol(aaa.dcq.filter)),rep("cont",ncol(cont.dcq.filter)))
groups.glob<-c(rep("aaa",ncol(aaa.glob.filter)),rep("cont",ncol(cont.glob.filter)))
#all.filter.dcq[is.na(all.filter.dcq)]<-0

#Function to replace the na with the mean
f1 <- function(vec) {
  m <- mean(vec, na.rm = TRUE)
  vec[is.na(vec)] <- m
  return(vec)
}

## Another way of dealing with NAs, by setting the missing values to the mean
## and then center the pca it is equal to setting the NAs to zero
#all.filter.mean <- apply(all.filter,1,f1)
#all.filter.mean<-t(all.filter.mean)

#If using f1 function set center=T in prcomp
pca1.dcq = prcomp(na.omit(all.filter.dcq))
#create a dataframe with the different loading vector from the pca
df<-as.data.frame(pca1.dcq$rotation[,1])
df$x<-pca1.dcq$rotation[,1]
df$y <- pca1.dcq$rotation[,2]
df$z<-pca1.dcq$rotation[,3]
df$groups1 <- groups.dcq # Add the group beloing for each sample
str(df) # so you can see the structure and names

pca1.glob = prcomp(na.omit(all.filter.glob))
#create a dataframe with the different loading vector from the pca
df2<-as.data.frame(pca1.glob$rotation[,1])
df2$x<-pca1.glob$rotation[,1]
df2$y <- pca1.glob$rotation[,2]
df2$z<-pca1.glob$rotation[,3]
df2$groups2 <- groups.glob # Add the group beloing for each sample
str(df) # so you can see the structure and names

#Plots the two biggest loading vectors
ggplot(df, aes(x,y)) + geom_point(aes(color=groups1), size=2) + scale_color_manual(values = c("orange", "purple"))
ggplot(df2, aes(x,y)) + geom_point(aes(color=groups2), size=2) + scale_color_manual(values = c("orange", "purple"))

#Plots the 3 biggest loading vectors
plot3d(df2$x, df2$y, df2$z,col=c(rep("red", 124), rep("purple", 84)), cex=20)
rgl.snapshot("plots///pca3D_glob.png", fmt="png") #Takes a snapshot of the current rgl view

#Plot the distribution of the data
plot(density(expr.raw,na.rm=T), main="Raw EC")
plot(density(expr.globRaw,na.rm=T,bw = 0.8229), main="Raw Global control")
plot(density(expr.dcq,na.rm=T), main="Endogenous control normalization")
plot(density(expr.dcqGlob,na.rm=T), main="Global normalization")

#Ttest
pval.glob<-rep(NA,nrow(aaa.glob.filter))
pval.dcq<-rep(NA,nrow(aaa.dcq.filter))

for( i in 1:nrow(aaa.glob.filter)){
  pval.glob[i]<-t.test(aaa.glob.filter[i,],cont.glob.filter[i,],na.action=na.omit)$p.val
}

for( i in 1:nrow(aaa.dcq.filter)){
  pval.dcq[i]<-t.test(aaa.dcq.filter[i,],cont.dcq.filter[i,], na.action = na.omit)$p.val
}

pval.adj.dcq<-p.adjust(pval.dcq,method = "fdr", n=length(pval.dcq))
pval.adj.glob<-p.adjust(pval.glob,method = "fdr", n=length(pval.glob))

names(pval.adj.dcq)<-rownames(aaa.dcq.filter)
names(pval.adj.glob)<-rownames(aaa.glob.filter)

which(pval.adj.dcq<0.05)
which(pval.adj.glob<0.05)

