#source("https://bioconductor.org/biocLite.R")
#library(HTqPCR)
library(ggplot2)
library(calibrate)
library(reshape)
# loading the package
library(dendextend)

setwd(".")

l<-read.csv2("..//export//SampleResults1.csv",skip=12, na.string="-", stringsAsFactors = FALSE)
l.glob<-read.csv2("..//export//SampleResultGlobnorm.csv",skip=12, na.string="-", stringsAsFactors = FALSE)


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

#Check the number of targets,
sum(l$Sample.Name%in%"Sample 44")

#Reformates the data with the mean Cq values
expr<-dataform(l,758,3,1,2)
expr.dcq<-dataform(l,758,6,1,2)
expr.dcqGlob<-dataform(l.glob,758,6,1,2)
expr.glob<-dataform(l.glob,758,3,1,2)





#To plot the data
#Creates a datafram
d <- data.frame(x = 1:171, y = expr.glob[3,], names = colnames(expr))
#Plots the dataframe with annotation to make it easier to read
ggplot(d, aes(x,y)) + geom_point(color=ifelse(expr[3,]<1,"red", "black"),size = 3) + 
  geom_text(aes(label=ifelse(y<1,colnames(expr), ""),hjust=-0.1,just=-2)) + xlab("Sample")+ylab("Mean Cq Endogenous Control") +
  theme(axis.title.x = element_text(face="bold", colour="black", size=15),
        axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y  = element_text(size=10,face="bold"))+scale_y_continuous(breaks=c(0, 10, 20,30),
                                                                             labels=c("N/A", "10", "20","30"))

#Boxplot of rawdata
p1<-ggplot(na.omit(melt(expr)), aes(as.factor(X2), value)) + geom_boxplot(color="black")+ylab("Raw mean Cq value")+xlab("Sample")+ggtitle("Raw data")
p2<-ggplot(na.omit(melt(expr.dcq)), aes(as.factor(X2), value)) + geom_boxplot(color="black")+ylab("Delta-Cq value (Normalized)")+xlab("Sample")+ ggtitle("Endogenous control")
p3<-ggplot(na.omit(melt(expr.dcqGlob)), aes(as.factor(X2), value)) + geom_boxplot(color="black")+ylab("Delta-Cq value (Normalized)")+xlab("Sample")+ggtitle("Global normalization")
multiplot(p1, p2, p3, cols=1)

# #Sets the NAs to 0
# expr[is.na(expr)]<-0
# expr.dcq[is.na(expr.dcq)]<-0

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
aaa<-expr.dcq[,colnames(expr.dcq)%in%Diseasamp$Sample.Name]
cont<-expr.dcq[,colnames(expr.dcq)%in%Contsamp$Sample.Name]

#Filter the data
filterData<-function(case,control){
tf.filt<-rep(FALSE, nrow(case))
for (i in 1:nrow(case)){
  tf.filt[i]<-(sum(is.na(case[i,])/ncol(case)) < 0.5 & sum(is.na(control[i,]))/ncol(control)<0.5)
}
tf.filt
}

tf.filter<-filterData(aaa,cont)
#Take only the genes which pass the test
aaa.filter<-aaa[tf.filter,]
cont.filter<-cont[tf.filter,]

all.filter<-as.matrix(cbind(aaa.filter,cont.filter))
groups<-c(rep("aaa",65),rep("cont",106))
all.filter[is.na(all.filter)]<-0

#Function to replace the na with the mean
f1 <- function(vec) {
  m <- mean(vec, na.rm = TRUE)
  vec[is.na(vec)] <- m
  return(vec)
}

# all.filter.mean <- apply(all.filter,1,f1)
# all.filter.mean<-t(all.filter.mean)

#If using f1 function set center=T in prcomp
pca1 = prcomp(na.omit(all.filter.mean))
#create a dataframe with the different loading vector from the pca
df<-as.data.frame(pca1$rotation[,1])
df$x<-pca1$rotation[,1]
df$y <- pca1$rotation[,2]
df$z<-pca1$rotation[,3]
df$groups <- groups # Add the group beloing for each sample
str(df) # so you can see the structure and names

#Plots the two biggest loading vectors
ggplot(df, aes(x,y)) + geom_point(aes(color=groups), size=2) + scale_color_manual(values = c("orange", "purple"))

#Plots the 3 biggest loading vectors
plot3d(df$x, df$y, df$z,col=c(rep("orange", 65), rep("purple", 103)), cex=10)
rgl.snapshot("plots///pca3D_EC.png", fmt="png") #Takes a snapshot of the current rgl view
