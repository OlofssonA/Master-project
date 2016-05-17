setwd("..//Master-project/")
load("AllData20160429.RData")
summary(all.data.sorted1)

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


batches<-read.csv2("..//Data//batch.csv")
removed<-cbind(read.delim("..//Data/removededs2.txt",header = F))


data<-dataform(allData.noDupl,818,4,5,3)

a<-data[,colnames(data)%in%batches[,1]]
b<-data[,colnames(data)%in%batches[,2]]
c<-data[,colnames(data)%in%batches[,3]]
d<-data[,colnames(data)%in%batches[,4]]
e<-data[,colnames(data)%in%batches[,5]]
f<-data[,colnames(data)%in%batches[,6]]
g<-data[,colnames(data)%in%batches[,7]]
h<-data[,colnames(data)%in%batches[,8]]
i<-data[,colnames(data)%in%batches[,9]]
j<-data[,colnames(data)%in%batches[,10]]
k<-data[,colnames(data)%in%batches[,11]]
l<-data[,colnames(data)%in%batches[,12]]
m<-data[,colnames(data)%in%batches[,13]]
n<-data[,colnames(data)%in%batches[,14]]
o<-data[,colnames(data)%in%batches[,15]]
p<-data[,colnames(data)%in%batches[,16]]
q<-data[,colnames(data)%in%batches[,17]]
r<-data[,colnames(data)%in%batches[,18]]
s<-data[,colnames(data)%in%batches[,19]]
t<-data[,colnames(data)%in%batches[,20]]
u<-data[,colnames(data)%in%batches[,21]]
v<-data[,colnames(data)%in%batches[,22]]
w<-data[,colnames(data)%in%batches[,23]]
x<-data[,colnames(data)%in%batches[,24]]

remove<-data[,colnames(data)%in%removed$V1]
used<-data[,!colnames(data)%in%removed$V1]

groups<-c(rep("Removed",ncol(remove)), rep("Used",ncol(used)))
pcaDat<-cbind(remove,used)

groups<-c(rep("a",ncol(a)),rep("b",ncol(b)),rep("c",ncol(c)),rep("d",ncol(d)),rep("e",ncol(e)),rep("f",ncol(f)),
             rep("g",ncol(g)),rep("h",ncol(h)), rep("i",ncol(i)),rep("j",ncol(j)),rep("k",ncol(k)),rep("l",ncol(l)),
             rep("m",ncol(m)),rep("n",ncol(n)), rep("o",ncol(o)), rep("p",ncol(p)), rep("q",ncol(q)) ,rep("r",ncol(r)), rep("s",ncol(s)),
             rep("t",ncol(t)), rep("u",ncol(u)),rep("v",ncol(v)),rep("w",ncol(w)),rep("x",ncol(x)))

pcaDat<-cbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x)

pca1.dcq = prcomp(pcaDat)


#create a dataframe with the different loading vector from the pca
df<-as.data.frame(pca1.dcq$rotation[,1])
df$x<-pca1.dcq$rotation[,1]
df$y <- pca1.dcq$rotation[,2]
df$z<-pca1.dcq$rotation[,3]
df$Batches <- groups # Add the group beloing for each sample

library(RColorBrewer)
nt <- 24
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(df, aes(x,y)) + geom_point(aes(color=Batches), size=3)+ scale_color_manual(values = col_vector)+
  ggtitle("Data coloured after the batches")


pie(rep(1,nt), col=sample(col_vector, nt))


hcd<-as.dendrogram(hclust(dist(t(data))))
plot(hcd)


op = par(mfrow = c(2, 1))
plot(cut(hcd, h = 400)$upper, main = "Upper tree of cut at h=75")
plot(cut(hcd, h = 400)$lower[[22]], main = "Second branch of lower tree with cut at h=75")

par(cex=0.5,font=1)
plot(hcd, main="Dendrogram of Ward's Method")

