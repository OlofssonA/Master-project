library(plyr) 

## This script is to be used with the batch script "unzipRename.sh"
## This script takes the analysis_result.txt file which has been unzipped by the batch script
## and then loads all the files and parse the well number, sample name, target miRNA and Raw Ct value (Cq)
## and then saves the file to a R

setwd('..//Data/')

filenames <- list.files(pattern="*.txt")

all.data<-matrix(rep(0,8),nrow = 2 ,ncol = 4)

for (i in 1:length(filenames)){
con<-file(filenames[i])
open(con);
results.list<-list()
current.line<-1
while(length(line<-readLines(con,n=1,warn=F)) > 0){
  results.list[[current.line]]<- unlist(strsplit(readLines(con), "\t"))
current.line<-current.line+1
}

close(con)


first<-results.list[[1]][16:30]
second<-results.list[[1]][(16+103):(30+103)]
third<-results.list[[1]][(16+206):(30+206)]


b <-seq(1L, length(results.list[[1]]), 103)



sample<-results.list[[1]][16+b]
well<-results.list[[1]][15+b]
target<-results.list[[1]][17+b]
ct<-results.list[[1]][19+b]

dataMerg<-na.omit(cbind(well,sample,target,ct))
all.data<-rbind(all.data,dataMerg)

}

all.data$Sample<-sub(" 0*", " ", all.data$Sample)
colnames(all.data)<-c("well", "Sample", "Target","Ct")
all.data<-all.data[-1:-2,]
all.data<-as.data.frame(all.data, stringsAsFactors=FALSE)

all.data<- data.frame(lapply(all.data, as.character), stringsAsFactors=FALSE)

all.data$Ct<-as.numeric(all.data$Ct)
all.data$Ct[all.data$Ct==100]<-NA

all.sum<-ddply(all.data, c("Sample","Target"),summarize, Ct= 
           mean(na.omit(Ct))) 
commonTargets<-c("RNU44_001094","RNU48_001006", "U6 rRNA_001973","ath-miR159a_000338")
tf.commonTarget<-all.data$Target%in%commonTargets
mean.commontargets<-mean(all.data$Ct[tf.commonTarget], na.rm = T)

tf.sumCommontargets<-all.sum$Target%in%commonTargets
mean.sumCommontarget<-mean(all.sum$Ct[tf.sumCommontargets], na.rm = T)



save(all.data,all.sum, file = "..//Master-project//AllDataQpcr.RData")
