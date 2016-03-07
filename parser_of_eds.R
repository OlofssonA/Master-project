setwd('..//EXjob/Data/')

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

colnames(all.data)<-c("well", "Sample", "Target","Ct")
