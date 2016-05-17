require(plyr) 
require(prodlim)

## This script is to be used with the batch script "unzipRename.sh"
## This script takes the analysis_result.txt file which has been unzipped by the batch script
## and then loads all the files and parse the well number, sample name, target miRNA and Raw Ct value (Cq)
## and then saves the file to a R

#The folder with the txt files
setwd('..//Data/')

# # Names of all the txt files in that folder
# filenamesAmp <- list.files(pattern="*.amp.txt")
# k<-read.delim("NQH27_OA_miRNA_Human.amp.txt", skip=1)
 filenamesCt<-list.files(pattern="*.raw.txt")
# 
# #Get the amplification score
# ampScore<-as.data.frame(cbind(0,0,0,0))
# colnames(ampScore)<-c("Well", "Target", "AmpScore", "Barcode")
# for(i in 1:length(filenamesAmp) ){
#   file<-read.delim(filenamesAmp[i], skip=1)
#   barcode<-sub(" *_OA_miRNA_Human.amp.txt", "", filenamesAmp[i])
#   amp1<-as.data.frame(cbind(file,rep(barcode,nrow(file))))
#   colnames(amp1)<-c("Well", "Target", "AmpScore", "Barcode")
#   ampScore<-rbind(ampScore,amp1)
#   
# }
# ampScore<-ampScore[-1,]


#creates a matrix to be able to store the data
all.data<-matrix(rep(0,10),nrow = 2 ,ncol = 5)

#For each txt file
for (i in 1:length(filenamesCt)){
  con<-file(filenamesCt[i])  #creates a connection
  open(con);              #Opens the connection
  results.list<-list()    #Creates a list to the connection in
  current.line<-1         # Counter for which line to read
  while(length(line<-readLines(con,n=1,warn=F)) > 0){ # Loop as long the line isnt empty
    results.list[[current.line]]<- unlist(strsplit(readLines(con), "\t")) #Split the line on tab and stores it in the list
    current.line<-current.line+1  # Next line
  }
  
  close(con) #Closes the connection
  
  # 
  # first<-results.list[[1]][16:30]
  # second<-results.list[[1]][(16+103):(30+103)]
  # third<-results.list[[1]][(16+206):(30+206)]
  
  
  b <-seq(1L, length(results.list[[1]]), 103) # The number of lines between each sample
  
  
  
  sample<-results.list[[1]][16+b] #All the sample names
  well<-results.list[[1]][15+b]  #All the wells
  target<-results.list[[1]][17+b] #All the targets
  ct<-results.list[[1]][19+b]     # The Ct values
  barcode2<-sub(" *_OA_miRNA_Human.raw.txt", "", filenamesCt[i])
  Barcode<-rep(barcode2, length(ct))
  
  
  dataMerg<-na.omit(cbind(well,sample,target,ct,Barcode)) # Binds all the data column wise
  all.data<-rbind(all.data,dataMerg) #Binds the data with the rest of the data rowwise
  
}


colnames(all.data)<-c("Well", "Sample", "Target","Ct", "Barcode") #Changes the name of the columns

all.data<-all.data[-1:-2,] #Removes the two first empty rows 
all.data<-as.data.frame(all.data, stringsAsFactors=FALSE) # Make it into dataframe
all.data$Ct<-as.numeric(all.data$Ct) #Makes the Ct values numeric
#all.data$Ct[all.data$Ct==100]<-NA #Sets all the Ct=100 to NA


#Changes the well information to numeric
all.data$Well <-as.numeric(all.data$Well)

#Sort the Ct values and ampscore values match
all.dataOrderd<-all.data[with(all.data,order(Barcode,Well)),]



#Sort the data on the sample name and barcode id
all.dataOrderd1<-all.data[with(all.data,order(Barcode,Sample)),]



#Remove all the 0 infront of a sample number, this to deal with replicates
all.dataOrderd1$Sample<-sub(" 0*", " ", all.dataOrderd1$Sample) 
all.dataOrderd1$Barcode<-sub(" ", "", all.dataOrderd1$Barcode) 




removeSamp<-read.csv2("rerunnedSamp.csv")
removeSamp$Barcode<-sub(" ", "", removeSamp$Barcode)



tf.Remove<-(all.dataOrderd1$Sample%in%removeSamp$Sample & all.dataOrderd1$Barcode%in%removeSamp$Barcode)


alldat<-as.data.frame(cbind(all.dataOrderd1$Sample,all.dataOrderd1$Barcode))

h<-row.match(alldat, removeSamp, nomatch = NA)
j<-is.equal(alldat,removeSamp)

m<-all.dataOrderd1[is.na(h),]

summ<-rep(0,nrow(removeSamp))
for(i in 1:nrow(removeSamp)){
  summ[i]<-sum(m$Barcode==removeSamp$Barcode[i] & m$Sample==removeSamp$Sample[i])
  
}





#Takes the mean of Ct values with the same target and sample id 
all.sum<-ddply(m, c("Sample","Target","Barcode"),summarize, Ct= 
                 mean(na.omit(Ct))) 


# all.sumAmpFilter<-ddply(all.data.AmpThre, c("Sample","Target"),summarize, Ct= 
#                           mean(na.omit(Ct))) 


# # Saves the data
#save(allData.noDupl, all.sum, file = "..//Master-project//AllData20160429.RData")
