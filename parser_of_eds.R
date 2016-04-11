require(plyr) 

## This script is to be used with the batch script "unzipRename.sh"
## This script takes the analysis_result.txt file which has been unzipped by the batch script
## and then loads all the files and parse the well number, sample name, target miRNA and Raw Ct value (Cq)
## and then saves the file to a R

#The folder with the txt files
setwd('..//Data/')

# Names of all the txt files in that folder
filenames <- list.files(pattern="*.txt")

#creates a matrix to be able to store the data
all.data<-matrix(rep(0,8),nrow = 2 ,ncol = 4)

#For each txt file
for (i in 1:length(filenames)){
  con<-file(filenames[i])  #creates a connection
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
  
  dataMerg<-na.omit(cbind(well,sample,target,ct)) # Binds all the data column wise
  all.data<-rbind(all.data,dataMerg) #Binds the data with the rest of the data rowwise
  
}


colnames(all.data)<-c("well", "Sample", "Target","Ct") #Changes the name of the columns
all.data<-all.data[-1:-2,] #Removes the two first empty rows 
all.data<-as.data.frame(all.data, stringsAsFactors=FALSE) # Make it into dataframe
all.data$Sample<-sub(" 0*", " ", all.data$Sample) #Remove all the 0 infront of a sample number



all.data$Ct<-as.numeric(all.data$Ct) #Makes the Ct values numeric
all.data$Ct[all.data$Ct==100]<-NA #Sets all the Ct=100 to NA

#Takes the mean of Ct values with the same target and sample id 
all.sum<-ddply(all.data, c("Sample","Target"),summarize, Ct= 
                 mean(na.omit(Ct))) 



# # Saves the data
#save(all.data,all.sum, file = "..//Master-project//AllDataQpcr.RData")
