
# This quick tool plots SST timeseries from BAS and LOTEK tags and uses cursor to identify the start/end dates of each birds' trip.

# Part B is an addendum that Fixes ALL Colony points to clean up light data during time of pre departure.



#################
# Part A
#################

## BAS -----------------------------------------------------------------------------------------------------------------------------

setwd("/Users/Melinda/Dropbox/GLS_working_directory/BAS_DeCompressed/2011/SST_txt/")

f <- list.files(getwd())

newm2<-matrix(NA,length(f),3)


for (i in 1:length(f)){
  
  rawts<-read.csv(f[i], header=FALSE)
  
  n<-dim(rawts)[1]
  
  plot(c(1:n),rawts$V4,'l')
  coor<-identify(c(1:n),rawts$V4)
  
  start<-strsplit(as.character(rawts[coor[1],2])," ")[[1]][1]
  end<-strsplit(as.character(rawts[coor[2],2])," ")[[1]][1]
  
  newm2[i,1]<-f[i]
  newm2[i,2]<-start
  newm2[i,3]<-end
  
}

colnames(newm2)<-c("file","start_day","end_day")
dp="/Users/Melinda/Dropbox/Chapter_03_Finalize/data/metadata/"
filename = "2011_start_end_BAS.csv"
write.csv(file=paste(dp,filename,sep=""),newm2)


## LTD1400 Data - CSV files -------------------------------------------------------------------------------------------------

setwd("/Users/Melinda/Dropbox/Chapter_03_Finalize/data")

f <- list.files(getwd())

newm<-matrix(NA,length(f),3)

for (i in 1:length(f)){
  
  rawts<-read.csv(f[i])
  
  n<-dim(rawts)[1]
  
  plot(c(1:n),rawts$Ext.Temp.deg.C,'l')
  coor<-identify(c(1:n),rawts$Ext.Temp.deg.C)
  
  start<-strsplit(as.character(rawts[coor[1],1])," ")[[1]][1]
  end<-strsplit(as.character(rawts[coor[2],1])," ")[[1]][1]
  
  newm[i,1]<-f[i]
  newm[i,2]<-start
  newm[i,3]<-end
  
}

colnames(newm)<-c("file","start_day","end_day")
dp="/Users/Melinda/Dropbox/Chapter_03_Finalize/data/metadata/"
filename = "2006_start_end.csv"
write.csv(file=paste(dp,filename,sep=""),newm)





## LAT2500 Data ------------------------------------------------------------------------------------------------------------------------------------

setwd("/Users/Melinda/Dropbox/Chapter_03_Finalize/data/Tern_GLS_allrawfiles/SST_TimeSeries/SST_2009/")

f <- list.files(getwd())

newm2<-matrix(NA,length(f),3)


for (i in 1:length(f)){

  rawts<-read.table(f[i], skip=3, sep=",") 
#    rawts$V3<-rawts$V3*.01
  n<-dim(rawts)[1]
  
  plot(c(1:n),rawts$V3,'l')
  coor<-identify(c(1:n),rawts$V3)
#   rawts<-rawts[1:coor,]
#   plot(c(1:coor),rawts$V3,'l')
#   coor<-identify(c(1:coor),rawts$V3)
  
  start<-strsplit(as.character(rawts[coor[1],2])," ")[[1]][12]
  end<-strsplit(as.character(rawts[coor[2],2])," ")[[1]][11]
  
  newm2[i,1]<-f[i]
  newm2[i,2]<-start
  newm2[i,3]<-end
}

colnames(newm2)<-c("file","start_day","end_day")
dp="/Users/Melinda/Dropbox/Chapter_03_Finalize/data/metadata/"
filename = "2009_start_Lotek.csv"
write.csv(file=paste(dp,filename,sep=""),newm2)




## ------------------------------------------------------------------------------------------------------------------------


#############
# Part B
#############

# Identify Colony Positions for Entire Geolocation track (to clean up trip estimation during chick-rear trips)

# setwd("/Users/Melinda/Dropbox/Chapter_03_Finalize/data/Tern_GLS_allrawfiles/SST_TimeSeries/SST_2009/")
# setwd("/Users/Melinda/Dropbox/Chapter_03_Finalize/data/Tern_GLS_allrawfiles/SST_TimeSeries/Tern_Lotek_040506/")
# setwd("/Users/Melinda/Dropbox/Chapter_03_Finalize/data/Tern_GLS_allrawfiles/SST_TimeSeries/Tern_Lotek_0810/")

setwd("/Users/Melinda/Dropbox/GLS_working_directory/BAS_DeCompressed/2012/SSTtxt/")

f <- list.files(getwd())

W <- list()

for (i in 1:length(f)){
  
  rawts<-read.table(f[i], sep=",",header=TRUE) 
  n<-dim(rawts)[1]
  
  plot(c(1:n),rawts[,7],'l')
  coor<-identify(c(1:n),rawts[,7])
  t<-as.character(rawts[coor,1])
  
  n<-nchar(f[i])
  
  W[["bird"]][i]<-substr(f[i],1,n-7)
  W[["timesAtTern"]][i]<-data.frame(t)
  
}

save(Wsub,file="BAS2012_TernSST_Flag.Rdata")








