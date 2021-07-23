## This script is used to identify, manually, the start and stop of post-breeding albatross migrations by inspecting SST timeseries. 
## Melinda Conners August 2003, adapted by Dallas Jordan Nov 2020

# When plot appears - the program will wait for you to click twice on the plot to indicate the start and stop indices of the migration. Click once for each of the two locations, and then click 'Finish' once done. 


# Required functions -----------------------------------------------------
clean_date<-function(date_messy, dividr) {
  # date_messy <- date in character format
  # dividr <- how the date units are separated (usually "/", but sometimes "-")
  t<-as.character(date_messy)
  lt<-strsplit(t,dividr)
  for (k in 1:length(lt)) {
    mt<-lt[[k]][1] #month
    lt[[k]][1]<-ifelse(nchar(mt)<2,paste("0",mt,sep=""),mt) # If month is one digit, make into two digit
    dt<-lt[[k]][2] #day
    lt[[k]][2]<-ifelse(nchar(dt)<2,paste("0",dt,sep=""),dt) # If day is one digit, make into two digit
  }
  
  tv<-matrix(NA, nrow=length(lt),ncol=1)
  tv<-as.data.frame(tv)
  for (k in 1:length(lt)) {
    tv[k,1]<-paste(lt[[k]][1],lt[[k]][2],lt[[k]][3], sep="/")
  }
  tv<-as.character(tv)
  return(tv)
}


## LAT2500 Data ------------------------------------------------------------------------------------------------

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation") # Directory to the SST Timeseries

f <- list.files("./data/daylog_work/sst_timeseries",".TXT") # List of .TXT files in directory - these should be your SST .txt logs

newm<-as.data.frame(matrix(NA,length(f),3)) # Create an empty matrix to fill 
colnames(newm)<-c("file","start_day","end_day")

for (i in 1:length(f)){
  
  # Read in file i , skip 3 lines due to headers
  rawts<-read.table(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_timeseries/",f[i],sep=""), skip=3, sep=",") 
  colnames(rawts) <- c("index", "date_time", "SST")
  
  # Plot SST Timeseries (use index (n) and not datetime on X axis, because it will save a lot of time)
  n<-dim(rawts)[1] 
  plot(c(1:n),rawts$SST,'l') 
  #Identify "Start" and "Stop" of PostBreeding Migration by clicking on timeseries. Press 'Finish' when done. This saves coordinates in coor matrix. 
  coor<-identify(c(1:n),rawts$SST) 
  
  # store start and end as datetime - requires a little bit of character-smithing
  startvec<-unlist(strsplit(as.character(rawts[coor[1],2])," "))
  startdate<-startvec[which(nchar(trimws(startvec))!=0)][1]
  
  endvec<-unlist(strsplit(as.character(rawts[coor[2],2])," "))
  enddate<-endvec[which(nchar(trimws(endvec))!=0)][1]
  
  cleanstart<-clean_date(startdate,"/")
  cleanend<-clean_date(enddate,"/")
  
  # Add to dataframe
  newm$file[i]<-f[i]
  newm$start_day[i]<-cleanstart
  newm$end_day[i]<-cleanend
}


dp="/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_timeseries_output/"
filename = "Midway_startstop_Lotek.txt"
write.csv(file=paste(dp,filename,sep=""),newm, row.names=FALSE)
write.table(newm, file = filename, sep = ',', row.names=FALSE)