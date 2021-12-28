# Extract chla data to points
# Dallas Jordan
# sometime in late oct 2021

# Setup -------------------------------------------------------------------

library(gamm4)
library(mgcv)
library(arm)
library(raster)
library(rgdal)
library(lubridate)

# data cleaning: 

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/LAALdata_midway_withTrackID_dtime.Rdata")
LAALmid <- LAAL
LAALmid$id <- paste0("lm",LAALmid$track)
LAALmid$dtime <- ymd_hms(LAALmid$dtime)
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/LAALdata_tern_withTrackID_dtime.Rdata")
LAALtern <- LAAL
LAALtern$id <- paste0("lt",LAALtern$track)
LAALtern <- transform(LAALtern, dtime=paste(day_gmt, time_gmt, sep=" "))
LAALtern <- LAALtern[,c("id","dtime","x","y","track")]
LAALtern$dtime <- parse_date_time(LAALtern$dtime, c("Ymd HM", "mdy HM"))

LAAL <- rbind(LAALmid, LAALtern)

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/BFALdata_midway_withTrackID_dtime.Rdata")
BFALmid <- BFAL
BFALmid$id <- paste0("bm",BFALmid$track)
BFALmid$dtime <- ymd_hms(BFALmid$dtime)
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/BFALdata_tern_withTrackID_dtime.Rdata")
BFALtern <- BFAL
BFALtern$id <- paste0("bt",BFALtern$track)
BFALtern <- transform(BFALtern, dtime=paste(day_gmt, time_gmt, sep=" "))
BFALtern <- BFALtern[,c("id","dtime","x","y","track")]
BFALtern$dtime <- parse_date_time(BFALtern$dtime, c("Ymd HM", "mdy HM"))

BFAL <- rbind(BFALmid, BFALtern)

# comparisons: allLAAL v allBFAL
#              ternLAAL v midwayLAAL
#              ternBFAL v midwayBFAL
#              ternLAAL v ternBFAL
#              midwayLAAL v midwayBFAL
all_data <- rbind(LAAL,BFAL)
all_data$months <- month(all_data$dtime)
counts_movement <- all_data %>% count(months,id)
save(all_data,file="all_tracks_master.Rdata")

# =======================================================================================

# START HERE: 
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

mid_data <- rbind(lm,bm)
tern_data <-rbind(lt,bt)

# fix some dates for later
mid_data$dtime <- as.Date(mid_data$dtime)
tern_data$dtime <- as.Date(tern_data$dtime)
lm$dtime <- as.Date(lm$dtime)
bm$dtime <- as.Date(bm$dtime)
lt$dtime <- as.Date(lt$dtime)
bt$dtime <- as.Date(bt$dtime)

convlon <- function(lon){
  ifelse(lon > 180, lon - 360, lon)
}

lm$x <- convlon(lm$x)
lt$x <- convlon(lt$x)
bm$x <- convlon(bm$x)
bt$x <- convlon(bt$x)

#### Do Midway first, then change to Tern
island <- "tern"
species <- "BFAL"
track.data <- bt

# ======================================================================================

## cycle through each year/ month / day of study period
## select track points that occur on each year/ month/ day and create a buffer around that point(s) 
## then select raster from that date and extract raster values to buffer (unlist values first)

for (yrs in 2007:2012){
  
  for (mths in 1:12) {
    ## this sets the number of days per month 
    if (mths == 1 | mths==3 | mths == 5 | mths == 7 | mths == 8 | mths == 10 | mths == 12) {numdays<-31}
    if (mths == 2) {numdays <- 28} 
    if (mths == 4 | mths==6 | mths == 9 | mths == 11) {numdays<-30}
    datesplit <- "20" 
    # yr<-strsplit(as.character(yrs),split=datesplit)[[1]][2]
    
    ##this puts a 0 before days and months <10 so that each month and day has two numbers 
    if (mths<10) {mths<-paste(0,as.character(mths),sep="")}
    
    chla_data <- data.frame(id=as.character(),
                            dtime=as.Date(character()), 
                            x=numeric(),
                            y=numeric(),
                            track=as.character(),
                            yrs=numeric(),
                            mths=numeric(),
                            dys=numeric(),
                            chla=numeric(),
                            stringsAsFactors=FALSE) 
    names(chla_data)<-names(track.data)
    n_ <- ncol(chla_data)
    names(chla_data)[n_]<-paste("chla")
    names(chla_data)[n_ - 1]<-paste("dys")
    names(chla_data)[n_ -2]<-paste("mths")
    names(chla_data)[n_ -3]<-paste("yrs")
    
    # don't need to loop through days since data is monthly-aggregated
    # for (dys in 1:numdays) {
    #   if (dys<10) {dys<-paste(0,as.character(dys),sep="")}

    dys <- "01"
    print(dys)
    print(mths)
    print(yrs)

    date.temp <- mdy(paste(mths,dys,yrs,sep=" "))
    points <- track.data
    points.temp <- points[which(floor_date(points$dtime,unit="month")==floor_date(date.temp,unit="month")),]
    points <- points.temp[,c(3,4)]
      
      if (nrow(points)>0){
        
        ## In your application, you'll want to make another loop so you're only creating a raster and extracting if there are points on that day (i.e., if nrow(temp)>0)
        
        ### This identifies a workspace for each year of data since rasters were organized by year in different folders
        if(yrs==2007){setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla/monthly/2007/")}
        if(yrs==2008){setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla/monthly/2008/")}
        if(yrs==2009){setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla/monthly/2009/")}
        if(yrs==2010){setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla/monthly/2010/")}
        if(yrs==2011){setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla/monthly/2011/")}
        if(yrs==2012){setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla/monthly/2012/")}
        
        filenm<-paste("CHL_",yrs,mths,dys,".img",sep="")
        
        ## this creates a raster from the filename you created above
        CHLAr <- raster(as.character(filenm))
        
        ## in this code I'm just extracting to points rather than a polygon
        ## I do have code somewhere that creates a polygon and extracts to that polygon rather than points-- can do more looking for it if need be
        ## the important thing is that you extract to the polygon and then unlist the points so you can append them to another dataframe
        chla <- extract(CHLAr, points,buffer=300000, fun=mean)
        
        temp2 <- cbind(points.temp,yrs, mths, dys, chla)
        
        ## here I'm appending the values for that track point to this dataframe that I created at the beginning of the loop
        
        chla_data<-rbind(chla_data,temp2)
        
      } 
    setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla_extracted")
    #output.name2 <-  paste("Daily_env_data_4km_",yrs, "to", mths,".txt",sep="")
    output.name2 <-  paste(island,species,"_chla_extracted",yrs,mths,".txt",sep="")
    write.csv(chla_data, as.character((output.name2)))
  }
}



