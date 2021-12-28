# adapting Lesley's code for creating a buffer around points and extracting environmental data

# things I wish I had but don't: 
  # SST for Melinda's data
  # Square buffers

##setwd("D:/Dropbox/Pilot whales/R_folder")

##setwd("~/Dropbox/Pilot whales/R_folder")

setwd("C:/Users/LT/Dropbox/Pilot whales/R_folder")

library(gamm4)
library(mgcv)
library(arm)
library(raster)
library(rgdal)


km4.data <- read.table("4km_grid_new_w_bathy.txt")


DOY <- read.csv("Dates_DOY.csv")
# DOY$Date <- as.POSIXlt(DOY$Date, format="%m/%d/%Y")
DOY$Date <- as.Date(DOY$Date, format="%m/%d/%Y")


SST_SLA_km4<-cbind(km4.data[1,],0,0,0,0,0,0,0)
names(SST_SLA_km4)<-names(km4.data)
n_ <- ncol(SST_SLA_km4)
names(SST_SLA_km4)[n_]<-paste("dfronts")
names(SST_SLA_km4)[n_ - 1]<-paste("SLA")
names(SST_SLA_km4)[n_ -2]<-paste("SST")
names(SST_SLA_km4)[n_ -3]<-paste("dys")
names(SST_SLA_km4)[n_ -4]<-paste("mths")
names(SST_SLA_km4)[n_ -5]<-paste("yrs")
names(SST_SLA_km4)[n_ -6]<-paste("date2")



## cycle through each year/ month / day of study period
## select track points that occur on each year/ month/ day and create a buffer around that point(s) 
## then select raster from that date and extract raster values to buffer (unlist values first)

for (yrs in 2008:2012)  {
  
  for (mths in 1:12) {
    ## this sets the number of days per month 
    if (mths == 1 | mths==3 | mths == 5 | mths == 7 | mths == 8 | mths == 10 | mths == 12) {numdays<-31}
    if (mths == 2) {numdays <- 28} 
    if (mths == 4 | mths==6 | mths == 9 | mths == 11) {numdays<-30}
    datesplit <- "20" 
    
    
    # yr<-strsplit(as.character(yrs),split=datesplit)[[1]][2]
    
    ##this puts a 0 before days and months <10 so that each month and day has two numbers 
    if (mths<10) {mths<-paste(0,as.character(mths),sep="")}
    
    for (dys in 1:numdays) {
      if (dys<10) {dys<-paste(0,as.character(dys),sep="")}
      
      
      
      
      
      ### NOTE FOR DALLAS
      ## ok this step is different than what you'll have to do-- you'll have to pick out points that occur on that particular date
      ## e.g., create a date.temp based on the m/d/year -- could use as.Date(x, format="%m/%d/%Y")
      ##then you can create a temp dataset with the points that occur on that date temp <- data[which(data$date== date.temp), ]
      
      ## in this particular application, I was extracting the raster to the same points every time but it's easy to pick points occurring on a particular date as above
      ## I have code that does picks out points from a date-- I can do some more digging for this if need be
      
      points <- cbind(km4.data$POINT_X,km4.data$POINT_Y)
      
      ## In your application, you'll want to make another loop so you're only creating a raster and extracting if there are points on that day (i.e., if nrow(temp)>0)
      
      
      ### This identifies a workspace for each year of data since rasters were organized by year in different folders
      if(yrs==2010){setwd("D:/JS_SST/2010")}
      if(yrs==2011){setwd("D:/JS_SST/2011")}
      if(yrs==2012){setwd("D:/JS_SST/2012")}
      if(yrs==2013){setwd("D:/JS_SST/2013")}
      if(yrs==2014){setwd("D:/JS_SST/2014G1SST")}
      if(yrs==2015){setwd("D:/JS_SST/2015G1SST")}
      if(yrs==2016){setwd("D:/JS_SST/2016G1SST")}
      
      if(yrs==2010|yrs==2011|yrs==2012|yrs==2013|yrs==2016){
        
        ## this creates a filename for the raster based on the date. I think in this case some years had different name formats so I had to distinguish between years when naming (prob not an issue for you)
        filenm<-paste(yrs,mths,dys,"120000-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST-analysed_sst.img",sep="")	
        
      } else {			
        filenm<-paste(yrs,mths,dys,".tif",sep="")			
      }
      
      ## this creates a raster from the filename you created above
      SSTr <- raster(as.character(filenm))
      
      ## in this code I'm just extracting to points rather than a polygon
      ## I do have code somewhere that creates a polygon and extracts to that polygon rather than points-- can do more looking for it if need be
      ## the important thing is that you extract to the polygon and then unlist the points so you can append them to another dataframe
      SST <- extract(SSTr, points)
      
      
      ## I can't remember why I had to do this-- probably not relevant for you 
      yr<-strsplit(as.character(yrs),split=datesplit)[[1]][2]
      dat <- paste(mths,dys,yr, sep="/")
      date2 <- as.Date(dat, format= "%m/%d/%y")
      DOYr <- DOY[which(DOY$Date==date2),2]
      
      if (DOYr<100&DOYr>9) {DOYr<-paste(0,as.character(DOYr),sep="")
      } else if (DOYr<10) {DOYr<-paste(0,0,as.character(DOYr),sep="")}  
      
      points2 <- points
      points2[,1] <- points2[,1]+360
      
      ## extract SLA
      
      ## this is extracting another type of enviro data-- you could extract all data for that day/ month at the same time in this loop / format
      if(yrs==2010){setwd("D:/SLA/Global/REP/sla/2010")}
      if(yrs==2011){setwd("D:/SLA/Global/REP/sla/2011")}
      if(yrs==2012){setwd("D:/SLA/Global/REP/sla/2012")}
      if(yrs==2013){setwd("D:/SLA/Global/REP/sla/2013")}
      if(yrs==2014){setwd("D:/SLA/Global/REP/sla/2014")}
      if(yrs==2015){setwd("D:/SLA/Global/REP/sla/2015")}
      if(yrs==2016){setwd("D:/SLA/Global/REP/sla/2016")}
      
      
      
      SLAname <- paste("sla_",yrs,DOYr,".img", sep="")
      SLAr <- raster(SLAname)
      SLA <- extract(SLAr, points2)
      
      setwd("D:/JS_SST/D_fronts")
      
      if (yrs == 2016) { dFrontsname <- paste("Dfr_ex_pr_Ex_pr_ex_ex_",yrs,mths,dys,"120000-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST-analysed_sst.img.tif.tif.tif", sep="") }
      ## Dfr_ex_pr_Ex_pr_ex_ex_20160101120000-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST-analysed_sst.img.tif.tif
      if (yrs < 2014) {  dFrontsname <- paste("Dfr_Ex_pr_",yrs,mths,dys,"120000-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST-analysed_sst.img", sep="")}
      if (yrs == 2014| yrs ==2015) {  dFrontsname <- paste("DfrEx_pr_ex_",yrs,mths,dys,".tif", sep="")}
      
      dFrontr <- raster(dFrontsname)
      
      sr <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      
      
      # Project Raster
      pr_dfrontr <- projectRaster(dFrontr, crs = sr)
      dfronts <- extract(pr_dfrontr,points )
      
      temp2 <- cbind(km4.data, date2, yrs, mths, dys, SST, SLA, dfronts)
      
      ## here I'm appending the values for that track point to this dataframe that I created at the beginning of the loop
      
      SST_SLA_km4<-rbind(SST_SLA_km4,temp2)
      
    }
    
    setwd("D:/Daily_input_output_4km")
    ##output.name2 <-  paste("Daily_env_data_4km_",yrs, "to", mths,".txt",sep="")
    output.name2 <-  paste("Daily_env_data_4km_",yrs, mths,".txt",sep="")
    
    write.csv(SST_SLA_km4, as.character((output.name2)))
    
  }
}