# probGLS implementation for Scott Shaffer's Midway albatross GLS data

# This script is used for 2011 and 2012 tag file downloads, which are in M/D/YYYY format and have split wet dry/temp logs
# This script is "pared down" from my other 01_ source code scripts for efficiency. Those scripts have many useful functions. 
# Dallas Jordan
# last edited: Nov. 2020
# Code citation: 
# [1] Geolocation Manual (Lisovski et al.,
# https://geolocationmanual.vogelwarte.ch/)

# For each tag:
# load in daylog and temperature log
# trim to first/last date as noted in data_inventory
# convert sunrise/sunset to HH:MM, 24-hour clock
# GMT-11 for all sunrise/sunset times for reference; probGLS uses GMT
# check that trn has full range of dates
# prep additional wetdry/immersion data (sen/act)
# run algorithm, plot, extract coords of geographic median

############################################################################################################

library(readr)
library(maptools)
library(GeoLocTools)
setupGeolocation()
library(mapdata)
library(rgeos)
library(raster)
library(chron)
library(suncalc)
library(ncdf4)
library(dplyr)
library(probGLS)
library(adehabitatHR)
library(lubridate)

############################################################################################################

# info per-tag - these entries require manual changes referencing the "Daylog Matching" tab of 
# the Data Inventory spreadsheet. This should be the only section that requires manual changing
# in between tag analyses

Species <- "BFAL"
Year <- "2111"
TrackNumber <- "07" 
HatchYear <- 2011

# first/last date as noted in data_inventory

# "MM/D/YYYY or MM/DD/YYYY", no 0 for e.g. 05, for older daylogs
actual_first <- "1/10/2011"
actual_last <- "12/30/2011"
start <- as.POSIXct("2011-01-09 23:59", tz="Pacific/Midway") # Time when tagged, 11:59pm of previous day if unknown
end   <- as.POSIXct("2011-12-31 00:01", tz="Pacific/Midway") # Time when tag was retrieved, 12:01am of next day if unknown 

# for all Midway tags 

lat.calib <- 28.2
lon.calib <- 182.65
wetdry.resolution <- 32 # sampling rate of wetdry log (immersion) in seconds, e.g. once every 5 min = 300 seconds

############################################################################################################

# CSVs loaded in here. Note: each of these csv files were converted from raw .bin files using Lotek 
# TagTalk software. The converted csv files had to rows above the data, stating "GMT Time Correction = 0" 
# and either "Lg00 (Day Log) Log:" or "Lg20 (Basic Log) Log:" depending on the type of log. These two rows
# have been deleted to make read-in easier. 

wd <- getwd()

# load in daylog
data <- read.csv(paste0("data/daylog_work/daylogs_matched/",Species,"/","dl_",Year,"_",TrackNumber, ".csv"),
                 sep=",", header = T,row.names = NULL, stringsAsFactors = F, skip = 2) # may need to change sep=","

# load in additional info log (wetdry/immersion, called "temperature logs")
# wetdry
tl_data <- read.csv(paste0("data/daylog_work/temperature_logs_matched/",Species,"/",HatchYear,"_preprocessing","/tl_",Year,"_",TrackNumber,"_wetdry.csv"),
                    sep=",", header = T,row.names = NULL, stringsAsFactors = F, skip = 2) # may need to change sep=","
IntTemp_data <- read.csv(paste0("data/daylog_work/temperature_logs_matched/",Species,"/",HatchYear,"_preprocessing","/tl_",Year,"_",TrackNumber,"_inttemp.csv"),
                         sep=",", header = T,row.names = NULL, stringsAsFactors = F, skip = 2) # may need to change sep=","
      IntTemp_data$Date <- strptime(IntTemp_data$Date,format="%m/%d/%Y", tz="GMT")
      dateFilter1 <- strptime(actual_first,format="%m/%d/%Y", tz="GMT")
      dateFilter2 <- strptime(actual_last,format="%m/%d/%Y", tz="GMT")
# trim to first/last date as noted in data_inventory

data <- data[min(which(grepl(actual_first,data$Date))):max(which(grepl(actual_last,data$Date))),]
tl_data <- tl_data[min(which(grepl(actual_first,tl_data$Date))):max(which(grepl(actual_last,tl_data$Date))),]
IntTemp_data <- IntTemp_data %>%
  filter(Date >= dateFilter1 & Date <= dateFilter2)

# visualize daylight period
RiseMinute <- data$Sunrise
RiseMinute <- RiseMinute[RiseMinute<2600] # excludes NA values
SetMinute <- data$Sunset
SetMinute <- SetMinute[SetMinute<2600] # excludes NA values
date_plot <- strptime(as.character(data$Date), format="%m/%d/%Y", tz="GMT") # might need to change the Y to a y if graph looks weird
date_plot <- as.POSIXct(date_plot)
plot(RiseMinute ~ date_plot, type = "n", ylim = range(c(RiseMinute, SetMinute)), ylab = "tag day-time")
polygon(c(date_plot, rev(date_plot)),c(RiseMinute, rev(SetMinute)), col = "yellow")
abline(h = 1440) # 24 hour period, 1440 minutes = 24 hours
title("daylight period (minute of the day in UTC)")

############################################################################################################

# lotek_to_dataframe to convert your calculated daylog into a format usable by probGLS
# necessary formats are as follows:
# date	column "TimeS" in daylog file as as.Date() [DONE]
# sunrise column "Sunrise" and "Sunset" in daylog file as as.character()
# daylog files have this as "UTC minute of the day". Must convert to HH:MM format!
# lon/lat is in numeric [DONE]

data$Date <- strptime(as.character(data$Date), format="%m/%d/%Y", tz="GMT")  # format may change depending on what daylog you are loading

FullRiseTime <- data$Date+minutes(data$Sunrise)
data$RiseTimeUTC <- paste(sprintf("%02d", hour(FullRiseTime)), 
                          sprintf("%02d", minute(FullRiseTime)), 
                          sep=":")
MidwayRiseTime <- FullRiseTime-hours(11)
data$RiseTimeMidway <- paste(sprintf("%02d", hour(MidwayRiseTime )), 
                             sprintf("%02d", minute(MidwayRiseTime )), 
                             sep=":")

FullSetTime <- data$Date+minutes(data$Sunset)
data$SetTimeUTC <- paste(sprintf("%02d", hour(FullSetTime)), 
                         sprintf("%02d", minute(FullSetTime)), 
                         sep=":")
MidwaySetTime <- FullSetTime-hours(11)
data$SetTimeMidway <- paste(sprintf("%02d", hour(MidwaySetTime)), 
                            sprintf("%02d", minute(MidwaySetTime)), 
                            sep=":")

data$Longitude..degs. <- as.numeric(data$Longitude..degs.)
data$Latitude..degs. <- as.numeric(data$Latitude..degs.)

########

# Run lotek_to_dataframe to convert your calculated riseset into a format usable by probGLS

trn <- lotek_to_dataframe( date = as.Date(as.character(data$Date), format = "%Y-%m-%d"), 
                           sunrise = data$RiseTimeUTC, 
                           sunset = data$SetTimeUTC, 
                           TRLon = data$Longitude..degs., 
                           TRLat = data$Latitude..degs.)
trn <- trn[!is.na(trn$tSecond),]
trn$keep <- loessFilter(trn, k = 3, plot = T) # run a Loess filter to take out unrealistic points. 
trn <- trn[which(trn$keep==T),]

############################################################################################################

# Plotting on-board algorithm calculated points - can compare final output to this as a sanity check

TRLon = data$Longitude..degs.
TRLat = data$Latitude..degs.
for (i in 1:length(TRLon)){
  if (TRLon[i]<0){
    TRLon[i] = TRLon[i]+360
  }
}
data(world2HiresMapEnv)
plot(TRLon,TRLat ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,240), ylim = c(-10, 75), bty = "n")
map('world2Hires', xlim= c(120,240), ylim = c(-10, 75),  col = "grey", add = T)
# connecting points - not necessarily helpful with visualization
# points(TRLon,TRLat, pch=16, col="cornflowerblue", type = "o")

############################################################################################################

# Prep additional data - wetdry and immersion temp ("sen" and "act" in algorithm)
# 'act' first

  ## datetime object needs to be in POSIXct, UTC time zone and must be called 'dtime'
  tl_data$dtime <- with(tl_data,paste0(Date,Time))
  # IMPORTANT! column (,5 or ,6) in the following line may change depending on what was recorded/programmed on the tag when deployed
  tl_data$dtime <- as.POSIXct(strptime(tl_data[,5], format="%m/%d/%Y %H:%M:%S"), tz="UTC")
  tl_data       <- tl_data[!is.na(tl_data$dtime),]
  act           <- tl_data
  
  ## wet dry data column must be called 'wetdry'; Lotek has 0 = wet and 1 = dry, but act expects
  ## 1 = wet and 0 = dry. 
  
  act$wetdry    <- 1-act$WetDryState
  act <- subset(act, select = c(dtime,wetdry))
  
  head(act)
  
  ## Can visualize proportion of time spent on water by day, code citation[1]
  ## 960 is the number of samples per day (40 90 sec intervals in an hour times 24 hours)
  ## Divide the tapply by the number of samples per day - 2700 32 sec intervals in 24 hours/112.5 32 sec intervals per hour)
  ## Theoretically will see a higher proportion after fledging
  
  divide_by <-(60*60*24)/wetdry.resolution
  plot(unique(as.Date(act$dtime)), tapply(act$wetdry,as.Date(act$dtime),sum)/divide_by,
       type="l",col=grey(0.5),ylab="daily proportion of saltwater immersion")
  points(unique(as.Date(act$dtime)),tapply(act$wetdry,as.Date(act$dtime),sum)/divide_by,
         pch=19,cex=0.8)

########

## Temperature data
## only keep temperature values when the logger was immersed in salt water 
## and if the 3 previous readings have been recorded while immersed in sea water as well
## citation[1]

  ## 2011 - tags appear to record a temperature every 160 seconds when "wet" is detected, for as long as wet is detected. 
  ## Will ingest all recorded values for sst_deduction, since it is likely this is capturing true wet periods. 
  
  IntTemp_data$dtime <- with(IntTemp_data,paste0(Date,Time))
  td <- IntTemp_data
  
  ## determine daily SST value recorded by the logger (takes the median temperature of all recordings in a day)
  sen           <- sst_deduction(datetime = td$dtime, temp = td$IntTemp, temp.range = c(-2,19))
  View(sen)
  
############################################################################################################

# Run algorithm

# download environmental data ----

# download yearly NetCDF files for (replace YEAR with appropriate number): 
# daily mean SST                   -> 'sst.day.mean.YEAR.v2.nc'
# daily SST error                  -> 'sst.day.err.YEAR.v2.nc'
# daily mean sea ice concentration -> 'icec.day.mean.YEAR.v2.nc'

# from:
# http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html
# and place all into the same folder

# Also, download the land mask file: 'lsmask.oisst.v2.nc' from the same directory 
# and place it in the same folder as all the other NetCDF files

# Specify final parameters for prob_algorithm 

tw    <- twilight_error_estimation() # Position estimates require error in twilight calculation; this baked-in function gives an error distribution for the algorithm 

pr   <- prob_algorithm(trn                         = trn, 
                       sensor                      = sen[sen$SST.remove==F,],
                       act                         = act, 
                       tagging.date                = start, 
                       retrieval.date              = end, 
                       loess.quartile              = NULL, # don't need to do this here because I did it above
                       tagging.location            = c(lon.calib,lat.calib), 
                       particle.number             = 800, 
                       iteration.number            = 60,
                       sunrise.sd                  = tw,
                       sunset.sd                   = tw,
                       range.solar                 = c(-7,-1),
                       boundary.box                = c(120,-120,-10,64), # 64 N for Laysan
                       days.around.spring.equinox  = c(10,10), 
                       days.around.fall.equinox    = c(10,10),
                       speed.dry                   = c(10,6,30), # Cite Weimerskirsch paper, Hyrenbach et al. 2002 for BFAL shows roughly similar to LAAL
                       speed.wet                   = c(1,1.3,5), 
                       sst.sd                      = 0.2, # logger-derived temperature accuracy    
                       max.sst.diff                = 3, # maximum discrepancy allowed between recorded and remote sensed SST        
                       east.west.comp              = F,
                       land.mask                   = T, 
                       ice.conc.cutoff             = 0, 
                       wetdry.resolution           = wetdry.resolution, # in seconds, i.e. 5 minutes = 300 seconds
                       NOAA.OI.location            = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/environment_data")

summary(pr)
plot_timeline(pr,degElevation = NULL, center.longitude =180)

plot_map(pr) # This function does not work for this data because BFAL and LAAL cross Pacific meridian (plot_map is Atlantic centric. Below is code for a Pacific-centric mapping). 
# Ultimately, unnecessary if importing into ArcMap, but can be helpful for data visualization. 

#######################################################################################################################################################

# Editing probGLS internal mapping function to display Pacific-centric view

# Change coordinates in probGLS output - only have to do this once before saving pr object (don't do it again on re-load)
convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}

pr[[1]]@coords[,1] <- convlon(pr[[1]]@coords[,1])
pr[[2]]@coords[,1] <- convlon(pr[[2]]@coords[,1])

p_map <- function (pr) 
{
  cc <- as.numeric(unlist(strsplit(as.character(pr[[4]][4,2]), "[ ]")))
  opar <- par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
  plot(pr[[1]], ylim=c(-10,75), xlim= c(120,240), col = "white", ylab = "Latitude", xlab = "Longitude")
  axis(1)
  axis(2)
  plot(pr[[1]], ylim=c(-10,75), xlim= c(120,240), col = colorRampPalette(c("grey90", "grey50"))(nrow(pr[[2]]))[pr[[1]]$step], 
       add = T, pch = 19, cex = 0.3)
  map('world2Hires', ylim=c(-10,75), xlim= c(120,240), add = T, col = 4, lwd = 0.5)
  mm2 <- pr[[2]]
  lines(mm2$lon, mm2$lat, col = 1)
  lines(c(mm2$lon[1], cc[1]), c(mm2$lat[1], cc[2]), lty = 3)
  lines(c(mm2$lon[nrow(mm2)], cc[1]), c(mm2$lat[nrow(mm2)], 
                                        cc[2]), lty = 3)
  points(mm2$lon, mm2$lat, cex = 0.7, pch = 21, bg = colorRampPalette(c("yellow", 
                                                                        "darkred"))(nrow(pr[[2]])))
  mm3 <- mm2[is.na(mm2$median.sun.elev), ]
  points(mm3$lon, mm3$lat, cex = 0.7, pch = 3)
  points(cc[1], cc[2], cex = 2.1, pch = 21, col = "violet", 
         lwd = 2)
  par(opar)
}

p_map(pr)

# save entire object so you can import and plot again later! Do this after convlon
save(pr, file = paste0(wd,"/data/daylog_work/spatial_points_dataframe_exports/pr_",Year,"_",
                       TrackNumber,".RData"))
# save image as 1440x900

############################################################################################################

## This script is used to identify, manually, the start and stop of post-breeding albatross migrations by inspecting SST timeseries. 
## Melinda Conners August 2003, adapted by Dallas Jordan Nov 2020
## This is a modified version of the 01b script that IDs all SST timeseries at once, for per-tag ID use

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
  
  f <- list.files(paste0("./data/daylog_work/temperature_logs_matched/",Species,"/",HatchYear,"_preprocessing"),pattern=paste0(Year,"_",TrackNumber,"_inttemp.csv")) # List of .csv files in directory - these should be your SST .txt logs, may need to change to .txt
  newm<-as.data.frame(matrix(NA,length(f),3)) # Create an empty matrix to fill 
  colnames(newm)<-c("file","start_day","end_day")
  
  for (i in 1:length(f)){
    
    # Read in file i , skip 3 lines due to headers
    # rawts<-read.table(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/temperature_logs_matched/",Species,"/",f[i],sep=""), skip=3, sep=",") 
    # colnames(rawts) <- c("index", "date","time", "SST","WetDry")
    # For 2008:
    rawts<-read.csv(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/temperature_logs_matched/",Species,"/",HatchYear,"_preprocessing/",f[1],sep=""), skip=3, sep=",", stringsAsFactors = F)
    colnames(rawts) <- c("index", "date","time", "SST")
    rawts$dtime <- with(rawts,paste0(date,time))
    rawts$dtime <- as.POSIXct(strptime(rawts[,5], format="%m/%d/%Y %H:%M:%S"), tz="UTC")
    rawts <- rawts %>%
      filter(dtime >= dateFilter1 & dtime <= dateFilter2)
    
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
    
    # Add to dataframe, dataframe exists because original code aggregated all sst .txt files in one dataframe
    newm$file[i]<-f[i]
    newm$start_day[i]<-cleanstart
    newm$end_day[i]<-cleanend
  }
  
  start_postbreeding <- cleanstart
  end_postbreeding <- cleanend
  
  # Grab dates to get postbreeding from pr object 
  # date(pr$`most probable track`$dtime) don't need this
  clean_start <- strptime(cleanstart,format="%m/%d/%Y", tz="UTC")
  clean_end <- strptime(cleanend,format="%m/%d/%Y", tz="UTC")
  # current_start <- date(clean_start)
  # current_end <- date(clean_end)
  
  # identify in plot_timeseries because tag SST memory filled up too soon
  current_start<-date("2011-07-06")
  current_end<-date("2011-11-24")
  

############################################################################################################
# Saving, exporting, and loading for re-plotting

# save geographic median .csv for all points in output object
most_probable <- pr[[2]]
most_probable_dtime <- as.data.frame(date(most_probable$dtime))
most_probable_lon <- as.data.frame(most_probable$lon)
most_probable_lat <- as.data.frame(most_probable$lat)
most_probable_export <- cbind(most_probable_dtime,most_probable_lon,most_probable_lat)
colnames(most_probable_export) <- c("dtime","Longitude", "Latitude")
# If you need to save this:
write.csv(most_probable_export,paste0(wd,"/data/daylog_work/geographic_median_exports/dl_",Year,"_",
                                      TrackNumber,".csv"), row.names = FALSE)    

# save geographic median .csv for just postbreeding points based on dates identified from SST analysis
pb_most_probable <- pr[[2]]
pb_most_probable_dtime <- as.data.frame(date(pb_most_probable$dtime))
pb_most_probable_lon <- as.data.frame(pb_most_probable$lon)
pb_most_probable_lat <- as.data.frame(pb_most_probable$lat)
pb_most_probable_export <- cbind(pb_most_probable_dtime,pb_most_probable_lon,pb_most_probable_lat)
colnames(pb_most_probable_export) <- c("dtime","Longitude", "Latitude")

pb_most_probable_export<- pb_most_probable_export %>%
  filter(dtime >= current_start & dtime <= current_end)

write.csv(pb_most_probable_export,paste0(wd,"/data/daylog_work/postbreeding_geographic_median_exports/pb_",Year,"_",
                                         TrackNumber,".csv"), row.names = FALSE)

# visualize postbreeding
plot(most_probable_export$Longitude,most_probable_export$Latitude ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,240), ylim = c(-10, 75), bty = "n", pch=16,cex=0.5)
points(pb_most_probable_export$Longitude,pb_most_probable_export$Latitude, pch=16, cex=0.5, col="red")
points(182.65,28.2,pch=7,cex=.8,col="darkblue")
map('world2Hires', xlim= c(120,240), ylim = c(-10, 75),  col = "grey", add = T)







# RE-LOADING TO PLOT/VIEW PROCESSED probGLS OBJECTS
# adds the object to the global environment, should be called "pr".

wd <- getwd()
Species <- "LAAL"
Year <- "1109"
TrackNumber <- "09" 
load(paste0(wd,"/data/daylog_work/spatial_points_dataframe_exports/",Species,"/pr_",Year,"_",
            TrackNumber,".RData"))

# check time span of pr object
head(pr$`most probable track`$dtime)
tail(pr$`most probable track`$dtime)
# In case you didn't do this before saving
convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}

pr[[1]]@coords[,1] <- convlon(pr[[1]]@coords[,1])
pr[[2]]@coords[,1] <- convlon(pr[[2]]@coords[,1])


p_map(pr) # may have to reload p_map function above!

