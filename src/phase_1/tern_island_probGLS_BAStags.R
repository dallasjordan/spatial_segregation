# probGLS implementation for Melinda Conners's Tern Island albatross GLS data
# This script is used to process BAStag (.lig/.trn formats of daylog and templogs)
# Dallas Jordan
# last edited: June 5 2022
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

remotes::install_github("benjamin-merkel/probGLS")

# Setup -------------------------------------------------------------------

library(readr)
library(maptools)
library(mapdata)
library(rgeos)
library(raster)
library(ncdf4)
library(dplyr)
library(probGLS)
library(adehabitatHR)
library(lubridate)
library(stringr)


# Set info per tag --------------------------------------------------------

# info per-tag - these entries require manual changes referencing the "Daylog Matching" tab of 
# the Data Inventory spreadsheet. This should be the only section that requires manual changing
# in between tag analyses

Species <- "28"
Year <- "11"
TrackNumber <- "05" 

# for all Tern tags, the coordinates of Tern island colony

lat.calib <- 23.87
lon.calib <- 193.72 # correct Tern coords already
wetdry.resolution <- NULL # sampling rate of Basic Log in seconds, e.g. once every 5 min = 300 seconds
# Note that the Mk7 and Mk19 have a different activity recording strategy to this and 
# give higher time resolution; the penalty is that the memory fills quicker than with the 
# normal algorithm. These devices record the exact time (within 3secs) a state change occurs; 
# but the new state is recorded only if it is sensed for 6secs or more. This means there
# is no consistent resolution; it varies

# Read in data ------------------------------------------------------------

# Files read in from unzipped BAS tag files. Should be more straightforward than Lotek tags - you need a .trn for the rise/set 
# times and the .deg file for the temperature data. [possibly another for immersion data]


setwd(paste0("/Users/dallasjordan/projects/spatial_segregation/files/data/daylog_work/Tern/tern_island_data_for_processing/load_in/20",Year,"/"))

# Run trn_to_dataframe to convert .trn files into a format usable by probGLS

trn <- trn_to_dataframe(paste0("./",Species,Year,TrackNumber,".trn"))
trn <- trn[!is.na(trn$tSecond),]
trn$keep <- loessFilter(trn, k = 3, plot = T) # run a Loess filter to take out unrealistic points. 
trn <- trn[which(trn$keep==T),]
head(trn)
tail(trn)

# first/last date as noted in conners_metdata.xlsx

# "MM/D/YYYY or MM/DD/YYYY", no 0 for e.g. 05, for older daylogs
start    <- as.Date("2011-01-19")
end      <- as.Date("2011-12-24")

############################################################################################################

# Prep additional data - wetdry and immersion temp ("sst" and "act" in algorithm)
# BIG NOTE - for BAS Tags (and possibly Lotek?) the sst ("sen") parameter is derived
# from the .act file. It's not separate. Find Wetdry data, then filter out for 
# 3 previous "wet" states, use that SST as sst/sen...or look at the Lotek manual
# to see the difference between IntTemp and SST1 (external stem readings vs internal
# sensor temp readings)
# 'act' first

## datetime object needs to be in POSIXct, UTC time zone and must be called 'dtime'
tl_data <- read_csv(paste0("./",Species,Year,TrackNumber,".act"), 
               col_names = FALSE)
tl_data <- tl_data[,c(2,4,5)]
colnames(tl_data) <- c('Date / Time', 'duration', "wet.dry")
tl_data$dtime <- as.POSIXct(strptime(tl_data$`Date / Time`, format="%d/%m/%y %H:%M:%S"), tz="UTC")
tl_data       <- tl_data[!is.na(tl_data$dtime),]
act           <- tl_data
act           <- act[,c(4,2,3)]
head(act) 

## wet dry data column must be called 'wetdry'; Lotek has 0 = wet and 1 = dry, but act expects
## 1 = wet and 0 = dry. 
## For BAS tags, need to convert 'wet' to 1 and 'dry' to 0

# act <- act %>% mutate(wet.dry = 
#   case_when(wet.dry=="wet" ~ 1,
#             wet.dry=="dry" ~ 0)
# )
# 
# head(act)

# try: 
act <- act[act$wet.dry=="wet",]
act$wetdry    <- act$duration
head(act)
act <- subset(act, select = c(dtime,wetdry))
head(act)


act           <- BBA_deg[BBA_deg$wet.dry=="wet",]
act$wetdry    <- act$duration
########

## Immersion temperature data
## there needs to be a temp recording somewhere to do this - unfortunately it
## is not present in the data provided to me and probably wasn't recorded. 
## sen=NULL since there is no SST data

# td <- tl_data
# 
# ## determine daily SST value recorded by the logger (takes the median temperature of all recordings in a day)
# sen           <- sst_deduction(datetime = td$dtime, temp = td$`IntTemp [C]`, temp.range = c(-2,19))
# View(sen)

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
                       sensor                      = NULL,
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
                       boundary.box                = c(120,-90,-10,64), # 64 N for Laysan
                       days.around.spring.equinox  = c(10,10), 
                       days.around.fall.equinox    = c(10,10),
                       speed.dry                   = c(10,6,50), # Cite Weimerskirsch paper, Hyrenbach et al. 2002 for BFAL shows roughly similar to LAAL
                       speed.wet                   = c(1,1.3,5), 
                       sst.sd                      = 0.2, # logger-derived temperature accuracy    
                       max.sst.diff                = 3, # maximum discrepancy allowed between recorded and remote sensed SST        
                       east.west.comp              = F,
                       land.mask                   = T, 
                       ice.conc.cutoff             = 0, 
                       wetdry.resolution           = wetdry.resolution, # in seconds, i.e. 5 minutes = 300 seconds
                       NOAA.OI.location            = "/Users/dallasjordan/projects/spatial_segregation/files/environment_data")

summary(pr)
# plot_timeline(pr,degElevation = NULL, center.longitude =180)

# plot_map(pr) # This function does not work for this data because BFAL and LAAL cross Pacific meridian (plot_map is Atlantic centric. Below is code for a Pacific-centric mapping). 
# Ultimately, unnecessary if importing into ArcMap, but can be helpful for data visualization. 

#######################################################################################################################################################

# Editing probGLS internal mapping function to display Pacific-centric view

# Change coordinates in probGLS output - only have to do this once before saving pr object (don't do it again on re-load)
# As of Feb 2021, I did have to run this again for LAAL 2008 - not sure why? Does it apply to others?
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

# save entire object so you can import and plot again later!
save(pr, file = paste0("/Users/dallasjordan/projects/spatial_segregation/files/data/daylog_work/Tern/tern_spatial_points_dataframe_exports/pr_",Species,Year,
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


## Identification of post-breeding period ------------------------------------------------------------------------------------------------

wd <- setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/Tern/tern_island_data_for_processing/temp_logs/")
f <- list.files(wd,pattern=paste0(Species,Year,TrackNumber)) # List of .csv files in directory - these should be your SST .txt logs, may need to change to .txt
newm<-as.data.frame(matrix(NA,length(f),3)) # Create an empty matrix to fill 
colnames(newm)<-c("file","start_day","end_day")

for (i in 1:length(f)){
  
  # Read in file i , skip 3 lines due to headers
  # rawts<-read.csv(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/temperature_logs_matched/",Species,"/",f[i],sep=""), skip=3, sep=",")
  # colnames(rawts) <- c("index", "date","time", "SST","WetDry")
  # For 2008:
  tl_data <- read_csv(paste0(Species,Year,TrackNumber,"_templog.TXT"), skip = 2)
  rawts<-read_csv(f[i], skip=2)
  colnames(rawts) <- c("index", "date", "SST","WetDry")
  # rawts<-rawts[min(which(grepl(actual_first,rawts$date))):max(which(grepl(actual_last,rawts$date))),]
  
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
current_start <- date(clean_start)
current_end <- date(clean_end)

############################################################################################################

# Saving, exporting, and loading for re-plotting

# save geographic median .csv for all points in output object
most_probable <- pr[[2]]
most_probable_dtime <- as.data.frame(most_probable$dtime)
most_probable_lon <- as.data.frame(most_probable$lon)
most_probable_lat <- as.data.frame(most_probable$lat)
most_probable_export <- cbind(most_probable_dtime,most_probable_lon,most_probable_lat)
colnames(most_probable_export) <- c("dtime","Longitude", "Latitude")
# If you need to save this:
write.csv(most_probable_export,paste0("/Users/dallasjordan/projects/spatial_segregation/files/data/daylog_work/Tern/tern_geographic_median_exports/t_dl_",Species,Year,
                                      TrackNumber,".csv"), row.names = FALSE)    

# save geographic median .csv for just postbreeding points based on dates identified from SST analysis
pb_most_probable <- pr[[2]]
pb_most_probable_dtime <- as.data.frame(pb_most_probable$dtime)
pb_most_probable_lon <- as.data.frame(pb_most_probable$lon)
pb_most_probable_lat <- as.data.frame(pb_most_probable$lat)
pb_most_probable_export <- cbind(pb_most_probable_dtime,pb_most_probable_lon,pb_most_probable_lat)
colnames(pb_most_probable_export) <- c("dtime","Longitude", "Latitude")

# filter by using known postbreeding dates that Melinda identified using the SST variability method instead of running ID portion of script again:
TrackNumber  
cleanstart <- "06/08/2011"
cleanend <- "11/24/2011" 
clean_start <- strptime(cleanstart,format="%m/%d/%Y", tz="UTC")
clean_end <- strptime(cleanend,format="%m/%d/%Y", tz="UTC")
current_start <- date(clean_start)
current_end <- date(clean_end)

pb_most_probable_export<- pb_most_probable_export %>%
  filter(dtime >= current_start & dtime <= current_end)

#if you need to save this: 
write.csv(pb_most_probable_export,paste0("/Users/dallasjordan/projects/spatial_segregation/files/data/daylog_work/Tern/tern_postbreeding_geographic_median_exports/t_pb_",Species,Year,
                                         TrackNumber,".csv"), row.names = FALSE)



# Visualize final data ----------------------------------------------------

# visualize full track 
plot(most_probable_export$Longitude,most_probable_export$Latitude ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,270), ylim = c(-10, 75), bty = "n", pch=16,cex=0.5)
points(most_probable_export$Longitude,most_probable_export$Latitude, pch=16, cex=0.5, col="red")
points(182.65,28.2,pch=7,cex=.8,col="darkblue")
map('world2Hires', xlim= c(120,270), ylim = c(-10, 75),  col = "grey", add = T)
# visualize postbreeding
plot(pb_most_probable_export$Longitude,pb_most_probable_export$Latitude ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,270), ylim = c(-10, 75), bty = "n", pch=16,cex=0.5)
points(pb_most_probable_export$Longitude,pb_most_probable_export$Latitude, pch=16, cex=0.5, col="red")
points(182.65,28.2,pch=7,cex=.8,col="darkblue")
map('world2Hires', xlim= c(120,270), ylim = c(-10, 75),  col = "grey", add = T)





# Other functions/processes that might be helpful -------------------------


# RE-LOADING TO PLOT/VIEW PROCESSED probGLS OBJECTS
# adds the object to the global environment, should be called "pr".

wd <- getwd()
Species <- "BFAL"
Year <- "2112"
TrackNumber <- "03" 
load(paste0(wd,"/data/daylog_work/spatial_points_dataframe_exports/",Species,"/pr_",Year,"_",
            TrackNumber,".RData"))

# Used for Feb 2021 work, where I needed to go back and capture full datetime for each postbreeding export:

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
pr[[1]]@coords[,1] <- convlon(pr[[1]]@coords[,1])
pr[[2]]@coords[,1] <- convlon(pr[[2]]@coords[,1])

cleanstart <- "06/24/2012"
cleanend <- "11/22/2012" # should be day after the last date listed
clean_start <- strptime(cleanstart,format="%m/%d/%Y", tz="UTC")
clean_end <- strptime(cleanend,format="%m/%d/%Y", tz="UTC")
current_start <- date(clean_start)
current_end <- date(clean_end)

# check time span of pr object
head(pr$`most probable track`$dtime)
tail(pr$`most probable track`$dtime)

p_map(pr) # may have to reload p_map function above!



# Going to run multi-year postbreeding tracks twice and save outputs differently!


# Edge case when you have tracks that span multiple postbreeding --------


# save geographic median .csv for all points in output object
most_probable <- pr[[2]]
most_probable_dtime <- as.data.frame(date(most_probable$dtime))
most_probable_lon <- as.data.frame(most_probable$lon)
most_probable_lat <- as.data.frame(most_probable$lat)
most_probable_export <- cbind(most_probable_dtime,most_probable_lon,most_probable_lat)
colnames(most_probable_export) <- c("dtime","Longitude", "Latitude")

most_probable_export<- most_probable_export %>%
  filter(dtime >= current_start & dtime < current_end)
# filter(dtime <= as.Date(paste0(year(current_start),"-12-15")))
# filter(dtime > as.Date(paste0((year(current_start)-1),"-12-15"))) # which filter function you need depends on if you are doing the first or second postbreeding of a multiyear tag
# If you need to save this (a or b track):
write.csv(most_probable_export,paste0(wd,"/data/daylog_work/geographic_median_exports/dl_",Year,"_",
                                      TrackNumber,"_",year(current_start),".csv"), row.names = FALSE) 

# save geographic median .csv for just postbreeding points based on dates identified from SST analysis
pb_most_probable <- pr[[2]]
pb_most_probable_dtime <- as.data.frame(date(pb_most_probable$dtime))
pb_most_probable_lon <- as.data.frame(pb_most_probable$lon)
pb_most_probable_lat <- as.data.frame(pb_most_probable$lat)
pb_most_probable_export <- cbind(pb_most_probable_dtime,pb_most_probable_lon,pb_most_probable_lat)
colnames(pb_most_probable_export) <- c("dtime","Longitude", "Latitude")

# pb_most_probable_export<- pb_most_probable_export %>%
#  filter(dtime >= current_start & dtime <= current_end)
pb_most_probable_export<- pb_most_probable_export %>%
  filter(dtime >= current_start & dtime <= as.Date(paste0((year(current_start)),"-12-15")))

#write.csv(pb_most_probable_export,paste0(wd,"/data/daylog_work/postbreeding_geographic_median_exports/pb_",Year,"_",
#                                         TrackNumber,"_",year(current_start),".csv"), row.names = FALSE)
write.csv(pb_most_probable_export,paste0(wd,"/data/daylog_work/postbreeding_geographic_median_exports/pb_",Year,"_",
                                         TrackNumber,"_",year(current_start),".csv"), row.names = FALSE)

# visualize postbreeding
plot(most_probable_export$Longitude,most_probable_export$Latitude ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,240), ylim = c(-10, 75), bty = "n", pch=16,cex=0.5)
points(pb_most_probable_export$Longitude,pb_most_probable_export$Latitude, pch=16, cex=0.5, col="red")
points(182.65,28.2,pch=7,cex=.8,col="darkblue")
map('world2Hires', xlim= c(120,240), ylim = c(-10, 75),  col = "grey", add = T)


# If you want to replot the comparison points:  ---------------------------
most_probable_export <- read.csv(paste0("data/daylog_work/geographic_median_exports/",Species,"/dl_",Year,"_",TrackNumber,".csv"))
pb_most_probable_export <- read.csv(paste0("data/daylog_work/postbreeding_geographic_median_exports/",Species,"/pb_",Year,"_",TrackNumber,".csv"))

plot(most_probable_export$Longitude,most_probable_export$Latitude ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,240), ylim = c(-10, 75), bty = "n", pch=16,cex=0.5)
points(pb_most_probable_export$Longitude,pb_most_probable_export$Latitude, pch=16, cex=0.5, col="red")
points(182.65,28.2,pch=7,cex=.8,col="darkblue")
map('world2Hires', xlim= c(120,240), ylim = c(-10, 75),  col = "grey", add = T)
