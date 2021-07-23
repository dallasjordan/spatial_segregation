# probGLS implementation for Melinda Conners's Tern Island albatross GLS data
# This script is used for older daylog files with "m/d/y" formatting. 
# Dallas Jordan
# last edited: July 2021
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
library(chron)
library(suncalc)
library(ncdf4)
library(dplyr)
library(probGLS)
library(adehabitatHR)
library(lubridate)


# Set info per tag --------------------------------------------------------

# info per-tag - these entries require manual changes referencing the "Daylog Matching" tab of 
# the Data Inventory spreadsheet. This should be the only section that requires manual changing
# in between tag analyses

Species <- "29"
Year <- "08"
TrackNumber <- "01" 

# first/last date as noted in data_inventory

# "MM/D/YYYY or MM/DD/YYYY", no 0 for e.g. 05, for older daylogs
start    <- as.Date("2012-07-11")
end      <- as.Date("2013-05-01")

# for all Tern tags, the coordinates of Tern island colony

lat.calib <- 23.87
lon.calib <- 193.72 # correct Tern coords already
wetdry.resolution <- 100 # sampling rate of Basic Log in seconds, e.g. once every 5 min = 300 seconds


# Read in data ------------------------------------------------------------

# CSVs loaded in here. Note: each of these csv files were converted from raw .bin files using Lotek 
# TagTalk software. The converted csv files had to rows above the data, stating "GMT Time Correction = 0" 
# and either "Lg00 (Day Log) Log:" or "Lg20 (Basic Log) Log:" depending on the type of log. These two rows
# have been deleted to make read-in easier. 

wd <- setwd("E:/project_data/spatial_segregation/tern_island_data/")

# load in daylog
data <- read.csv(paste0("data/daylog_work/daylogs_matched/",Species,"/","dl_",Year,"_",TrackNumber, ".csv"),
                 sep="", header = T,row.names = NULL, stringsAsFactors = F, skip = 2) # may need to change sep=","

# For 2008 data: 
colnames(data) <- c("Rec","Date","Time","Sunrise","Sunset","WetDryChange","SST1[C]","Latitude..degs.","Longitude..degs.")
data <- data[,1:9]

# load in additional info log (wetdry/immersion, called "temperature logs")
tl_data <- read.csv(paste0("data/daylog_work/temperature_logs_matched/",Species,"/tl_",Year,"_",TrackNumber, ".csv"),
                    sep="", header = T,row.names = NULL, stringsAsFactors = F, skip = 2) # may need to change sep=","

# For 2008 data: 
colnames(tl_data) <- c("Rec","Date","Time","IntTemp[C]","WetDryState")
tl_data <- tl_data[,1:5]

# trim to first/last date as noted in data_inventory

data <- data[min(which(grepl(actual_first,data$Date))):max(which(grepl(actual_last,data$Date))),]
####!!!!!!!!#######!!!!!!!######### DO NOT RUN IF ABOVE LINE WORKS
# special case if the time series ended early (tag wasn't retrieved within a year and
# recorded too much)/if last date recorded is before last day of tag's deployment (e.g. memory is full)
data <- data[min(which(grepl(actual_first,data$Date))):nrow(data),]


tl_data <- tl_data[min(which(grepl(actual_first,tl_data$Date))):max(which(grepl(actual_last,tl_data$Date))),]
####!!!!!!!!#######!!!!!!!######### DO NOT RUN IF ABOVE LINE WORKS
# special case if the time series ended early (tag wasn't retrieved within a year and
# recorded too much)/if last date recorded is before last day of tag's deployment (e.g. memory is full)
tl_data <- tl_data[min(which(grepl(actual_first,tl_data$Date))):nrow(tl_data),]

# MAYBE DONT NEED THIS remove rows with NA (200) values for longitude
# data <- data[-which(grepl(200,data$Longitude..degs.)),]

# visualize daylight period
RiseMinute <- data$Sunrise
RiseMinute <- RiseMinute[RiseMinute<2500] # excludes NA values
SetMinute <- data$Sunset
SetMinute <- SetMinute[SetMinute<2500] # excludes NA values
date_plot <- strptime(as.character(data$Date), format="%m/%d/%Y", tz="GMT") # might need to change the Y to a y if graph looks weird
date_plot <- as.POSIXct(date_plot)
plot(RiseMinute ~ date_plot, type = "n", ylim = range(c(RiseMinute, SetMinute)), ylab = "tag day-time")
polygon(c(date_plot, rev(date_plot)),c(RiseMinute, rev(SetMinute)), col = "yellow")
abline(h = 1440) # 24 hour period, 1440 minutes = 24 hours
title("daylight period (minute of the day in UTC)")


# Lotek to dataframe ------------------------------------------------------

# lotek_to_dataframe to convert your calculated daylog into a format usable by probGLS by back-calculating sunrise and sunset times. The Lotek
# hard-coded sun elevation angle is -3.44. It uses this to generate sunrise/sunset times. 
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

# For daylogs that already have times in HH:MM:
# data$Sunrise <- gsub("(\\d\\d)(\\d\\d)", "\\1:\\2", data$Sunrise)
# data$Sunrise <- gsub("(\\d)(\\d\\d)", "0\\1:\\2", data$Sunrise)
# data$Sunset <- gsub("(\\d\\d)(\\d\\d)", "\\1:\\2", data$Sunset)
# data$Sunset <- gsub("(\\d)(\\d\\d)", "0\\1:\\2", data$Sunset)

# Run lotek_to_dataframe to convert your calculated riseset into a format usable by probGLS

trn <- lotek_to_dataframe( date = as.Date(as.character(data$Date), format = "%Y-%m-%d"), 
                           sunrise = data$RiseTimeUTC, 
                           sunset = data$SetTimeUTC, 
                           TRLon = data$Longitude..degs., 
                           TRLat = data$Latitude..degs.)
trn <- trn[!is.na(trn$tSecond),]
trn$keep <- loessFilter(trn, k = 3, plot = T) # run a Loess filter to take out unrealistic points. 
trn <- trn[which(trn$keep==T),]


# Plot lotek-calculated points --------------------------------------------

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


