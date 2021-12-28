# SST boxplots, taken from basic logs of data from Scott
# Dallas Jordan
# last updated oct 21 2021

# new challenge: have to load in postbreeding tracks, get date range, cull any SST data outside that date range


# Setup -------------------------------------------------------------------

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

read.file <- function (file.name) {
  file <- try(read.table(file.name))
  if (class(file) == 'try-error') {
    cat('Caught an error during read.table, trying read_csv.\n')
    file <- as.data.table(read_csv(file.name, sep = " ", quote = ""))
  }
  file
}


# LOTEK -------------------------------------------------------------------

# load in data 

# info per-tag - these entries require manual changes 

Species <- "BFAL"
Year <- "2008"

setwd(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',Year,'/',Species,"/"))

file.list <- list.files(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',Year,"/",Species,"/"))

# for non-2008 data: 
N <- length(file.list)
big.df <- data.frame()
pb.file.list <- list.files(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/',Species,"/",Year))
for (i in 1:length(file.list)){
  add.file <- read_csv(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',Year,'/',Species,'/',file.list[i]), skip=2)
  pb.find.dates <- read_csv(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/',Species,"/",Year,"/",pb.file.list[i]))
  min.date <- min(pb.find.dates$dtime)
  max.date <- max(pb.find.dates$dtime)
  add.file$Date <- parse_date_time(add.file$Date, c("mdY", "mdy","dmY","dmy"))
  add.file <- add.file %>% filter(Date>min.date)
  add.file <- add.file %>% filter(Date<max.date)
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}
big.df <- big.df[,c(1,2,3,7)]
colnames(big.df) <- c("Rec","date","time","SST")
big.df$SST <- as.numeric(big.df$SST)
big.df <- na.omit(big.df)

#for 2008 data (this one is really messed up, have to change for LAAL and BFAL):
N <- length(file.list)
big.df <- data.frame()
pb.file.list <- list.files(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/',Species,"/",Year))
for (i in 1:length(file.list)){
  add.file <- read_csv("dl_2108_01.csv",skip = 2)
  colnames(add.file) <- c("Rec","Date","Time","Sunrise","Sunset","WetDryChange","SST1[C]","Latitude..degs.","Longitude..degs.")
  # add.file <- add.file[,1:9]
  pb.find.dates <- read_csv(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/',Species,"/",Year,"/",pb.file.list[i]))
  min.date <- min(pb.find.dates$dtime)
  max.date <- max(pb.find.dates$dtime)
  add.file$Date <- parse_date_time(add.file$Date, c("mdY", "mdy", "mdY HM", "mdy HM"))
  add.file <- add.file %>% filter(Date>min.date)
  add.file <- add.file %>% filter(Date<max.date)
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}
big.df <- big.df[,c(1,2,3,6)]
colnames(big.df) <- c("Rec","date","time","SST")


# resume: 

# remove NA values, '1290.70' values which indicate no reading in Lotek tags. <35 and >-3 was set by visually inspecting
# logs and seeing what "error" values are - there were many "38.58" values indicating this was an error code
big.df <- big.df %>% filter(SST<35) # filter SST to be less than 35C
big.df <- big.df %>% filter(SST>-3) # filter SST to be greater than -3

write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/cleaned/filtered_for_postbreeding/",Species,"/",Year,'_',Species,".csv"),row.names = FALSE)



