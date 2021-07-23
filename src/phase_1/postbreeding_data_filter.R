# Script for loading in Midway/Tern Island data and selecting only post-breeding data (July 1 - November 1)
# Dallas Jordan
# last edited: August 2020

############################################################################################################

# if you need to reinstall packages
install.packages("suncalc")
install.packages("remotes")
remotes::install_github("benjamin-merkel/probGLS", force= TRUE)
install.packages("mapdata")
library(devtools)
install_github("SLisovski/GeoLocTools")

# start here

library(maptools)
library(GeoLocTools)
setupGeolocation()
library(mapdata)
library(rgeos)
library(raster)
library(chron)
library(suncalc)
library(ncdf4)
library(lubridate)
library(dplyr)

# if you aren't working in your project 'spatial segregation'
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation")
filenames<- list.files()

# load in master raw data
raw_data_midway <- read.csv("./data/raw_data.csv", stringsAsFactors = F)
raw_data_tern <- read.csv("./data/Conners_TEIS_allmodeledGLS.csv", stringsAsFactors = F)
unique_midway <- unique(raw_data_midway$Bird.ID)
unique_tern <- unique(raw_data_tern$tripid)

# create list of partitioned bird_id dataframes

partition_midway <- vector("list", length(unique_midway))
for (i in 1:length(unique_midway)) {
  add <- as.data.frame(dplyr::filter(raw_data_midway, Bird.ID==unique_midway[i]))
  partition_midway[[i]] <- add
}

partition_tern <- vector("list", length(unique_tern))
for (i in 1:length(unique_tern)) {
  add <- as.data.frame(dplyr::filter(raw_data_tern, tripid==unique_tern[i]))
  partition_tern[[i]] <- add
}
############################################################################################################

# convert all times in all datasets to POSIX

for (i in 1:length(partition_tern)) {
  save <- mdy_hm(partition_tern[[i]]$GMT)
  partition_tern[[i]]$GMT <- save
}

for (i in 1:length(partition_midway)) {
  save <- mdy_hm(partition_midway[[i]]$Date.and.Time)
  partition_midway[[i]]$Date.and.Time <- save
}

# create postbreeding (July 1 - November 1) datasets

for (i in 1:length(partition_tern)) {
  partition_tern[[i]] <- partition_tern[[i]] %>% filter(month(partition_tern[[i]]$GMT) %in% (7:10))
}

for (i in 1:length(partition_midway)) {
  partition_midway[[i]] <- partition_midway[[i]] %>% filter(month(partition_midway[[i]]$Date.and.Time) %in% (7:10))
}

# export coordinate tracks
# make sure tern/midway_export folders are empty before you do this

for (i in 1:length(partition_tern)) {
  export <- data.frame(partition_tern[[i]]$tripid, partition_tern[[i]]$GMT, partition_tern[[i]]$xm.lon, partition_tern[[i]]$xm.lat)
  colnames(export) <- c("id","datetime","lon","lat")
  write.csv(x = export, file = 
              paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
                     paste0("tern_",export$id[1],"_",
                            year(export$datetime[1]),".csv")), row.names = FALSE)
}

for (i in 1:length(partition_midway)) {
  export <- data.frame(partition_midway[[i]]$Bird.ID, partition_midway[[i]]$Date.and.Time, 
                       partition_midway[[i]]$Raw.Longitude, partition_midway[[i]]$Raw.Latitude)
  colnames(export) <- c("id","datetime","lon","lat")
  write.csv(x = export, file = 
              paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/", 
                     paste0("midway_",export$id[1],"_",
                            year(export$datetime[1]),".csv")), row.names = FALSE)
}









































