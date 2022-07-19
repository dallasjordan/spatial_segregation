# Code to take revised Tern Island data and my previous Midway postbreeding data
# and reorganize, preserving trip ID and date. The idea is the same as data_reorganization.R
# but I'm a much better coder now
# Dallas Jordan
# July 18 2022
# Last updated July 18 2022


# Setup -------------------------------------------------------------------

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Fix Melinda's files that were not able to be processed in probGLS -------
# there were 3 of them...not making a reprex here, just trying to get them cleaned
# setwd("~/projects/spatial_segregation/files/data/final_reanalysis/postbreeding_tracks/csv/tern")
# fix <- read_csv("~/projects/spatial_segregation/files/data/final_reanalysis/postbreeding_tracks/csv/tern/t_pb_291204_EDIT.csv")
# 
# colnames(fix) <- c("tripid","year","date","time","Longitude","Latitude")
# fix <- fix %>% unite(dtime, c(date,time),sep=" ")
# fix$dtime <- as.POSIXct(strptime(fix$dtime, format="%Y-%m-%d %H:%M:%S"), tz="GMT")
# fix <- fix[,c(3:5)]
# write_csv(fix, file="t_pb_291204.csv")


# Midway ------------------------------------------------------------------
#COMBINE ALL INDIVIDUALS IN A YEAR OF BOTH SPECIES, PRESERVE INDIVIDUAL ID AND DTIME:

island = "midway"

setwd(paste0("/Users/dallasjordan/projects/spatial_segregation/files/data/final_reanalysis/postbreeding_tracks/csv/",island,"/"))

file_list <- list.files(path=getwd())
file_names <- str_extract(file_list, pattern="[^.csv]")
ids <- sub(".csv.*","",file_list)
ids <- str_extract(file_list,pattern="(\\d)+")
load_files <- function(file_list){
  data <- read_csv(file_list)
}
loaded <- lapply(file_list,load_files)
loaded_ids <-Map(cbind,loaded,tripid=ids)
midway_df <- bind_rows(loaded_ids)

island = "tern"

setwd(paste0("/Users/dallasjordan/projects/spatial_segregation/files/data/final_reanalysis/postbreeding_tracks/csv/",island,"/"))

file_list <- list.files(path=getwd())
file_names <- str_extract(file_list, pattern="[^.csv]")
ids <- sub(".csv.*","",file_list)
ids <- str_extract(file_list,pattern="(\\d)+")
load_files <- function(file_list){
  data <- read_csv(file_list)
}
loaded <- lapply(file_list,load_files)
loaded_ids <-Map(cbind,loaded,tripid=ids)
tern_df <- bind_rows(loaded_ids)

both_df <- rbind(midway_df,tern_df)


setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/")
save(LAAL, file="LAALdata_midway_withTrackID_dtime.Rdata")
save(BFAL, file="BFALdata_midway_withTrackID_dtime.Rdata")

# Tern --------------------------------------------------------------------


