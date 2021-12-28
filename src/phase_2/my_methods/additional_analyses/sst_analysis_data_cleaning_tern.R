# SST boxplots, taken from basic logs of data from Melinda
# Dallas Jordan
# last updated oct 21 2021

# tasks
# pull up data inventory, check the google drive logs vs whats in data inventory
# two steps - lotek, then BAS
# load in the right file for each Tern bird
# parse out sst
# save to a new data frame, preserve year and species and island


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


setwd(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/sst_matched/',Year,'/',Species,"/"))

file.list <- list.files(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/sst_matched/',Year,"/",Species,"/"))

# N <- length(file.list)
# data.list <- vector("list",N)
# for (i in 1:length(file.list)){
#   add.file <- read.csv(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/',Year,"/",Species,"/",file.list[i]), header=FALSE, stringsAsFactors=TRUE)
#   data.list[[i]] <-add.file
# }

N <- length(file.list)
big.df <- data.frame()
for (i in 1:length(file.list)){
  # for BFAL 2010
  #add.file <- read.table(paste0("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/sst_matched/",Year,"/",Species,"/",file.list[i]), quote="\"", comment.char="", skip=3)
  # not sure for what: 
  # add.file <- read.table(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/sst_matched/',Year,"/",Species,"/",file.list[i]), header=FALSE, stringsAsFactors=FALSE, skip=3)
  # for everything else: 
  add.file <- read.table(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/sst_matched/',Year,"/",Species,"/",file.list[i]),quote="\"",comment.char="",stringsAsFactors=FALSE, skip=3)
  add.file$id <- paste0("bt",Year,"_",i)
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

# for 2008 BFAL: 
colnames(big.df) <- c("Rec","date","time","Sunrise","Sunset","WetDryChange","SST","Latitude","Longitude","id")
big.df$time <- gsub(",","",big.df$time)
big.df$SST <- gsub(",","",big.df$SST)
big.df$Sunrise <- gsub(",","",big.df$Sunrise)
big.df$Sunset <- gsub(",","",big.df$Sunset)
big.df$WetDryChange <- gsub(",","",big.df$WetDryChange)
big.df$Latitude <- gsub(",","",big.df$Latitude)
big.df$SST <- as.numeric(big.df$SST)

#for 2008 LAAL: 
# colnames(big.df) <- c("Rec","date","time","Sunrise","Sunset","WetDryChange","SST","Latitude","Longitude","id")
# big.df$time <- gsub(",","",big.df$time)
# big.df$SST <- gsub(",","",big.df$SST)
# big.df$Sunrise <- gsub(",","",big.df$SST)
# big.df$Sunset <- gsub(",","",big.df$Sunset)
# big.df$WetDryChange <- gsub(",","",big.df$WetDryChange)
# big.df$Latitude <- gsub(",","",big.df$Latitude)
# big.df$SST <- as.numeric(big.df$SST)

# for 2009 LAAL:
# big.df <- big.df[,-c(6:15)]
# colnames(big.df) <- c("Rec","date","time","Sunrise","Sunset","WetCounter","SST","Latitude","Longitude","id")
# big.df$SST <- gsub(",","",big.df$SST)
# big.df$Sunrise <- gsub(",","",big.df$SST)
# big.df$Sunset <- gsub(",","",big.df$Sunset)
# big.df$WetCounter <- gsub(",","",big.df$WetCounter)
# big.df$Latitude <- gsub(",","",big.df$Latitude)
# big.df$SST <- as.numeric(big.df$SST)
# class(big.df$SST)

# for 2009 BFAL:
# big.df <- big.df[,-c(5:15)]
# colnames(big.df) <- c("Rec","dtime","Sunrise","Sunset","SST","WetCounter","Latitude","Longitude","id")
# big.df <- big.df[,c("Rec","dtime","Sunrise","Sunset","WetCounter","SST","Latitude","Longitude","id")]
# big.df$SST <- gsub(",","",big.df$SST)
# big.df$Sunrise <- gsub(",","",big.df$SST)
# big.df$Sunset <- gsub(",","",big.df$Sunset)
# big.df$WetCounter <- gsub(",","",big.df$WetCounter)
# big.df$Latitude <- gsub(",","",big.df$Latitude)
# big.df$SST <- as.numeric(big.df$SST)
# class(big.df$SST)

# for 2010: 
# colnames(big.df) <- c("Rec","date","time","Sunrise","Sunset","WetDryChange","SST","Latitude","Longitude","id")
# big.df$SST <- gsub(",","",big.df$SST)
# big.df$Sunrise <- gsub(",","",big.df$SST)
# big.df$Sunset <- gsub(",","",big.df$Sunset)
# big.df$WetDryChange <- gsub(",","",big.df$WetDryChange)
# big.df$Latitude <- gsub(",","",big.df$Latitude)
# big.df$SST <- as.numeric(big.df$SST)

# remove NA values, '1290.70' values which indicate no reading in Lotek tags. <35 and >-3 was set by visually inspecting
# logs and seeing what "error" values are - there were many "38.58" values indicating this was an error code
big.df <- big.df %>% filter(SST<35) # filter SST to be less than 35C
big.df <- big.df %>% filter(SST>-3) # filter SST to be greater than -3

# for LAAL 2009 i'm removing 0 values - they correlate to tag failures and look highly suspicious
# big.df <- big.df %>% filter(SST!=0)

# LAAL 2010 SST recording failed; not processed

write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/collated/",Year,'_',Species,".csv"),row.names = FALSE)

      
# manual addition for problem file - 2010 BFAL
add_manual <- read_csv("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/sst_matched/2010/BFAL/gls_2010_2810001.csv",
                                           col_names = FALSE, skip = 2)
add_manual <- add_manual[-c(1),]
add_manual$filler <- NA
add_manual <- add_manual[,c(1,2,9,3,4,5,6,7,8)]
add_manual$id <- "bt2010_5"
colnames(add_manual) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","id")
big.df<- rbind(big.df,add_manual)


# BAS ---------------------------------------------------------------------

# load in data 

# info per-tag - these entries require manual changes 

Species <- "BFAL"
Year <- "2012"


setwd(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/',Year,'/',Species,"/"))

file.list <- list.files(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/',Year,"/",Species,"/"))

# N <- length(file.list)
# data.list <- vector("list",N)
# for (i in 1:length(file.list)){
#   add.file <- read.csv(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/',Year,"/",Species,"/",file.list[i]), header=FALSE, stringsAsFactors=TRUE)
#   data.list[[i]] <-add.file
# }

N <- length(file.list)
big.df <- data.frame()
for (i in 1:length(file.list)){
  add.file <- read.csv(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/',Year,"/",Species,"/",file.list[i]), header=FALSE)
  add.file$id <- add.file$id <- paste0("bt",Year,"_",i)
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

big.df <- big.df[,c(1,2,4,5)]
colnames(big.df) <- c("Rec","dtime","SST","id")

# remove NA values, '1290.70' values which indicate no reading in Lotek tags. <35 and >-3 was set by visually inspecting
# logs and seeing what "error" values are - there were many "38.58" values indicating this was an error code
big.df <- big.df %>% filter(SST<35) # filter SST to be less than 35C
big.df <- big.df %>% filter(SST>-3) # filter SST to be greater than -3

write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/collated/",Year,'_',Species,".csv"),row.names = FALSE)

# Merging files and standardizing -----------------------------------------
# Tern 
Species <- "BFAL"
years <- c(2008,2009,2010,2011,2012)
setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/collated/',Species,"/"))

file.list <- list.files(getwd())

csvList <- list()
for (i in 1:length(file.list)){
  add<-read.csv(file.list[i])
  csvList[[i]]<-add
}

bt_2008<-csvList[[1]]
bt_2009<-csvList[[2]]
bt_2010<-csvList[[3]]
bt_2011<-csvList[[4]]
bt_2012<-csvList[[5]]

bt_2008 <- bt_2008[,c(1,2,7,10)]
bt_2009 <- bt_2009[,c(1,2,6,9)]
bt_2010 <- bt_2010[,c(1,2,7,10)]

colnames(bt_2008)<- c("Rec","dtime","SST","id")
colnames(bt_2010)<- c("Rec","dtime","SST","id")

all_BFAL_SST_Tern <- rbind(bt_2008,bt_2009,bt_2010,bt_2011,bt_2012)

Species <- "LAAL"
years <- c(2008,2009,2010,2011,2012)
setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/collated/',Species,"/"))

file.list <- list.files(getwd())

csvList <- list()
for (i in 1:length(file.list)){
  add<-read.csv(file.list[i])
  csvList[[i]]<-add
}

lt_2008<-csvList[[1]]
lt_2009<-csvList[[2]]
lt_2010<-csvList[[3]]
lt_2011<-csvList[[4]]
lt_2012<-csvList[[5]]

lt_2008 <- lt_2008[,c(1,2,7,10)]
lt_2009 <- lt_2009[,c(1,2,7,10)]
lt_2010 <- lt_2010[,c(1,2,7,10)]

colnames(lt_2008)<- c("Rec","dtime","SST","id")
colnames(lt_2009)<- c("Rec","dtime","SST","id")
colnames(lt_2010)<- c("Rec","dtime","SST","id")

all_LAAL_SST_Tern <- rbind(lt_2008,lt_2009,lt_2010,lt_2011,lt_2012)

# old, structure of collated files changed, rendering this useless
# big.df <- data.frame()
# for (i in 1:length(file.list)){
#   add.file <- read_csv(file.list[i])
#   add.file$year <- years[i] 
#   add.file$island <- "Tern"
#   add.file$spp <- Species
#   if (i==1){big.df <- add.file
#   } else {big.df <- rbind(big.df, add.file)}
# }

write.csv(all_LAAL_SST_Tern,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/collated/allLAAL_SST_Tern.csv"),row.names = FALSE)
