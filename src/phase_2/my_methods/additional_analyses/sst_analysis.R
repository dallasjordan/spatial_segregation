# SST boxplots, taken from basic logs of data from Scott
# Dallas Jordan
# last updated Nov 12 2021

# new challenge: cycle by day, average SST into a Daily Average, append to movement data


# Setup -------------------------------------------------------------------

library(readr)
library(maptools)
# library(GeoLocTools)
# setupGeolocation()
library(mapdata)
library(rgeos)
library(raster)
# library(chron)
# library(suncalc)
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


# data cleaning of Midway sst files (tern done in a seperate script) ------------------------------------
# do once! You won't want to repeat this anyway, 
# you manually did elements of a for-loop because none of the yearly files match columns...

#### Create one big csv of all LAAL midway daylogs, drop invalid values

Species <- "LAAL"
Year <- 2009

big.df <- data.frame()
sst.file.list <- list.files(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',
                                   Year,"/",Species,"/"))
for (i in 1:length(sst.file.list)){
    add.file <- read_csv(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',Year,'/',Species,'/',sst.file.list[i]), skip=2)
    add.file$Date <- parse_date_time(add.file$Date, c("mdY", "mdy","dmY","dmy"))
    add.file$id <- paste0("lm",Year,"_",i)
    if (i==1){big.df <- add.file
    } else {big.df <- rbind(big.df, add.file)}
}

# for 2008 files: 
Species <- "LAAL"
Year <- 2008

big.df <- data.frame()
sst.file.list <- list.files(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',
                                   Year,"/",Species,"/"))
for (i in 1:length(sst.file.list)){
  print(i)
  add.file <- read.delim(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',Year,'/',Species,'/',sst.file.list[i]), 
                         header=F,skip=3)
  add.file$V2 <- parse_date_time(add.file$V2, c("mdY", "mdy", "mdY HM", "mdy HMS", "mdY HMS"))
  add.file$id <- paste0("lm",Year,"_",i)
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

big.df <- big.df[,c(1,2,3,6,9)]
colnames(big.df) <- c("Rec","date","time","SST", "id")

big.df.2012 <- big.df
big.df.2011 <- big.df
big.df.2010 <- big.df
big.df.2009 <- big.df
big.df.2008 <- big.df

# rbind all years then clean: 
big.df <- rbind(big.df.2008,big.df.2009,big.df.2010,big.df.2011,big.df.2012)
big.df <- na.omit(big.df)
save(big.df, file="checkpoint")
big.df <- big.df %>% filter(SST>-4)
big.df <- big.df %>% filter(SST<35)

all_LAAL_midway_SST <- big.df
save(all_LAAL_midway_SST, file="/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/cleaned/all_LAAL_midway_SST.Rdata")

#### Repeat for BFAL at Midway: 

Species <- "BFAL"
Year <- 2009

big.df <- data.frame()
sst.file.list <- list.files(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',
                                   Year,"/",Species,"/"))
for (i in 1:length(sst.file.list)){
  add.file <- read_csv(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',Year,'/',Species,'/',sst.file.list[i]), skip=2)
  add.file$Date <- parse_date_time(add.file$Date, c("mdY", "mdy","dmY","dmy"))
  add.file$id <- paste0("bm",Year,"_",i)
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

# for 2008 files: 
Species <- "BFAL"
Year <- 2008

big.df <- data.frame()
sst.file.list <- list.files(paste0('/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',
                                   Year,"/",Species,"/"))
for (i in 1:length(sst.file.list)){
  print(i)
  if (i==1){add.file <- read.csv(paste0("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/",Year,"/",Species,"/",sst.file.list[i]), 
                                 skip=3,header=FALSE)} else {add.file <- read.delim(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/',Year,'/',Species,'/',sst.file.list[i]), 
                                                                                    header=F,skip=3) }
  
  add.file$V2 <- parse_date_time(add.file$V2, c("mdY", "mdy", "mdY HM", "mdy HMS", "mdY HMS","mdy HM"))
  add.file$id <- paste0("bm",Year,"_",i)
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

big.df <- big.df[,c(1,2,3,6,9)]
colnames(big.df) <- c("Rec","date","time","SST","id")

big.df.2012 <- big.df
big.df.2011 <- big.df
big.df.2010 <- big.df
big.df.2009 <- big.df
big.df.2008 <- big.df

# rbind all years then clean: 
big.df <- rbind(big.df.2008,big.df.2009,big.df.2010,big.df.2011,big.df.2012)
big.df <- na.omit(big.df)
save(big.df, file="checkpoint")
big.df <- big.df %>% filter(SST>-4)
big.df <- big.df %>% filter(SST<35)

all_BFAL_midway_SST <- big.df
save(all_BFAL_midway_SST, file="/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/cleaned/all_BFAL_midway_SST.Rdata")

#### Merge and standardize (is this section still needed?)
# Midway
# Species <- "BFAL"
# years <- c(2008,2009,2010,2011,2012)
# setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/cleaned/filtered_for_postbreeding/',Species,"/"))
# 
# file.list <- list.files(getwd())
# 
# big.df <- data.frame()
# for (i in 1:length(file.list)){
#   add.file <- read_csv(file.list[i])
#   add.file$year <- years[i] 
#   add.file$island <- "Midway"
#   add.file$spp <- Species
#   if (i==1){big.df <- add.file
#   } else {big.df <- rbind(big.df, add.file)}
# }
# 
# # write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/cleaned/filtered_for_postbreeding/allBFAL_SST_Midway.csv"),row.names = FALSE)
# 

# Intermediate step - merge 4 SST files -----------------------------------

# You now have 4 SST CSVs that were created from the daylogs from each species - Midway LAAL/BFAL created above here, 
# and the Tern LAAL/BFAL SST CSVs were created in a previous version of this script, but also pulled from temp logs in a long process during 
# the script sst_analysis_data_cleaning_tern.R. 
# Load all 4 in, then standardize and merge them into a master SST CSV. 

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/cleaned/all_LAAL_midway_SST.Rdata")
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/daylog_work/sst_data/cleaned/all_BFAL_midway_SST.Rdata")
all_LAAL_tern_SST <- read_csv("/users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/collated/allLAAL_SST_Tern.csv")
all_BFAL_tern_SST <- read_csv("/users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/melinda_data/sst_data/collated/allBFAL_SST_Tern.csv")

# standardize Midway: 
all_LAAL_midway_SST$island <- "Midway"
all_LAAL_midway_SST$spp <- "LAAL"
all_LAAL_midway_SST$time <- year(all_LAAL_midway_SST$date)
all_LAAL_midway_SST <- all_LAAL_midway_SST %>% rename(year=time)
all_LAAL_midway_SST <- all_LAAL_midway_SST[,c(1,2,4,3,5,6,7)]

all_BFAL_midway_SST$island <- "Midway"
all_BFAL_midway_SST$spp <- "BFAL"
all_BFAL_midway_SST$time <- year(all_BFAL_midway_SST$date)
all_BFAL_midway_SST <- all_BFAL_midway_SST %>% rename(year=time)
all_BFAL_midway_SST <- all_BFAL_midway_SST[,c(1,2,4,3,5,6,7)]

# standardize Tern, have to add Island Year Spp:

all_LAAL_tern_SST <- all_LAAL_tern_SST %>% rename(date=dtime)
all_BFAL_tern_SST <- all_BFAL_tern_SST %>% rename(date=dtime)
all_LAAL_tern_SST$date <- parse_date_time(all_LAAL_tern_SST$date, c("mdy HMS", "dmy HMS","mdy","dmy", "mdy HM", "mdY HMS", "dmY HMS", "dmY HM", "mdY HM"))
all_BFAL_tern_SST$date <- parse_date_time(all_BFAL_tern_SST$date, c("mdy HMS", "dmy HMS","mdy","dmy", "mdy HM", "mdY HMS", "dmY HMS", "dmY HM", "mdY HM"))
all_LAAL_tern_SST$year <- year(all_LAAL_tern_SST$date)
all_BFAL_tern_SST$year <- year(all_BFAL_tern_SST$date)
all_LAAL_tern_SST$spp <- "LAAL"
all_BFAL_tern_SST$spp <- "BFAL"
all_LAAL_tern_SST$island <- "Tern"
all_BFAL_tern_SST$island <- "Tern"
all_LAAL_tern_SST <- all_LAAL_tern_SST[,c(1,2,3,5,4,7,6)]
all_BFAL_tern_SST <- all_BFAL_tern_SST[,c(1,2,3,5,4,7,6)]

all_SST <- rbind(all_LAAL_midway_SST,all_LAAL_tern_SST,all_BFAL_midway_SST,all_BFAL_tern_SST)
all_SST <- na.omit(all_SST)

all_SST$month <- month(all_SST$date)
all_SST <- all_SST %>% mutate(month = month.name[month])

all_SST$date <- floor_date(all_SST$date, unit="day")

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/SST/")
#save(all_SST, file = "all_SST_master.Rdata")


# =======================================================================================

# START HERE: 
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/SST/")
load("all_SST_master.Rdata")
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
Island <- "Tern"
species <- "BFAL"

sst_loop <- all_SST %>% filter(spp==species)
sst_loop <- sst_loop %>% filter(island==Island)
track.data <- bt

# Cycle through each point, average matching SST measurements -----------------------------------------------------

for (i in 1:nrow(track.data)){
  date.temp <- track.data$dtime[i]
  sst.temp1 <- sst_loop[which(sst_loop$date==date.temp),]
  sst.temp2 <- sst.temp1[which(sst.temp1$id==track.data$id[i]),]
  if (nrow(sst.temp2)>=1){avg_loop <- mean(sst.temp2$SST)
  track.data$DailyMeanSST[i] <- avg_loop} else {track.data$DailyMeanSST[i]<-NA}
}

track.data<-na.omit(track.data)
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/SST/")

# save(track.data, file = "ternBFAL_DailyAverageSST.Rdata")



