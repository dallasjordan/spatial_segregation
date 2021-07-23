# Sample size analysis
# Dallas Jordan, Feb 1 2021
# This script takes my modeled postbreeding tracks for LAAL and BFAL from Midway and Tern, conducts KDE using adehabitatHR, and calculates the area for each 95% UD.
# The first half of this script is a modified version of my script 'calculate_UD_area_midway' and 'calculate_UD_area_tern'. Modifications include saving the areas calculated
# to an ID (year and colony) in a new data frame. This data frame for each island is then run through the bootstrapping procedures described by Lascelles et al. 2016. A nonlinear
# regression is run on each plot to determine the asymptote (.... finish this)

library(spatstat)
library(dplyr)
library(maptools)
library(rgdal)
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)
library(geosphere)
library(fields)

###############################################################################

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
years <- c("2008","2009","2010","2011","2012")
pdc_mercator_proj<-sf::st_crs(3832) # used towards end of script

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  loaded <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(loaded)){
    loaded[[j]]$ind_id <- paste0(j)
  }
  data_loop <- do.call(rbind, loaded)
  colnames(data_loop)= c("dtime","x","y", "ind_id")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","x","y","ind_id")]


# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  loaded <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(loaded)){
    loaded[[j]]$ind_id <- paste0(j)
  }
  data_loop <- do.call(rbind, loaded)
  colnames(data_loop)= c("dtime","x","y", "ind_id")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","x","y", "ind_id")]

# Potentially useful for later: 
          # combine into 1 dataframe, preserve id
          
          data.sp <- rbind(data1.sp,data2.sp)
          
          data.spL <- split(data.sp, data.sp$id)
          
          resample_2008 <- rbind(data.spL[["LAAL_2008"]],data.spL[["BFAL_2008"]])
          resample_2009 <- rbind(data.spL[["LAAL_2009"]],data.spL[["BFAL_2009"]])
          resample_2010 <- rbind(data.spL[["LAAL_2010"]],data.spL[["BFAL_2010"]])
          resample_2011 <- rbind(data.spL[["LAAL_2011"]],data.spL[["BFAL_2011"]])
          resample_2012 <- rbind(data.spL[["LAAL_2012"]],data.spL[["BFAL_2012"]])
          
          # Replace stuff
          
          resample_list_2008 <- split(resample_2008, resample_2008$id)
          current_LAAL <- resample_list_2008[["LAAL_2008"]]
          current_BFAL <- resample_list_2008[["BFAL_2008"]]

# Sample based on ID numbers - then take that sample, calculate an area, and save it. Then, sample 2 ID numbers - take those samples, calculate an area, and save it, etc. 


current_LAAL <- data1.sp
coordinates(current_LAAL) <- c("x","y")
proj4string(current_LAAL) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
current_LAAL <- spTransform(current_LAAL,pdc_mercator_proj$proj4string)

coordinates(current_BFAL) <- c("x","y")
proj4string(current_BFAL) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
current_BFAL<- spTransform(current_BFAL,pdc_mercator_proj$proj4string)

kernel.ref <- kernelUD(current_LAAL, same4all = T, grid=30)
data.kernel.poly <- getverticeshr(kernel.ref, percent = 90, unin = "m", unout = "km2") 
data.kernel.poly$id
data.kernel.poly$area

kernel.ref <- kernelUD(current_BFAL, same4all = T, grid=30)
data.kernel.poly <- getverticeshr(kernel.ref, percent = 90, unin = "m", unout = "km2") 
data.kernel.poly$id
data.kernel.poly$area



# Question of "which years we need to throw out" appears to be an issue - the sample sizes aren't large enough to hit an asymptote. 
# First one I want to do is all LAAL from Midway. the bootstrap is used to do the nonlinear regression. 

