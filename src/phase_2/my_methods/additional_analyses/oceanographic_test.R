# Script to do oceanographic analysis of points and create boxplots

# Tasks for today 
  # download 1 month of wind data, get that working
  # test with data for that month
  # create a boxplot of that month


# Testing from Tziporah



# Setup -------------------------------------------------------------------

library(tidyverse) # data cleanup
library(sf) # spatial
library(rgeos)
library(raster) # convert netcdf to raster
library(ncmeta)
library(ncdf4)
library(RNetCDF)
library(nngeo)
library(lubridate) # Clean up dates
library(rgdal)
library(stars) 
library(lattice)
library(rasterVis) # vectorplot
library(DescTools) #closest
library(imputeTS) #na.interpolation
library(swfscMisc)

###############################################################################
# The purpose of this code is to extract MERRA-2 wind data (u and v components 
# for 2m, 10 m, and 50m) for albatross GPS tracks  
###############################################################################

##########################################################################################
# AUXILIARY FUNCTIONS
# run before analyses to calculate wind speed 
##########################################################################################
# -----------------------------------------------------------------
# uv2ddff
# -----------------------------------------------------------------
uv2ddff <- function(u, v = NULL, rad = FALSE){
  # if input u is zoo or data.frame
  zoo_index <- NULL # Default
  if (inherits(u, c("zoo", "data.frame"))) {
    # If input 'u' is zoo: keep index
    if (inherits(u, "zoo")) zoo_index <- index(u)
    if (!all(c("u", "v") %in% names(u)))
      stop("necessary colums \"u\" and/or \"v\" missing")
    # If "v" is set in addition: warn
    if (!is.null(v)) {
      warning(sprintf("input \"u\" to uv2ddff is \"%s\":", class(u)),
              "\"v\" specified as well but will be ignored!")
    }
    v = as.numeric(u$v)
    u = as.numeric(u$u)
    # if u has 2 columns the second column is taken as v
  } else if (NCOL(u) == 2) {
    v <- u[,2]
    u <- u[,1]
  } else {
    if (is.null(v)) stop("input \"v\" missing")
    # If lenths do not match and none of them is of length 1: stop.
    if (!(length(u) == length(v)) & !any(c(length(u), length(v)) == 1L)) {
      stop("Length of \"u\" and \"v\" not identical")
      # Else recycle the one with length one to the length of the other one.
      # Damn it, so much one's in one sentence!
    } else if (length(u) == 1) {
      u <- rep(u, length(v))
    } else if (length(v) == 1) {
      v <- rep(v, length(u))
    }
  }
  # polar coordinates:
  ff <- sqrt(u^2 + v^2)
  dd <- atan(v/u) + (u < 0) * pi
  # Only non-na combis
  idx <- which(!is.na(dd) & !is.na(ff));   dd[idx] <- dd[idx] + 2 * pi
  # convert angle to meteorological convention
  dd <- 3 * pi / 2 - dd
  idx <- which(!is.na(dd) & !is.na(ff));   dd[idx] <- dd[idx] + 2 * pi
  # if rad (radiants) = F we have to convert to degrees.
  if (!rad) dd <- dd * 180 / pi
  res <- data.frame(dd, ff)
  if (is.null(zoo_index)) return(res) else return(zoo(res, zoo_index))
}

wrap360 = function(lon) {lon360<-ifelse(lon<0,lon+360,lon);return(lon360)}
Lon360to180 <- function(lon){
  ((lon + 180) %% 360) - 180
}

#################################################################################
# Download and Import data
#################################################################################
# extent tziporah used
# 80 N 80 S 120 W -120 E
# -110, -70, 0, -30 & 0, -70, 10, -30 Extent for Bird Island -110, -70, 10, -30)
# 165, 20, 180, 40 & -180, 20, -170, 40 Extent for Midway (40N, 20S, 165W, -170E) NOTE:  cannot cross anti-meridian, and raster was not correct from -180 to 180. Have to download two files for every date.
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
#  Useful Links
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
# File specification: https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
# How to download the data: https://daac.gsfc.nasa.gov/information/howto?title=How%20to%20Download%20MERRA-2%20Daily%20Mean%20Data
# How to calculate and plot wind speed using MERRA-2 wind component: https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20calculate%20and%20plot%20wind%20speed%20using%20MERRA-2%20wind%20component%20data%20using%20Python
# Hourly Data: https://disc.gsfc.nasa.gov/datasets/M2I1NXASM_5.12.4/summary

nc_dir<-"E:/project_data/spatial_segregation/data/oceanographic"
setwd(nc_dir)

chla_files <- list.files(pattern='*.nc') 

for (i in 1:length(chla_files)) {
  mi<-read_ncdf(chla_files[i])
  times <- st_get_dimension_values(mi, "time")
  if (i==1) {
    times_all<-times
  }else{
    times_all<-c(times_all,times)
  }
}

 # have to download in two parts then merge
# download 120 to 180 and -180 to -120, and for latitude both 0 to 75

chla <- read_ncdf(a)
times <- st_get_dimension_values(mi, "time")
all_times<-as.tibble(times)
all_times.name<-as.character(all_times)
all_times_num<-as.numeric(unlist(all_times))

# isolate chla 
chla.df.left <- raster::as.data.frame(mi, xy=TRUE)




# setwd("/Users/tziporahserota/Desktop/Midway2/")
# b <- list.files("/Users/tziporahserota/Desktop/Midway2/", pattern = ".nc")
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
#  Create list of times
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
for (i in 1:length(a)) {
  mi<-read_ncdf(a[i])
  times <- st_get_dimension_values(mi, "time")
  if (i==1) {
    times_all<-times
  }else{
    times_all<-c(times_all,times)
  }
}

all_times<-as.tibble(times_all)
all_times.name<-as.character(all_times)
all_times_num<-as.numeric(unlist(all_times))

# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
#  Isolate U and V components for 2 and 10 m 
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
# U-component, 2 m 
for (i in 1:length(a)) {
  mi<-read_ncdf(a[i])
  wind_u<-as(mi[2,,,], "Raster")
  wind_u.df= raster::as.data.frame(wind_u, xy = TRUE)
  if (i==1) {
    wind_U2M<-wind_u.df
  }else{
    wind_U2M<-cbind(wind_U2M,wind_u.df[, -c(1,2)])
  }
}

# V-component, 2 m 
for (i in 1:length(a)) {
  mi<-read_ncdf(a[i])
  wind_v<-as(mi[5,,,], "Raster")
  wind_v.df= raster::as.data.frame(wind_v, xy = TRUE)
  if (i==1) {
    wind_V2M<-wind_v.df
  }else{
    wind_V2M<-cbind(wind_V2M,wind_v.df[, -c(1,2)])
  }
}


# U-component, 10 m 
for (i in 1:length(a)) {
  mi<-read_ncdf(a[i])
  wind_u<-as(mi[1,,,], "Raster")
  wind_u.df= raster::as.data.frame(wind_u, xy = TRUE)
  if (i==1) {
    wind_U10M<-wind_u.df
  }else{
    wind_U10M<-cbind(wind_U10M,wind_u.df[, -c(1,2)])
  }
}

# V-component, 10 m 
for (i in 1:length(a)) {
  mi<-read_ncdf(a[i])
  wind_v<-as(mi[4,,,], "Raster")
  wind_v.df= raster::as.data.frame(wind_v, xy = TRUE)
  if (i==1) {
    wind_V10M<-wind_v.df
  }else{
    wind_V10M<-cbind(wind_V10M,wind_v.df[, -c(1,2)])
  }
}


# U-component, 50 m 
for (i in 1:length(a)) {
  mi<-read_ncdf(a[i])
  wind_u<-as(mi[3,,,], "Raster")
  wind_u.df= raster::as.data.frame(wind_u, xy = TRUE)
  if (i==1) {
    wind_U50M<-wind_u.df
  }else{
    wind_U50M<-cbind(wind_U50M,wind_u.df[, -c(1,2)])
  }
}

# V-component, 50 m 
for (i in 1:length(a)) {
  mi<-read_ncdf(a[i])
  wind_v<-as(mi[6,,,], "Raster")
  wind_v.df= raster::as.data.frame(wind_v, xy = TRUE)
  if (i==1) {
    wind_V50M<-wind_v.df
  }else{
    wind_V50M<-cbind(wind_V50M,wind_v.df[, -c(1,2)])
  }
}

# have to repeat with second data file
setwd("/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/Merra-2Datasets/Hourly/Dec2019-Feb2020_WAALinc/R/")
b <- list.files("/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/Merra-2Datasets/Hourly/Dec2019-Feb2020_WAALinc/R/", pattern = ".nc")

# U-component, 2 m 
for (i in 1:length(b)) {
  mi<-read_ncdf(b[i])
  wind_u2<-as(mi[2,,,], "Raster")
  wind_u.df2= raster::as.data.frame(wind_u2, xy = TRUE)
  if (i==1) {
    wind_U2M2<-wind_u.df2
  }else{
    wind_U2M2<-cbind(wind_U2M2,wind_u.df2[, -c(1,2)])
  }
}

# V-component, 2 m 
for (i in 1:length(b)) {
  mi<-read_ncdf(b[i]) 
  wind_v2<-as(mi[5,,,], "Raster")
  wind_v.df2= raster::as.data.frame(wind_v2, xy = TRUE)
  if (i==1) {
    wind_V2M2<-wind_v.df2
  }else{
    wind_V2M2<-cbind(wind_V2M2,wind_v.df2[, -c(1,2)])
  }
}


# U-component, 10 m 
for (i in 1:length(b)) {
  mi<-read_ncdf(b[i])
  wind_u2<-as(mi[1,,,], "Raster")
  wind_u.df2= raster::as.data.frame(wind_u2, xy = TRUE)
  if (i==1) {
    wind_U10M2<-wind_u.df2
  }else{
    wind_U10M2<-cbind(wind_U10M2,wind_u.df2[, -c(1,2)])
  }
}

# V-component, 10 m 
for (i in 1:length(b)) {
  mi<-read_ncdf(b[i])
  wind_v2<-as(mi[4,,,], "Raster")
  wind_v.df2= raster::as.data.frame(wind_v2, xy = TRUE)
  if (i==1) {
    wind_V10M2<-wind_v.df2
  }else{
    wind_V10M2<-cbind(wind_V10M2,wind_v.df2[, -c(1,2)])
  }
}


# U-component, 50 m 
for (i in 1:length(b)) {
  mi<-read_ncdf(b[i])
  wind_u2<-as(mi[3,,,], "Raster")
  wind_u.df2= raster::as.data.frame(wind_u2, xy = TRUE)
  if (i==1) {
    wind_U50M2<-wind_u.df2
  }else{
    wind_U50M2<-cbind(wind_U50M2,wind_u.df2[, -c(1,2)])
  }
}

# V-component, 50 m 
for (i in 1:length(b)) {
  mi<-read_ncdf(b[i])
  wind_v2<-as(mi[6,,,], "Raster")
  wind_v.df2= raster::as.data.frame(wind_v2, xy = TRUE)
  if (i==1) {
    wind_V50M2<-wind_v.df2
  }else{
    wind_V50M2<-cbind(wind_V50M2,wind_v.df2[, -c(1,2)])
  }
}

wind_U2M <- rbind(wind_U2M, wind_U2M2)
wind_U10M <- rbind(wind_U10M, wind_U10M2)
wind_U50M <- rbind(wind_U50M, wind_U50M2)
wind_V2M <- rbind(wind_V2M, wind_V2M2)
wind_V10M <- rbind(wind_V10M, wind_V10M2)
wind_V50M <- rbind(wind_V50M, wind_V50M2)

# convert to rasters
u_stack2M<-rasterFromXYZ(wind_U2M)
v_stack2M<-rasterFromXYZ(wind_V2M)
u_stack10M<-rasterFromXYZ(wind_U10M)
v_stack10M<-rasterFromXYZ(wind_V10M)
u_stack50M<-rasterFromXYZ(wind_U50M)
v_stack50M<-rasterFromXYZ(wind_V50M)




























# Read chla in first
setwd("E:/project_data/spatial_segregation/data/oceanographic")

chla_files <- list.files(pattern='*.nc') 

chla <- read_ncdf(chla_files[2])
chla # Make sure you got the right stuff!

times <- st_get_dimension_values(chla, "time")
chla_value <- as(chla[1,,,], "Raster")

# Stack all wind datasets.
# Each stack will be associated with a date and hourly timestamp (24 per day)
nc_dir<-'C:/Users/dalla/Downloads/'
setwd(nc_dir)

wind_files <- list.files(pattern='*.nc') 

wind_t1 <- read_ncdf(wind_files[3])
wind_t1 # Make sure you got the right stuff!
times_t1<-st_get_dimension_values(wind_t1, "time")  # stores times from each file 
wind_t1_u<-as(wind_t1[1,,,], "Raster") # refer to U component
wind_t1_v<-as(wind_t1[2,,,], "Raster") # v component 


wind_t2 <- read_ncdf(wind_files[2])
wind_t2
times_t2<-st_get_dimension_values(wind_t2, "time")  
wind_t2_u<-as(wind_t2[1,,,], "Raster")
wind_t2_v<-as(wind_t2[2,,,], "Raster")


# Stack rasters  ....................
u_stack<-stack(wind_t1_u, wind_t2_u)
v_stack<-stack(wind_t1_v, wind_t2_v)

all_times<-rbind(as.tibble(times_t1), as.tibble(times_t2))


# -----------------------------------------------------------------
# uv2ddff
# -----------------------------------------------------------------
uv2ddff <- function(u, v = NULL, rad = FALSE){
  # if input u is zoo or data.frame
  zoo_index <- NULL # Default
  if (inherits(u, c("zoo", "data.frame"))) {
    # If input 'u' is zoo: keep index
    if (inherits(u, "zoo")) zoo_index <- index(u)
    if (!all(c("u", "v") %in% names(u)))
      stop("necessary colums \"u\" and/or \"v\" missing")
    # If "v" is set in addition: warn
    if (!is.null(v)) {
      warning(sprintf("input \"u\" to uv2ddff is \"%s\":", class(u)),
              "\"v\" specified as well but will be ignored!")
    }
    v = as.numeric(u$v)
    u = as.numeric(u$u)
    # if u has 2 columns the second column is taken as v
  } else if (NCOL(u) == 2) {
    v <- u[,2]
    u <- u[,1]
  } else {
    if (is.null(v)) stop("input \"v\" missing")
    # If lenths do not match and none of them is of length 1: stop.
    if (!(length(u) == length(v)) & !any(c(length(u), length(v)) == 1L)) {
      stop("Length of \"u\" and \"v\" not identical")
      # Else recycle the one with length one to the length of the other one.
      # Damn it, so much one's in one sentence!
    } else if (length(u) == 1) {
      u <- rep(u, length(v))
    } else if (length(v) == 1) {
      v <- rep(v, length(u))
    }
  }
  # polar coordinates:
  ff <- sqrt(u^2 + v^2)
  dd <- atan(v/u) + (u < 0) * pi
  # Only non-na combis
  idx <- which(!is.na(dd) & !is.na(ff));   dd[idx] <- dd[idx] + 2 * pi
  # convert angle to meteorological convention
  dd <- 3 * pi / 2 - dd
  idx <- which(!is.na(dd) & !is.na(ff));   dd[idx] <- dd[idx] + 2 * pi
  # if rad (radiants) = F we have to convert to degrees.
  if (!rad) dd <- dd * 180 / pi
  res <- data.frame(dd, ff)
  if (is.null(zoo_index)) return(res) else return(zoo(res, zoo_index))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # testing with melinda's code
  # Setup -------------------------------------------------------------------
  
  require(sf)
  require(tidyverse)
  require(raster)
  library(nngeo)
  library(lubridate)
  library(rgdal)
  library(stars)
  library(ncmeta)
  library(ncdf4)
  library(RNetCDF)
  library(lattice)
  library(rCAT) #rad2deg()
  library(rasterVis) #vectorplot
  library(colorRamps) #matlab.like
  library(viridisLite)
  library(colorspace)
  library(DescTools) #closest
  library(imputeTS) #na.interpolation
  library(swfscMisc)
  require(ggplot2)
  require(RColorBrewer)
  
  wrap360 = function(lon) {lon360<-ifelse(lon<0,lon+360,lon);return(lon360)}
  Lon360to180 <- function(lon){
    ((lon + 180) %% 360) - 180
  }
  
  
  # Note: Important: Run code at bottom of script for uv2ddff function
  
  ########################################################################
  # 1. Precursor to Analysis
  # Obtain wind dataset from ECMWF ERA5 at sea level (1000 hPa)
  # Grid Extent: -30N -70S -80E -20E
  # Here: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form
  ########################################################################
  
  ########################################################################
  # 1. Import Wind Data (Downloaded as NetCDF grids from ECMWF)
  ########################################################################
  # Stack all wind datasets.
  # Each stack will be associated with a date and hourly timestamp (24 per day)
  
  nc_dir<-"E:/project_data/spatial_segregation/data/oceanographic"
  setwd(nc_dir)
  
  wind_files <- list.files(pattern='*.nc') 
  
  wind_t1 <- read_ncdf(wind_files[1])
  wind_t1 # Make sure you got the right stuff!
  times_t1<-st_get_dimension_values(wind_t1, "time")  # stores times from each file 
  wind_t1_u<-as(wind_t1[1,,,], "Raster") # refer to U component
  wind_t1_v<-as(wind_t1[2,,,], "Raster") # v component 
  
  
  wind_t2 <- read_ncdf(wind_files[2])
  wind_t2
  times_t2<-st_get_dimension_values(wind_t2, "time")  
  wind_t2_u<-as(wind_t2[1,,,], "Raster")
  wind_t2_v<-as(wind_t2[2,,,], "Raster")
  
  
  # Stack rasters  ....................
  u_stack<-stack(wind_t1_u, wind_t2_u)
  v_stack<-stack(wind_t1_v, wind_t2_v)
  
  all_times<-rbind(as.tibble(times_t1), as.tibble(times_t2))
  all_times_num<-as.numeric(unlist(all_times)) # append times: one big columns of time
  
}