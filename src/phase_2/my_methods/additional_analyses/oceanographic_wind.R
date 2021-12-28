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
library(dplyr)

###############################################################################
# The purpose of this code is to extract environmental data (wind u and v components 
# for ___ m resolution) for albatross GLS tracks  
###############################################################################


# AUXILIARY FUNCTIONS - run before analyses to calculate wind speed --------

# uv2ddff

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



# Download and import data ------------------------------------------------
# extent tziporah used
# 80 N 80 S 120 W -120 E
# -110, -70, 0, -30 & 0, -70, 10, -30 Extent for Bird Island -110, -70, 10, -30)
# 165, 20, 180, 40 & -180, 20, -170, 40 Extent for Midway (40N, 20S, 165W, -170E) NOTE:  cannot cross anti-meridian, and raster was not correct from -180 to 180. Have to download two files for every date.

#  Useful Links

# File specification: https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
# How to download the data: https://daac.gsfc.nasa.gov/information/howto?title=How%20to%20Download%20MERRA-2%20Daily%20Mean%20Data
# How to calculate and plot wind speed using MERRA-2 wind component: https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20calculate%20and%20plot%20wind%20speed%20using%20MERRA-2%20wind%20component%20data%20using%20Python
# Hourly Data: https://disc.gsfc.nasa.gov/datasets/M2I1NXASM_5.12.4/summary

nc_dir<-"E:/project_data/spatial_segregation/data/oceanographic"
setwd(nc_dir)

# for mac
# setwd("/Volumes/Samsung_T5/project_data/spatial_segregation/data/oceanographic")
setwd("/Users/dallasjordan/Desktop/ocean_data/2008/jan2/left")
wind_files_left <- list.files(pattern='*.nc') 
setwd("/Users/dallasjordan/Desktop/ocean_data/2008/jan2/right")
wind_files_right <- list.files(pattern='*.nc')

# list of dates that you have
#can just read in left to get times, since they are the same times in both files
for (i in 1:length(wind_files_left)) {
  mi<-read_ncdf(wind_files_left[i])
  times <- st_get_dimension_values(mi, "time")
  if (i==1) {
    times_all<-times
  }else{
    times_all<-c(times_all,times)
  }
}

all_times<-as_tibble(times_all)
all_times.name<-as.character(all_times)
all_times_num<-as.numeric(unlist(all_times))

# have to download in two parts then merge
# download 120 to 180 and -180 to -120, and for latitude both 0 to 75


# Isolate wind for left of antimeridian -----------------------------------

# U-component, 2 m 
for (i in 1:length(wind_files_left)) {
  mi<-read_ncdf(wind_files_left[i])
  wind_u<-as(mi[1,,,], "Raster")
  wind_u.df= raster::as.data.frame(wind_u, xy = TRUE)
  if (i==1) {
    wind_U2M<-wind_u.df
  }else{
    wind_U2M<-cbind(wind_U2M,wind_u.df[, -c(1,2)])
  }
}

# V-component, 2 m 
for (i in 1:length(wind_files_left)) {
  mi<-read_ncdf(wind_files_left[i])
  wind_v<-as(mi[3,,,], "Raster")
  wind_v.df= raster::as.data.frame(wind_v, xy = TRUE)
  if (i==1) {
    wind_V2M<-wind_v.df
  }else{
    wind_V2M<-cbind(wind_V2M,wind_v.df[, -c(1,2)])
  }
}


# U-component, 10 m 
for (i in 1:length(wind_files_left)) {
  mi<-read_ncdf(wind_files_left[i])
  wind_u<-as(mi[6,,,], "Raster")
  wind_u.df= raster::as.data.frame(wind_u, xy = TRUE)
  if (i==1) {
    wind_U10M<-wind_u.df
  }else{
    wind_U10M<-cbind(wind_U10M,wind_u.df[, -c(1,2)])
  }
}

# V-component, 10 m 
for (i in 1:length(wind_files_left)) {
  mi<-read_ncdf(wind_files_left[i])
  wind_v<-as(mi[4,,,], "Raster")
  wind_v.df= raster::as.data.frame(wind_v, xy = TRUE)
  if (i==1) {
    wind_V10M<-wind_v.df
  }else{
    wind_V10M<-cbind(wind_V10M,wind_v.df[, -c(1,2)])
  }
}


# U-component, 50 m 
for (i in 1:length(wind_files_left)) {
  mi<-read_ncdf(wind_files_left[i])
  wind_u<-as(mi[5,,,], "Raster")
  wind_u.df= raster::as.data.frame(wind_u, xy = TRUE)
  if (i==1) {
    wind_U50M<-wind_u.df
  }else{
    wind_U50M<-cbind(wind_U50M,wind_u.df[, -c(1,2)])
  }
}

# V-component, 50 m 
for (i in 1:length(wind_files_left)) {
  mi<-read_ncdf(wind_files_left[i])
  wind_v<-as(mi[2,,,], "Raster")
  wind_v.df= raster::as.data.frame(wind_v, xy = TRUE)
  if (i==1) {
    wind_V50M<-wind_v.df
  }else{
    wind_V50M<-cbind(wind_V50M,wind_v.df[, -c(1,2)])
  }
}








u_stack2M<-rasterFromXYZ(wind_U2M)

# Isolate wind for right of antimeridian ----------------------------------

# wind right
chla_files_right <- chla_files[grep("right",chla_files)]
for (i in 1:length(chla_files_right)) {
  mi<-read_ncdf(chla_files_right[i])
  chla_right.df= raster::as.data.frame(mi, xy = TRUE)
  if (i==1) {
    chla_right<-chla_right.df
  }else{
    chla_right<-rbind(chla_right,chla_right.df)
  }
}
# combine left and right
chla_all <- rbind(chla_left, chla_right)

# convert to rasters
chla_stack<-rasterFromXYZ(chla_left)

# You have a df of everything - so now get your points, create a buffer around each point, then extract to that polygon


# Load in GLS points ------------------------------------------------------

# for (i in 1:length(files)) {
#   bi<-read.csv(files[i])
#   if (i==1) {
#     b<-bi
#   }else{
#     b<-rbind(b,bi)
#   }
# } 

# Start with Midway LAAL 2008

points <- read_csv("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/LAAL/2008/pb_1108_01.csv")

# read in netcdf - pull out a specific date according to what you need

# chla_files_left <- chla_files[grep("left",chla_files)]
# test_mi <- NA
# for (i in 1:length(chla_files_left)) {
#   test_mi[i] = brick(chla_files_left[i], xy = TRUE)
# }
# 
# chla_all
# temp <- chla_all[which(chla_all$date== date.temp), ]
# 
# test_1 <- read_ncdf(chla_files_left[1])
# chla.df = data.table(as.data.frame(test_1, xy = TRUE))
# 
# # rasterize based by month
# test_January <- chla.df %>% subset(month(chla.df$time)==1)
# raster_January <- rasterFromXYZ(jan)
# 
# mi<-read_ncdf(a[i])
# 
# 
# test_2 <- brick(test_1)
# chla_raster <- brick(chla_files_left[1], varname="chlorophyll")
# chla_raster; class(tmp_raster)
# 
# raster

# load in January data set

# jan_chla_km4 <- read_csv("~/Downloads/chla_jan.csv", skip=1)
# jan_chla_km4<-cbind(km4.data[1,],0,0,0,0,0,0,0)
# july_c <- read_csv("~/Downloads/july_c.csv", skip=1)
# july_c <- july_c[,c(1,3,2,4)]
# 
# 
# colnames(july_c) <- c('UTC', 'y', 'x', 'value')
# e <- extent(july_c[,2:3])
# r <- raster(e, ncol=8640, nrow=4320)
# x <- rasterize(july_c[, 2:3], r, july_c[,4], fun=mean)
# 
# test_1 <- rasterFromXYZ(july_c)

# read netCDF in 
# convert netCDF to raster


projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
df <- st_as_sf(x = points,                         
               coords = c("Longitude", "Latitude"),crs=projcrs)

test <- read_stars("~/Downloads/jul2008chla.nc")

test<- t(test)
plot(test)
test= raster::as.data.frame(test, xy = TRUE)
diff_long<-diff(test$x)
diff_long<- as.data.frame(diff_long)
diff_long <- rbind(0,diff_long)
test <- cbind(test,diff_long)

counts<- test %>% count(diff_long)
counts

test<- read_csv("~/Downloads/test_csv.csv")
test<-test[-1,]
test$latitude<-as.numeric(test$latitude)
test$longitude<-as.numeric(test$longitude)
test$latitude<-round(test$latitude,digits=2)
test$longitude<-round(test$longitude,digits=2)
diff_lat <- diff(test$latitude)
diff_lat <- as.data.frame(diff_lat)
diff_lat <- rbind(0,diff_lat)
test <- cbind(test,diff_lat)
unique(test$diff_lat)

rast <- rasterFromXYZ(test)

test<- test[,-5]
colnames(test) <- c('UTC', 'y', 'x', 'value')
e <- extent(test[,2:3])
r <- raster(e, ncol=8640, nrow=4320)
test <- test[,c("UTC",'x','y','value')]
x <- rasterize(test[, 2:3], r, test[,4], fun=mean)




test_rast <- rasterFromXYZ(test)
st_crs(test)<- projcrs

poly = st_buffer(df, dist = 3)
x = aggregate(test, poly, mean)
st_as_sf(x)

extracted <- st_extract(test,df, time_column = "dtime", interpolate_time = FALSE)
# This any day before the 16th of the month gets rounded to the previous month



# year<-2008
# month<-7
# if (month(df$dtime)==7){
#   df$dtime[month(df$dtime==7),]<-ymd(20080716)
# }
# df$dtime[month(df$dtime==7),]<-ymd(20080716)






# this is overflowing. so, need to chop it up into smaller rasters. Subset the dataframe by lat lon parameters and make rasters out of those

chla_small_left <- chla_left %>% filter(latitude>70)

# chla_small_left worked. Took a long time to process but filtered to >70. 

radius <- 150000 # in meters, 150,000 m is 150km radius, makes 300km x 300km box

library(geosphere)

# distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine) to calc meter distance between two points at any lat

makeBox <- fun(lat,lon){
  yPlus <- lat+radius # this is a degree plus a meter - fix this
  xPlus <- lon+radius
  yMinus <- lat-radius
  xMinus <- lon-radius
}
# make smaller dataframes. make those rasters. Divide up points into necessary rasters. 































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