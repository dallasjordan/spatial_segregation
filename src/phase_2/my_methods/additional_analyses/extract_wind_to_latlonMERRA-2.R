#################################################################################
# Load Libraries 
#################################################################################
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
# -110, -70, 0, -30 & 0, -70, 10, -30 Extent for Bird Island -110, -70, 10, -30)
# 165, 20, 180, 40 & -180, 20, -170, 40 Extent for Midway (40N, 20S, 165W, -170E) NOTE:  cannot cross anti-meridian, and raster was not correct from -180 to 180. Have to download two files for every date.
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
#  Useful Links
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
# File specification: https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
# How to download the data: https://daac.gsfc.nasa.gov/information/howto?title=How%20to%20Download%20MERRA-2%20Daily%20Mean%20Data
# How to calculate and plot wind speed using MERRA-2 wind component: https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20calculate%20and%20plot%20wind%20speed%20using%20MERRA-2%20wind%20component%20data%20using%20Python
# Hourly Data: https://disc.gsfc.nasa.gov/datasets/M2I1NXASM_5.12.4/summary

setwd("/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/Merra-2Datasets/Hourly/Dec2019-Feb2020_WAALinc/L/")
a <- list.files("/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/Merra-2Datasets/Hourly/Dec2019-Feb2020_WAALinc/L/", pattern = ".nc")

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

####################################################################################################
# 2. Read in hourly bird locations and gather wind data
####################################################################################################
int_now <- 3600
setwd(paste0('/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/WAAL_L1_interpolated/WAAL_2019-2020_incubation/', int_now, 's/'))
files<-list.files(pattern = ".csv")

# Append all individual bird lat lon files to one file ------
for (i in 1:length(files)) {
  bi<-read.csv(files[i])
  if (i==1) {
    b<-bi
  }else{
    b<-rbind(b,bi)
  }
} 

# Add speed, distance, bearing columns to m -----------------
b$ground_speed_kmHr <- NA
b$distanceij_km <- NA
b$bearingij <- NA
b$U2M<- NA
b$V2M<- NA
b$U10M<- NA
b$V10M<- NA
b$U50M<- NA
b$V50M<- NA
b$wind_speed2M <- NA
b$wind_speed10M <- NA
b$wind_speed50M <- NA


trips<- unique(b$tripID)
for (i in 1:length(trips)) {
  tripi<-b[b$tripID==trips[i],]
  
  for (j in 1:length(tripi$id)-1) {
    hour_int<- int_now/3600 # int_now is in seconds (3600 seconds in one hour)
    tripi$distanceij_km[j]<-swfscMisc::distance(tripi$lat[j],tripi$lon[j],tripi$lat[j+1],tripi$lon[j+1], units = "km", method = "haversine")
    tripi$ground_speed_kmHr[j]<-tripi$distanceij_km[j]/hour_int    
    tripi$bearingij[j]<-as.numeric(bearing(tripi$lat[j],tripi$lon[j],tripi$lat[j+1],tripi$lon[j+1])[1])
  }
  
  b[b$tripID==trips[i],]<-tripi
  
}

# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
# Loop through m and add wind information: u, v, velocity for 2M and 10M
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
for ( j in 1:length(b$id)) {
  
  timej <- as.POSIXct(b$datetime[j], format = "%Y-%m-%d %H:%M:%S" , tz = "UTC")
  timej_num <- as.numeric(timej)
  # Find index of current_time in all times. Use that index to pull out relevant raster layer. 
  raster_dt_index <- as.numeric(which(abs(all_times_num-timej_num) == min(abs(all_times_num-timej_num))))
  
  # Isolate u and v rasters at time j
  ustack_timej2 <- subset(u_stack2M, raster_dt_index, drop=TRUE)
  vstack_timej2 <- subset(v_stack2M, raster_dt_index, drop=TRUE)
  ustack_timej10 <- subset(u_stack10M, raster_dt_index, drop=TRUE)
  vstack_timej10 <- subset(v_stack10M, raster_dt_index, drop=TRUE)
  ustack_timej50 <- subset(u_stack50M, raster_dt_index, drop=TRUE)
  vstack_timej50 <- subset(v_stack50M, raster_dt_index, drop=TRUE)
  
  # isolate coordinates
  xy_j <- as.data.frame(cbind(b$lon[j], b$lat[j]))
  colnames(xy_j) <- c("lon","lat")
  xy_j$lon <- Lon360to180(xy_j$lon) # not necessary for midway
  
  # Extract u and v components for time j at location x and y
  u_j2 <- extract(ustack_timej2, xy_j)
  v_j2 <- extract(vstack_timej2, xy_j)
  b$U2M[j]<- u_j2
  b$V2M[j]<- v_j2
  
  u_j10 <- extract(ustack_timej10, xy_j)
  v_j10 <- extract(vstack_timej10, xy_j)
  b$U10M[j]<- u_j10
  b$V10M[j]<- v_j10
  
  u_j50 <- extract(ustack_timej50, xy_j)
  v_j50 <- extract(vstack_timej50, xy_j)
  b$U50M[j]<- u_j50
  b$V50M[j]<- v_j50
  # -----------------------------------------------------------------------------
  # Get Wind Direction and Wind Velocity from U and V Components:
  # -----------------------------------------------------------------------------
  # Requires mathematical to meteorological adjustment:
  # Important: our wind vectors were given in mathematical notation, not in meteorological convention. 
  # Need to adjust by 270* - this is what uv2ddff does.
  # http://colaweb.gmu.edu/dev/clim301/lectures/wind/wind-uv
  # from uv2ddff: Wind direction is either in meteorological degrees (0 from North, from 90 East,
  # 180 from South, and 270 from West) or in mathematical radiant if input \code{rad = TRUE}.
  
  res2M<-uv2ddff(u_j2, v_j2)
  b$wind_speed2M[j] <- res2M$ff
  
  res10M<-uv2ddff(u_j10, v_j10)
  b$wind_speed10M[j] <- res10M$ff
  
  res50M<-uv2ddff(u_j50, v_j50)
  b$wind_speed50M[j] <- res50M$ff
}

spvec <- substr(b$id,1,4)
b<-b %>% mutate(species=spvec)

write_csv(b, '/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/wind_paired-data/MERRA-2/L0_paired/2019-2020_WAAL_incubation/hourly/allbirds_hourly.csv')


####################################################################################################
# 3. Interpolate wind data for 30 sec data
####################################################################################################

# Import 30s lat lon data and append individual bird files into one file. 
int_now <- 30
setwd(paste0('/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/WAAL_L1_interpolated/WAAL_2019-2020_incubation/', int_now, 's/'))
files<-list.files(pattern=".csv")

# Append all individual bird lat lon files to one file ------
for (i in 1:length(files)) {
  mi<-read.csv(files[i]) 
  if (i==1) {
    m_hi<-mi
  }else{
    m_hi<-rbind(m_hi,mi)
  }
} 

# placeholders 
m_hi$distanceij_km <- NA
m_hi$ground_speed_kmHr <- NA
m_hi$bearingij <- NA
m_hi$U2M<- NA
m_hi$V2M<- NA
m_hi$U10M<- NA
m_hi$V10M<- NA
m_hi$U50M<- NA
m_hi$V50M<- NA
m_hi$wind_speed2M <- NA
m_hi$wind_speed10M <- NA
m_hi$wind_speed50M <- NA
m_hi$wind_dir3602M <- NA
m_hi$wind_dir36010M <- NA
m_hi$wind_dir36050M <- NA
m_hi$windshear10m <- NA
m_hi$windshear50m <- NA
m_hi$bwa2M<- NA
m_hi$bwa_class2M<- NA
m_hi$bwa10M<- NA
m_hi$bwa_class10M<- NA

# -----------------------------------------------------------------------------
# Add Bird Speed, Distance, Bearing
# -----------------------------------------------------------------------------

trips<- unique(m_hi$tripID)
for (i in 1:length(trips)) {
  tripi<-m_hi[m_hi$tripID==trips[i],]
  
  for (j in 1:length(tripi$id)-1) {
    hour_int<- int_now/3600 # int_now is in seconds (3600 seconds in one hour)
    tripi$distanceij_km[j]<-swfscMisc::distance(tripi$lat[j],tripi$lon[j],tripi$lat[j+1],tripi$lon[j+1], units = "km", method = "haversine")
    tripi$ground_speed_kmHr[j]<-tripi$distanceij_km[j]/hour_int
    tripi$bearingij[j]<-as.numeric(bearing(tripi$lat[j],tripi$lon[j],tripi$lat[j+1],tripi$lon[j+1])[1])
  }
  
  m_hi[m_hi$tripID==trips[i],]<-tripi
  
}

# -----------------------------------------------------------------------------
# Add Wind Metrics
# -----------------------------------------------------------------------------

# Import hourly lat lon with wind u and v
m_lo <- b
# m_lo <- read.csv('/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/wind_paired-data/L0_paired/2019-2020/hourly/allbirds_hourly.csv')
birds<-unique(m_lo$id)



for (i in 1:length(birds)) {
  birdi<-birds[i]
  
  mi_lo = m_lo %>% filter(id == birds[i])
  mi_hi = m_hi %>% filter(id == birds[i])
  
  # Match times in datasets and extract u and v into hi res dataset every hour from hourly dataset
  hi_res_timematch<-as.numeric(as.POSIXct(as.character(mi_hi$datetime)))
  
  for (j in 1:length(mi_lo$id)) {
    j_match<-as.numeric(as.POSIXct(as.character(mi_lo$datetime[j])))
    match_ix<-Closest( hi_res_timematch, j_match, which=TRUE) # Which in m_hi is nearest j_match (m$datetime[j])
    mi_hi$U2M[match_ix] <- mi_lo$U2M[j]
    mi_hi$V2M[match_ix] <- mi_lo$V2M[j]
    mi_hi$U10M[match_ix] <- mi_lo$U10M[j]
    mi_hi$V10M[match_ix] <- mi_lo$V10M[j]
    mi_hi$U50M[match_ix] <- mi_lo$U50M[j]
    mi_hi$V50M[match_ix] <- mi_lo$V50M[j]
  }
  rm(j)
  # Interpolate between hourly u and v
  # !!!! Here, used linear, but may want to explore other options:
  # such as Guassian processes (recommended by Petar Djuric) >> ask Levi, Andy?
  mi_hi$U2M<-na_interpolation(mi_hi$U2M, option = "linear")
  mi_hi$V2M<-na_interpolation(mi_hi$V2M, option = "linear")
  mi_hi$U10M<-na_interpolation(mi_hi$U10M, option = "linear")
  mi_hi$V10M<-na_interpolation(mi_hi$V10M, option = "linear")
  mi_hi$U50M<-na_interpolation(mi_hi$U50M, option = "linear")
  mi_hi$V50M<-na_interpolation(mi_hi$V50M, option = "linear")
  
  for (j in 1:length(mi_hi$id)) {
    # Calculate Wind Speed and DIrection from u and v and BWA, BWAClass
    # uv2ddff for wind direction and velocity
    ddff2M <- uv2ddff(mi_hi$U2M[j],mi_hi$V2M[j])
    mi_hi$wind_speed2M[j] <- ddff2M$ff
    mi_hi$wind_dir3602M[j]<-ddff2M$dd
    
    ddff10M <- uv2ddff(mi_hi$U10M[j],mi_hi$V10M[j])
    mi_hi$wind_speed10M[j] <- ddff10M$ff
    mi_hi$wind_dir36010M[j]<-ddff10M$dd
    
    ddff50M <- uv2ddff(mi_hi$U50M[j],mi_hi$V50M[j])
    mi_hi$wind_speed50M[j] <- ddff50M$ff
    mi_hi$wind_dir36050M[j]<-ddff50M$dd
    
    # calculate windshear
    mi_hi$windshear10m[j]<-  abs(mi_hi$wind_speed10M[j]- mi_hi$wind_speed2M[j])
    mi_hi$windshear50m[j]<-  abs(mi_hi$wind_speed50M[j]- mi_hi$wind_speed10M[j])
    
    # bird-wind-angle
    bird_wind_angle2M<-abs(Lon360to180(mi_hi$bearingij[j]-mi_hi$wind_dir3602M[j]))
    bird_wind_angle10M<-abs(Lon360to180(mi_hi$bearingij[j]-mi_hi$wind_dir36010M[j]))
    
    # bwa class
    if (is.na(bird_wind_angle2M)) {
      mi_hi$bwa_class2M[j] <- NA
    }else{
      if (bird_wind_angle2M < 50) {
        mi_hi$bwa_class2M[j] <- "Head-Wind"
      }else if (bird_wind_angle2M > 130) {
        mi_hi$bwa_class2M[j] <- "Tail-Wind"
      }else{
        mi_hi$bwa_class2M[j] <- "Cross-Wind"
      }
    }
    mi_hi$bwa2M[j] <- bird_wind_angle2M
    
    if (is.na(bird_wind_angle10M)) {
      mi_hi$bwa_class10M[j] <- NA
    }else{
      if (bird_wind_angle10M < 50) {
        mi_hi$bwa_class10M[j] <- "Head-Wind"
      }else if (bird_wind_angle10M > 130) {
        mi_hi$bwa_class10M[j] <- "Tail-Wind"
      }else{
        mi_hi$bwa_class10M[j] <- "Cross-Wind"
      }
    }
    mi_hi$bwa10M[j] <- bird_wind_angle10M
  }
  
  spvec<-substr(birdi, 1, 4)
  mi_hi$species<-spvec
  
  
  # write individual bird file
  dir_i<-paste0('/Volumes/GoogleDrive/My Drive/THORNE_LAB/Data/Feldman_Analysis/Wind/wind_paired-data/MERRA-2/L0_paired/2019-2020_WAAL_incubation/indiv_birds/', birdi, '_wind_bwa.csv')
  write_csv(mi_hi, dir_i)
  
  rm(list=ls()[! ls() %in% c("wrap360", "uv2ddff", "Lon360to180", "birds", "m_lo", "m_hi", "i", "int_now")])
  
}