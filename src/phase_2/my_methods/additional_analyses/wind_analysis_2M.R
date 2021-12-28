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

getNcTime <- function(nc) {
  require(lubridate)
  ncdims <- names(nc$dim) #get netcdf dimensions
  timevar <- ncdims[which(ncdims %in% c("time", "Time", "datetime", "Datetime", "date", "Date"))[1]] #find time variable
  times <- ncvar_get(nc, timevar)
  if (length(timevar)==0) stop("ERROR! Could not identify the correct time variable")
  timeatt <- ncatt_get(nc, timevar) #get attributes
  timedef <- strsplit(timeatt$units, " ")[[1]]
  timeunit <- timedef[1]
  tz <- timedef[5]
  timestart <- strsplit(timedef[4], ":")[[1]]
  if (length(timestart) != 3 || timestart[1] > 24 || timestart[2] > 60 || timestart[3] > 60 || any(timestart < 0)) {
    cat("Warning:", timestart, "not a valid start time. Assuming 00:00:00\n")
    warning(paste("Warning:", timestart, "not a valid start time. Assuming 00:00:00\n"))
    timedef[4] <- "00:00:00"
  }
  if (! tz %in% OlsonNames()) {
    cat("Warning:", tz, "not a valid timezone. Assuming UTC\n")
    warning(paste("Warning:", timestart, "not a valid start time. Assuming 00:00:00\n"))
    tz <- "UTC"
  }
  timestart <- ymd_hms(paste(timedef[3], timedef[4]), tz=tz)
  f <- switch(tolower(timeunit), #Find the correct lubridate time function based on the unit
              seconds=seconds, second=seconds, sec=seconds,
              minutes=minutes, minute=minutes, min=minutes,
              hours=hours,     hour=hours,     h=hours,
              days=days,       day=days,       d=days,
              months=months,   month=months,   m=months,
              years=years,     year=years,     yr=years,
              NA
  )
  suppressWarnings(if (is.na(f)) stop("Could not understand the time unit format"))
  timestart + f(times)
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

# =======================================================================================

# START HERE: 
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
island <- "tern"
species <- "BFAL"
track.data <- bt

# ======================================================================================

## cycle through each year/ month / day of study period
## select track points that occur on each year/ month/ day and create a buffer around that point(s) 
## then select raster from that date and extract raster values to buffer (unlist values first)
library(stars)
library(lubridate)

for (yrs in 2008:2012){
  
  for (mths in 1:12) {
    ## this sets the number of days per month 
    if (mths == 1 | mths==3 | mths == 5 | mths == 7 | mths == 8 | mths == 10 | mths == 12) {numdays<-31}
    if (mths == 2) {numdays <- 28} 
    if (mths == 4 | mths==6 | mths == 9 | mths == 11) {numdays<-30}
    datesplit <- "20" 
    # yr<-strsplit(as.character(yrs),split=datesplit)[[1]][2]
    
    ##this puts a 0 before days and months <10 so that each month and day has two numbers 
    if (mths<10) {mths<-paste(0,as.character(mths),sep="")}
    
    wind_data <- data.frame(id=as.character(),
                            dtime=as.Date(character()), 
                            x=numeric(),
                            y=numeric(),
                            track=as.character(),
                            yrs=numeric(),
                            mths=numeric(),
                            dys=numeric(),
                            wind_u=numeric(),
                            wind_v=numeric(),
                            wind_speed2M=numeric(),
                            wind_direction2M=numeric(),
                            stringsAsFactors=FALSE) 
    names(wind_data)<-names(track.data)
    n_ <- ncol(wind_data)
    names(wind_data)[n_]<-paste("wind_direction2M")
    names(wind_data)[n_ -1]<-paste("wind_speed2M")
    names(wind_data)[n_ -2]<-paste("daily_avg_v2m")
    names(wind_data)[n_ -3]<-paste("daily_avg_u2m")
    names(wind_data)[n_ -4]<-paste("dys")
    names(wind_data)[n_ -5]<-paste("mths")
    names(wind_data)[n_ -6]<-paste("yrs")
    
    for (dys in 1:numdays) {
      if (dys<10) {dys<-paste(0,as.character(dys),sep="")}
      
    print(dys)
    print(mths)
    print(yrs)
    
    date.temp <- mdy(paste(mths,dys,yrs,sep=" "))
    points <- track.data
    points.temp <- points[which(points$dtime==date.temp),]
    points <- points.temp[,c(3,4)]
    
    if (nrow(points)>0){

      setwd("/Volumes/Samsung_T5/project_data/spatial_segregation/data/oceanographic/wind_MERRA_raw/left_all/")
      if (yrs==2008 | yrs==2009 | yrs==2010){filename_left<-paste("MERRA2_300.tavg1_2d_slv_Nx.",yrs,mths,dys,".SUB.nc",sep="")}
      if (yrs==2011 | yrs==2012){filename_left<-paste("MERRA2_400.tavg1_2d_slv_Nx.",yrs,mths,dys,".SUB.nc",sep="")}
    
      ## this reads a netcdf, converts to raster, converts to dataframe
      wind_ncdf <- read_ncdf(filename_left)
      wind_u2<-as(wind_ncdf[2,,,], "Raster")
      wind_u.df= raster::as.data.frame(wind_u2, xy = TRUE)
      wind_U2M <- wind_u.df
      
      # same thing again but right side of antimeridian
      setwd("/Volumes/Samsung_T5/project_data/spatial_segregation/data/oceanographic/wind_MERRA_raw/right_all/")
      if (yrs==2008 | yrs==2009 | yrs==2010){filename_right<-paste("MERRA2_300.tavg1_2d_slv_Nx.",yrs,mths,dys,".SUB.nc",sep="")}
      if (yrs==2011 | yrs==2012){filename_right<-paste("MERRA2_400.tavg1_2d_slv_Nx.",yrs,mths,dys,".SUB.nc",sep="")}
      
      wind_ncdf2 <- read_ncdf(filename_right)
      wind_u2v2<-as(wind_ncdf2[2,,,], "Raster")
      wind_u.df2= raster::as.data.frame(wind_u2v2, xy = TRUE)
      wind_U2M2 <- wind_u.df2
      
      #merge both sides, convert to raster
      wind_U2M <- rbind(wind_U2M, wind_U2M2) 
      u_stack2M<-rasterFromXYZ(wind_U2M)

      # extract from this raster, which should just include points that are present in "points" (movement data)
      wind_extracted <- extract(u_stack2M, points,buffer=300000, fun=mean)
      
      # bind extracted values
      temp2 <- cbind(points.temp,yrs, mths, dys, wind_extracted)
      
      # average 24 hours into one daily average
      
      temp.means <- rowMeans(temp2[,9:32])
      temp2$daily_avg_u2m <- temp.means
      temp2 <- temp2[,c(1:8,33)]
      
      # REPEAT AGAIN TO GET V2
      
      setwd("/Volumes/Samsung_T5/project_data/spatial_segregation/data/oceanographic/wind_MERRA_raw/left_all/")
      if (yrs==2008 | yrs==2009 | yrs==2010){filename_left<-paste("MERRA2_300.tavg1_2d_slv_Nx.",yrs,mths,dys,".SUB.nc",sep="")}
      if (yrs==2011 | yrs==2012){filename_left<-paste("MERRA2_400.tavg1_2d_slv_Nx.",yrs,mths,dys,".SUB.nc",sep="")}
      
      ## this reads a netcdf, converts to raster, converts to dataframe
      wind_ncdf <- read_ncdf(filename_left)
      wind_v2<-as(wind_ncdf[5,,,], "Raster")
      wind_v.df= raster::as.data.frame(wind_v2, xy = TRUE)
      wind_V2M <- wind_v.df
      
      # same thing again but right side of antimeridian
      setwd("/Volumes/Samsung_T5/project_data/spatial_segregation/data/oceanographic/wind_MERRA_raw/right_all/")
      if (yrs==2008 | yrs==2009 | yrs==2010){filename_right<-paste("MERRA2_300.tavg1_2d_slv_Nx.",yrs,mths,dys,".SUB.nc",sep="")}
      if (yrs==2011 | yrs==2012){filename_right<-paste("MERRA2_400.tavg1_2d_slv_Nx.",yrs,mths,dys,".SUB.nc",sep="")}
      
      wind_ncdf2 <- read_ncdf(filename_right)
      wind_v2v2<-as(wind_ncdf2[5,,,], "Raster")
      wind_v2.df2= raster::as.data.frame(wind_v2v2, xy = TRUE)
      wind_V2M2 <- wind_v2.df2
      
      #merge both sides, convert to raster
      wind_V2M <- rbind(wind_V2M, wind_V2M2) 
      v_stack2M<-rasterFromXYZ(wind_V2M)
      
      # extract from this raster, which should just include points that are present in "points" (movement data)
      wind_extracted <- extract(v_stack2M, points,buffer=300000, fun=mean)
      
      # bind extracted values
      temp3 <- cbind(points.temp,yrs, mths, dys, wind_extracted)
      
      # calculate daily mean from 24 hourly readings
      temp.means <- rowMeans(temp3[,9:32])
      temp2$daily_avg_v2m <- temp.means
      
      # now calc wind speed (ff) and direction (dd), add to main df
      holder_dd <- vector(mode="numeric")
      holder_ff <- vector(mode="numeric")
      for (i in 1:nrow(temp2)){
        add_temp <- uv2ddff(temp2[i,9], temp2[i,10])
        add_dd <- add_temp[1,1]
        add_ff <- add_temp[1,2]
        holder_dd <- append(holder_dd,add_dd)
        holder_ff <- append(holder_ff,add_ff)
      }
      temp2$wind_speed2M <- holder_ff
      temp2$wind_direction2M <- holder_dd
      wind_data <- rbind(wind_data,temp2)
      
      # IMPORTANT - LAYER 1-24 REPRESENTS HOURS OF THE DAY, YOU NEED TO AVERAGE INTO A SINGLE VALUE!!! 
     }
    } 
    setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/")
    #output.name2 <-  paste("Daily_env_data_4km_",yrs, "to", mths,".txt",sep="")
    output.name2 <-  paste(island,species,"_wind2M_extracted",yrs,mths,".txt",sep="")
    write.csv(wind_data, as.character((output.name2)))
  }
}






















 # load files back in
setwd("/Volumes/Samsung_T5/project_data/spatial_segregation/data/oceanographic/wind")
a <- list.files(getwd())
for (i in 1:length(a)){
  load(a[i])
}

load(a[4])
load(a[5])
load(a[10])
load(a[11])

wind_U2M <- rbind(wind_U2M, wind_U2M2) ###
save(wind_U2M,file="wind_U2M_final.Rdata")
load("wind_U2M_final.Rdata")
wind_U10M <- rbind(wind_U10M, wind_U10M2)
wind_U50M <- rbind(wind_U50M, wind_U50M2)

wind_V2M <- rbind(wind_V2M, wind_V2M2) ###
wind_V10M <- rbind(wind_V10M, wind_V10M2)
wind_V50M <- rbind(wind_V50M, wind_V50M2)

# convert to rasters
u_stack2M<-brick(wind_U2M)
v_stack2M<-rasterFromXYZ(wind_V2M)
u_stack10M<-rasterFromXYZ(wind_U10M)
v_stack10M<-rasterFromXYZ(wind_V10M)
u_stack50M<-rasterFromXYZ(wind_U50M)
v_stack50M<-rasterFromXYZ(wind_V50M)

save(u_stack2M,file="u_stack2M.Rdata")
save(v_stack2M,file="v_stack2M.Rdata")

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