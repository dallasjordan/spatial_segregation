######################################
# Set Environment
######################################
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
nc_dir<-'/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Analysis/Wind/wind_datasets/BirdIsland20182021/ERA5_SingleLevels_10m/'
setwd(nc_dir)

wind_files <- list.files(pattern='*.nc') 

wind_t1 <- read_ncdf(wind_files[1])
wind_t1 # Make sure you got the right stuff!
times_t1<-st_get_dimension_values(wind_t1, "time")  # stores times from each file 
lat_t1 <- st_get_dimension_values(wind_t1, "latitude")
lon_t1 <- st_get_dimension_values(wind_t1, "longitude")
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

####################################################################################################
# 2. Read in hourly bird locations and gather wind data
####################################################################################################

int_now <- 3600 # seconds 

setwd(paste0('/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Bird_Island/2019_2020/Tag_Data/L1_cleaned-data/GPS_L1_3_interpolated/', int_now, 's/'))
files<-list.files() # GPS files 

# Append all individual bird lat lon files to one file ------
for (i in 1:length(files)) {
  mi<-read.csv(files[i])
  if (i==1) {
    m<-mi
  }else{
    m<-rbind(m,mi)
  }
} 

# Add speed, distance, bearing columns to m -----------------
m$ground_speed_kmHr <- NA
m$distanceij_km <- NA
m$bearingij <- NA
m$u<- NA
m$v<- NA
m$wind_speed <- NA
m$wind_dir360 <- NA
m$bwa <- NA
m$bwa_class <- NA # head tail cross 

trips<- unique(m$tripID)
for (i in 1:length(trips)) {
  tripi<-m[m$tripID==trips[i],]
  
  for (j in 1:length(tripi$id)-1) {
    hour_int<- int_now/3600 # int_now is in seconds (3600 seconds in one hour)
    tripi$distanceij_km[j]<-swfscMisc::distance(tripi$lat[j],tripi$lon[j],tripi$lat[j+1],tripi$lon[j+1], units = "km", method = "haversine")
    tripi$ground_speed_kmHr[j]<-tripi$distanceij_km[j]/hour_int    
    tripi$bearingij[j]<-as.numeric(bearing(tripi$lat[j],tripi$lon[j],tripi$lat[j+1],tripi$lon[j+1])[1])
  }
  
  m[m$tripID==trips[i],]<-tripi
  
}  #speed, distance, bearing between points 

# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
# Loop through m and add wind information: u, v, velocity, direction, Bird-Wind-Angle
# ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
for ( j in 1:length(m$id)) {
  
  timej <- as.POSIXct(m$datetime[j], format = "%Y-%m-%d %H:%M:%S" , tz = "UTC")
  timej_num <- as.numeric(timej)
  # Find index of current_time in all times. Use that index to pull out relevant raster layer. 
  raster_dt_index <- which(abs(all_times_num-timej_num) == min(abs(all_times_num-timej_num))) # taking difference between all gps points and take the min: tells which layer to pull out and isolate 
  
  # Isolate u and v rasters at time j
  ustack_timej <- subset(u_stack, raster_dt_index, drop=TRUE) 
  vstack_timej <- subset(v_stack, raster_dt_index, drop=TRUE)
  
  # isolate coordinates
  xy_j <- as.data.frame(cbind(m$lon[j], m$lat[j]))
  colnames(xy_j) <- c("lon","lat")
  xy_j$lon <- Lon360to180(xy_j$lon) 
  
  # Extract u and v components for time j at location x and y
  u_j <- extract(ustack_timej, xy_j)
  v_j <- extract(vstack_timej, xy_j)
  m$u[j]<- u_j
  m$v[j]<- v_j
  # -----------------------------------------------------------------------------
  # Get Wind Direction and Wind Velocity from U and V Components:
  # -----------------------------------------------------------------------------
  # Requires mathematical to meteorological adjustment:
  # Important: our wind vectors were given in mathematical notation, not in meteorological convention. 
  # Need to adjust by 270* - this is what uv2ddff does.
  # http://colaweb.gmu.edu/dev/clim301/lectures/wind/wind-uv
  # from uv2ddff: Wind direction is either in meteorological degrees (0 from North, from 90 East,
  # 180 from South, and 270 from West) or in mathematical radiant if input \code{rad = TRUE}.
  
  # HIGHLIGHT AUXILLARY FUNCTION AND RUN 
  res<-uv2ddff(u_j, v_j) 
  m$wind_speed[j] <- res$ff
  m$wind_dir360[j]<-res$dd
  
  # BIRD WIND ANGLE
  # note!!! Bearing is on same compass as wind direction (0-360) BUT it's going IN the direction (while wind is coming FROM that direction). In other words, wind with a direction of 2 degrees is coming from the north. Whereas a bird with a bearing of 2 degrees is going TO the north. 
  bird_bearing <- m$bearingij[j]
  bird_wind_angle<-abs(Lon360to180(bird_bearing-m$wind_dir360[j])) # bird bearing 0-180, wind 0-360
  
  # Definitions from Spear and Ainley 1997
  # Headwind  = abs(bird-wind): 0-59 
  # Crosswind = abs(bird-wind): 60-119 
  # Tailwind =  abs(bird-wind): 120-180 
  
  if (is.na(bird_wind_angle)) {
    m$bwa_class[j] <- NA
  }else{
    if (bird_wind_angle < 50) {
      bwaClass <- "Head-Wind"
    }else if (bird_wind_angle > 130) {
      bwaClass <- "Tail-Wind"
    }else{
      bwaClass <- "Cross-Wind"
    }
  }

  m$bwa[j] <- bird_wind_angle
  m$bwa_class[j] <- bwaClass
  
}
spvec <- substr(m$id,1,4)
m<-m %>% mutate(species=spvec)

newm<-m %>% dplyr::select(c("species",colnames(m)[c(1:14)]))

write_csv(m, '/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Analysis/Wind/wind_analysis/bird-latlon_with-wind-bwa/hourly/allbirds_hourly_bwa_ERA5-Single-Levels-10m.csv')


plot.windrose(spd = m$wind_speed[m$species=="BBAL"],
              dir = m$bwa[m$species=="BBAL"])

# 

####################################################################################################
# 2. Interpolate wind data for 30 sec data
####################################################################################################

# Import 30s lat lon data and append individual bird files into one file. 
int_now <- 30

setwd(paste0('/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Bird_Island/2019_2020/Tag_Data/L1_cleaned-data/GPS_L1_3_interpolated/', int_now, 's/'))
files<-list.files()

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
m_hi$u <- NA
m_hi$v <- NA
m_hi$wind_speed<- NA
m_hi$wind_dir360<- NA
m_hi$bwa<- NA
m_hi$bwa_class<- NA

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
m_lo <- read.csv('/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Analysis/Wind/wind_analysis/bird-latlon_with-wind-bwa/hourly/allbirds_hourly_bwa_ERA5-SingleLevels-10m.csv')

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
    mi_hi$u[match_ix] <- mi_lo$u[j]
    mi_hi$v[match_ix] <- mi_lo$v[j]
  }
  rm(j)
  # Interpolate between hourly u and v
  # !!!! Here, used linear, but may want to explore other options:
  # such as Guassian processes (recommended by Petar Djuric) >> ask Levi, Andy?
  mi_hi$u<-na_interpolation(mi_hi$u, option = "linear") # linear interpolation 
  mi_hi$v<-na_interpolation(mi_hi$v, option = "linear")
  
  
  for (j in 1:length(mi_hi$id)) {
    # Calculate Wind Speed and DIrection from u and v and BWA, BWAClass
    # uv2ddff for wind direction and velocity
    ddff <- uv2ddff(mi_hi$u[j],mi_hi$v[j])
    mi_hi$wind_speed[j] <- ddff$ff
    mi_hi$wind_dir360[j]<-ddff$dd
    
    # bird-wind-angle
    bird_wind_angle<-abs(Lon360to180(mi_hi$bearingij[j]-mi_hi$wind_dir360[j]))
    
    # bwa class
    if (is.na(bird_wind_angle)) {
      mi_hi$bwa_class[j] <- NA
    }else{
      if (bird_wind_angle < 50) {
        mi_hi$bwa_class[j] <- "Head-Wind"
      }else if (bird_wind_angle > 130) {
        mi_hi$bwa_class[j] <- "Tail-Wind"
      }else{
        mi_hi$bwa_class[j] <- "Cross-Wind"
      }
    }
    
    mi_hi$bwa[j] <- bird_wind_angle
    
  }
  
  col_list<-colnames(m_lo)[c(1:14)]
  mi_hi<-mi_hi %>% dplyr::select(col_list)
  spvec<-substr(birdi, 1, 4)
  mi_hi$species<-spvec
  

  # write individual bird file
  dir_i<-paste0('/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Analysis/Wind/wind_analysis/bird-latlon_with-wind-bwa/30s/indiv_birds/', birdi, '_wind_bwa.csv')
  write_csv(mi_hi, dir_i)
  
  rm(list=ls()[! ls() %in% c("wrap360", "uv2ddff", "plot.windrose", "Lon360to180", "birds", "m_lo", "m_hi", "i", "int_now")])
  
}












##########################################################################################
# AUXULIARY FUNCTIONS
##########################################################################################



plot.windrose <- function(data,
                          spd,
                          dir,
                          spdres = 2,
                          dirres = 30,
                          spdmin = 2,
                          spdmax = 20,
                          spdseq = NULL,
                          palette = "YlGnBu",
                          countmax = NA,
                          debug = 0){
  
  
  # Look to see what data was passed in to the function
  if (is.numeric(spd) & is.numeric(dir)){
    # assume that we've been given vectors of the speed and direction vectors
    data <- data.frame(spd = spd,
                       dir = dir)
    spd = "spd"
    dir = "dir"
  } else if (exists("data")){
    # Assume that we've been given a data frame, and the name of the speed 
    # and direction columns. This is the format we want for later use.    
  }  
  
  # Tidy up input data ----
  n.in <- NROW(data)
  dnu <- (is.na(data[[spd]]) | is.na(data[[dir]]))
  data[[spd]][dnu] <- NA
  data[[dir]][dnu] <- NA
  
  # figure out the wind speed bins ----
  if (missing(spdseq)){
    spdseq <- seq(spdmin,spdmax,spdres)
  } else {
    if (debug >0){
      cat("Using custom speed bins \n")
    }
  }
  # get some information about the number of bins, etc.
  n.spd.seq <- length(spdseq)
  n.colors.in.range <- n.spd.seq - 1
  
  # create the color map
  spd.colors <- colorRampPalette(brewer.pal(min(max(3,
                                                    n.colors.in.range),
                                                min(9,
                                                    n.colors.in.range)),                                               
                                            palette))(n.colors.in.range)
  
  if (max(data[[spd]],na.rm = TRUE) > spdmax){    
    spd.breaks <- c(spdseq,
                    max(data[[spd]],na.rm = TRUE))
    spd.labels <- c(paste(c(spdseq[1:n.spd.seq-1]),
                          '-',
                          c(spdseq[2:n.spd.seq])),
                    paste(spdmax,
                          "-",
                          max(data[[spd]],na.rm = TRUE)))
    spd.colors <- c(spd.colors, "grey50")
  } else{
    spd.breaks <- spdseq
    spd.labels <- paste(c(spdseq[1:n.spd.seq-1]),
                        '-',
                        c(spdseq[2:n.spd.seq]))    
  }
  data$spd.binned <- cut(x = data[[spd]],
                         breaks = spd.breaks,
                         labels = spd.labels,
                         ordered_result = TRUE)
  # clean up the data
  data. <- na.omit(data)
  
  # figure out the wind direction bins
  dir.breaks <- c(-dirres/2,
                  seq(dirres/2, 360-dirres/2, by = dirres),
                  360+dirres/2)  
  dir.labels <- c(paste(360-dirres/2,"-",dirres/2),
                  paste(seq(dirres/2, 360-3*dirres/2, by = dirres),
                        "-",
                        seq(3*dirres/2, 360-dirres/2, by = dirres)),
                  paste(360-dirres/2,"-",dirres/2))
  # assign each wind direction to a bin
  dir.binned <- cut(data[[dir]],
                    breaks = dir.breaks,
                    ordered_result = TRUE)
  levels(dir.binned) <- dir.labels
  data$dir.binned <- dir.binned
  
  # Run debug if required ----
  if (debug>0){    
    cat(dir.breaks,"\n")
    cat(dir.labels,"\n")
    cat(levels(dir.binned),"\n")       
  }  
  
  # deal with change in ordering introduced somewhere around version 2.2
  if(packageVersion("ggplot2") > "2.2"){    
    cat("Hadley broke my code\n")
    data$spd.binned = with(data, factor(spd.binned, levels = rev(levels(spd.binned))))
    spd.colors = rev(spd.colors)
  }
  
  # create the plot ----
  p.windrose <- ggplot(data = data,
                       aes(x = dir.binned,
                           fill = spd.binned)) +
    geom_bar() + 
    scale_x_discrete(drop = FALSE,
                     labels = waiver()) +
    coord_polar(start = -((dirres/2)/360) * 2*pi) +
    scale_fill_manual(name = "Wind Speed (m/s)", 
                      values = spd.colors,
                      drop = FALSE) +
    #theme_bw() +
    theme(axis.title.x = element_blank(),
          #panel.border = element_rect(colour = "blank"),
          panel.grid.major = element_line(colour="grey65"))
  
  # adjust axes if required
  if (!is.na(countmax)){
    p.windrose <- p.windrose +
      ylim(c(0,countmax))
  }
  
  # print the plot
  print(p.windrose)  
  
  # return the handle to the wind rose
  return(p.windrose)
}









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


