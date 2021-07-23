# probGLS implementation for Scott Shaffer's Midway albatross GLS data
# Dallas Jordan
# last edited: Sept. 2020

# 9/22 TO-DO
  # Fix algo inputs
  # Check scale up process
  # Input data
  # refined Bayesian points or daylogs? Daylogs.


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
library(dplyr)

# if you aren't working in your project 'spatial segregation'
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation")
filenames<- list.files()

# load in master raw data
raw_data <- read.csv("./data/raw_data.csv", stringsAsFactors = F)
unique <- unique(raw_data$Bird.ID)

# create list of partitioned bird_id dataframes

partition <- vector("list", length(unique))
for (i in 1:length(unique)) {
  add <- as.data.frame(dplyr::filter(raw_data, Bird.ID==unique[i]))
  partition[[i]] <- add
}
############################################################################################################

# Load in data set from partitioned lists for first tag, will have to repeat manually for each tag. ## ADD COMMENT ON BACKCALCULATIONS
# Difficult to throw this into a for-loop, therefore I change line 55 and 69 to the for each individual tag ID in my parititioned list. 

# START HERE for each bird.
bc <- dplyr::select(partition[[20]], -Bird.ID, -SST) # Change partition [[ ]] for each individual bird. 
bc <- bc[,c(1,3,2)] # reorders columns for neatness
bc$Date.and.Time <- strptime(bc$Date.and.Time, format="%m/%d/%y", tz="Pacific/Midway") # Parses string as a date format
bc$Date.and.Time <- as.Date(bc$Date.and.Time, format= "%F", tz="Pacific/Midway") # Changes format of the date to a format usable by rise_set
bc <- bc %>% rename(date = Date.and.Time,lat = Raw.Latitude,lon = Raw.Longitude)
rise_set <- getSunlightTimes(data=bc, keep= c("sunrise", "sunset","nauticalDawn"), tz="Pacific/Midway")
bc$sunrise <-rise_set[,4]
bc$sunset <- rise_set[,5]

# Filter out any places where no data. Some tags recorded events with no location estimates. These lines are filtered out here. 
bc <- bc[!is.na(bc$sunrise),]
bc <- bc[!is.na(bc$sunset),]

# Extract SST and datetime from paritioned dataframe lists, creates a new SST dataframe for use in algorithm. Dates are reformatted to be usable in probGLS. 
SST <- dplyr::select(partition[[3]], -Bird.ID, -Raw.Latitude, -Raw.Longitude) # Change partition [[ ]] for each individual bird. 
SST   <- SST[complete.cases(SST),]
SST$dtime <- as.POSIXct(strptime(SST[,1], format = "%m/%d/%y"),tz='UTC')
SST       <- SST[!is.na(SST$dtime),]
SST$Date.and.Time <- strptime(SST$Date.and.Time, format="%m/%d/%y", tz="Pacific/Midway")
SST$Date.and.Time <- as.Date(SST$Date.and.Time, format= "%F", tz="Pacific/Midway")



# JULY 30 REVERT TO THIS IF YOUR SCALE UP DOESNT WORK - DON'T JUST DELETE, SAVE THIS

# # Load in trimmed data set for first tag, will have to repeat manually for each tag. 
# bc_bird110801 <- read.csv("./data/1108_01.csv", stringsAsFactors = F)
# bc_bird110801 <- dplyr::select(bc_bird110801, -Bird.ID, -SST)
# bc_bird110801 <- bc_bird110801[,c(1,3,2)]
# bc_bird110801$Date.and.Time <- strptime(bc_bird110801$Date.and.Time, format="%m/%d/%y", tz="Pacific/Midway")
# bc_bird110801$Date.and.Time <- as.Date(bc_bird110801$Date.and.Time, format= "%F", tz="Pacific/Midway")
# bc_bird110801 <- bc_bird110801 %>% rename(date = Date.and.Time,lat = Raw.Latitude,lon = Raw.Longitude)
# rise_set <- getSunlightTimes(data=bc_bird110801, keep= c("sunrise", "sunset","nauticalDawn"), tz="Pacific/Midway")
# bc_bird110801$sunrise <-rise_set[,4]
# bc_bird110801$sunset <- rise_set[,5]
# 
# # Filter out any places where no data
# bc_bird110801 <- bc_bird110801[!is.na(bc_bird110801$sunrise),]
# bc_bird110801 <- bc_bird110801[!is.na(bc_bird110801$sunset),]
# 
# # Create SST and datetime dataframe
# SST_bird110801 <- read.csv("./1108_01.csv", stringsAsFactors = F)
# SST_bird110801 <- dplyr::select(SST_bird110801, -Bird.ID, -Raw.Latitude, -Raw.Longitude)
# SST_bird110801$Date.and.Time <- strptime(SST_bird110801$Date.and.Time, format="%m/%d/%y", tz="Pacific/Midway")
# SST_bird110801$Date.and.Time <- as.Date(SST_bird110801$Date.and.Time, format= "%F", tz="Pacific/Midway")

############################################################################################################


# lotek_to_dataframe to convert your calculated riseset into a format usable by probGLS
# necessary formats are as follows:
# date	column "TimeS" in daylog file as as.Date() [DONE]
# sunrise column "Sunrise" in daylog file as as.character()
# sunset column "Sunset" in daylog file as as.character()
# lon/lat is in numeric [DONE]

# Fill empty cells with NA
# bc_bird110801 <- as.data.frame(apply(bc_bird110801, 2, function(x) gsub("^$|^ $", NA, x)))

# Filter out rows containing all NAs
# bc_bird110801 <- bc_bird110801 %>% filter_all(any_vars(!is.na(.)))

# ADD COMMENTARY - changed rise_set data to necessary format for probGLS (as.character())
sunrise <- as.character(sapply(strsplit(as.character(bc$sunrise), split = " "), "[", 2))
sunrise.1 <- sapply(strsplit(sunrise, split = ":"), "[", 1)
sunrise.2 <-  sapply(strsplit(sunrise, split = ":"), "[", 2)
sunrise <- paste(sunrise.1, sunrise.2, sep=":")
sunset <- as.character(sapply(strsplit(as.character(bc$sunset), split = " "), "[", 2))
sunset.1 <- sapply(strsplit(sunset, split = ":"), "[", 1)
sunset.2 <-  sapply(strsplit(sunset, split = ":"), "[", 2)
sunset <- paste(sunset.1, sunset.2, sep=":")
date <- as.Date(as.character(bc$date), format = "%Y-%m-%d")

TRLat <- as.numeric(as.character(bc$lat))
TRLon <- as.numeric(as.character(bc$lon))


#######

# can plot to visualize onboard algorithm calculated points

for (i in 1:length(TRLon)){
  if (TRLon[i]<0){
    TRLon[i] = TRLon[i]+360
  }
}

data(world2HiresMapEnv)
plot(TRLon,TRLat ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(150,230), ylim = c(0, 70), bty = "n")
map('world2Hires', xlim= c(150,230), ylim = c(0, 70),  col = "grey", add = T)
points(TRLon,TRLat, pch=16, col="cornflowerblue", type = "o")

########

# Run lotek_to_dataframe to convert your calculated riseset into a format usable by probGLS

trn <- lotek_to_dataframe( date = date, 
                           sunrise = sunrise, 
                           sunset = sunset, 
                           TRLon = TRLon, 
                           TRLat = TRLat)
trn <- trn[!is.na(trn$tSecond),]
trn$keep <- loessFilter(trn, k = 3, plot = T) # run a Loess filter to take out unrealistic poiints. 


############################################################################################################

# Run algorithm

# download environmental data ----

# download yearly NetCDF files for (replace YEAR with appropriate number): 
# daily mean SST                   -> 'sst.day.mean.YEAR.v2.nc'
# daily SST error                  -> 'sst.day.err.YEAR.v2.nc'
# daily mean sea ice concentration -> 'icec.day.mean.YEAR.v2.nc'

# from:
# http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html
# and place all into the same folder

# Also, download the land mask file: 'lsmask.oisst.v2.nc' from the same directory 
# and place it in the same folder as all the other NetCDF files

# Specify parameters for prob_algorithm 
start <- as.POSIXct("2007-12-28 00:00", tz="Pacific/Midway") # Time when tagged
end   <- as.POSIXct("2008-12-04 08:55", tz="Pacific/Midway") # Time when tag was retrieved. 
midway.coords <- c(182.65,28.2)
tw    <- twilight_error_estimation() # Position estimates require error in twilight calculation; this baked-in function gives an error distribution for the algorithm 
sen   <- sst_deduction(datetime=SST$dtime, temp=SST$SST, temp.range=(c(-1.52,36.3))) # calculates proportion of time spent in water as well as creating a dataframe of recorded SST. 
abline(h=c(-2,19),lty=2,col="orange") # visualizes range of SST 
sst.sd<- sd(sen$sst) # takes standard deviation of recorded SST
rownames(sen)<-NULL # gets rid of rownames 

pr   <- prob_algorithm(trn                         = trn, 
                       sensor                      = NULL,
                       act                         = NULL, 
                       tagging.date                = start, 
                       retrieval.date              = end, 
                       loess.quartile              = NULL, 
                       tagging.location            = midway.coords, 
                       particle.number             = 800, 
                       iteration.number            = 60,
                       sunrise.sd                  = c(2.49,0.94,4.98),
                       sunset.sd                   = c(2.49,0.94,4.98),
                       range.solar                 = c(-10,3),
                       boundary.box                = c(150,-140,10,60),
                       days.around.spring.equinox  = c(10,10), 
                       days.around.fall.equinox    = c(10,10),
                       speed.dry                   = c(10,6,25), # Cite Weimerskirsch paper
                       speed.wet                   = c(1,1.3,5), 
                       sst.sd                      = sst.sd,  # cannot be NULL - weird       
                       max.sst.diff                = 3,          
                       east.west.comp              = F,
                       land.mask                   = T, 
                       ice.conc.cutoff             = 0, 
                       wetdry.resolution           = 1,
                       NOAA.OI.location            = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/environment_data")


plot_map(pr) # This function does not work for this data because BFAL and LAAL cross Pacific meridian (plot_map is Atlantic centric. Below is code for a Pacific-centric mapping). 
# Ultimately, unnecessary if importing into ArcMap, but can be helpful for data visualization. 

#######################################################################################################################################################

# Unfinished attempt to edit probGLS internal mapping function to display Pacific-centric view
# trace("plot_map", edit=T) to open editor

test_coords <- as.data.frame(pr[[1]]@coords)
for (i in 1:length(test_coords$lon)){
  if (test_coords$lon[i]<0){
    test_coords$lon[i] = test_coords$lon[i]+360
  }
}

data(world2HiresMapEnv)
map('world2Hires', ylim=c(10,70), xlim= c(130,235))
plot(test_coords, ylab = "Latitude", xlab = "Longitude", add=T)

### 

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}

pr[[1]]@coords[,1] <- convlon(pr[[1]]@coords[,1])
pr[[2]]@coords[,1] <- convlon(pr[[2]]@coords[,1])
#newextent<-extent(130,135,10,70)
#extent(pr[[1]])<-newextent

p_map <- function (pr) 
{
  cc <- as.numeric(unlist(strsplit(as.character(pr[[4]][4,2]), "[ ]")))
  opar <- par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
  plot(pr[[1]], ylim=c(10,70), xlim= c(130,235), col = "white", ylab = "Latitude", xlab = "Longitude")
  axis(1)
  axis(2)
  plot(pr[[1]], ylim=c(10,70), xlim= c(130,235), col = colorRampPalette(c("grey90", "grey50"))(nrow(pr[[2]]))[pr[[1]]$step], 
       add = T, pch = 19, cex = 0.3)
  map('world2Hires', ylim=c(10,70), xlim= c(130,235), add = T, col = 4, lwd = 0.5)
  mm2 <- pr[[2]]
  lines(mm2$lon, mm2$lat, col = 1)
  lines(c(mm2$lon[1], cc[1]), c(mm2$lat[1], cc[2]), lty = 3)
  lines(c(mm2$lon[nrow(mm2)], cc[1]), c(mm2$lat[nrow(mm2)], 
                                        cc[2]), lty = 3)
  points(mm2$lon, mm2$lat, cex = 1, pch = 21, bg = colorRampPalette(c("yellow", 
                                                                      "darkred"))(nrow(pr[[2]])))
  mm3 <- mm2[is.na(mm2$median.sun.elev), ]
  points(mm3$lon, mm3$lat, cex = 0.7, pch = 3)
  points(cc[1], cc[2], cex = 2.1, pch = 21, col = "violet", 
         lwd = 2)
  par(opar)
}

p_map(pr)

############################################################################################################

# Extract coords. Replot using Pacific view mapper. 

xy <- coordinates(pr$`most probable track`)
lon_map<-xy[,1]
lat_map<-xy[,2]
export <- as.data.frame(coordinates(pr$`most probable track`))
write.csv(export,"export_track3.csv", row.names = FALSE)

# AT THIS POINT, PROBABLY WANT TO EXPORT AND DO IN ARCMAP
# END PRE-PROCESSING FOR INDIVIDUAL BIRD; RESTART AT TOP FOR NEXT BIRD. 

Lon180to360 <- function(lon){
  lon %% 360
}

plot(Lon180to360(lon_map), Lon180to360(lat_map)) # you will need to define max and min lons and lats for your extent
plot(wrld_simpl, add = TRUE, col="grey")
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE, col="grey")












































############################################################################################################

#NOT SURE IF I NEED THIS 7/14
#If you confirm that Lotek tag recorded times are in GMT, then apply sunrise/sunset time fix. Calculated sunrise/sunset are in UTC. -11 hours to convert to SST (Midway time)

utcsunrise<-chron(times=paste(as.character(rawd$sunrise),":00",sep=""))
utcsunset<-chron(times=paste(as.character(rawd$sunset),":00",sep=""))

#change sunrise sunset minutes to hst , so that theyre like LAT2500.
#turn into numeric hst minutes (like LAT2500)

#Rise:
sst.rise<-utcsunrise-chron(times="11:00:00")
r.hrs<-hours(sst.rise)
r.min<-minutes(sst.rise)
hourminrise<-sprintf("%02d%02d", r.hrs, r.min) 
RiseNumSST<-as.numeric(hourminrise)

#Set:w
sst.set<-utcsunset-chron(times="10:00:00")
s.hrs<-hours(hst.set)
s.min<-minutes(hst.set)
hourminset<-sprintf("%02d%02d", s.hrs, s.min) 
SetNumHST<-as.numeric(hourminset)