# Lotek twilight times with SGAT
# Michael Sumner

# 2014-09-01

# The SGAT package provides a range of methods for estimating geolocation given either estimated times of twilight or the light profiles from each twilight. In particular, SGAT provides
# the simple threshold method implemented in the GeoLight package
# the template or curve fitting method implemented in the tripEstimation package, and
# !!!! a new threshold based variant !!!! of the algorithm from the tripEstimation package.
# This third method is better suited to some tag data than the original tripEstimation, and especially here since we only have times of Sunrise and Sunset. Ultimately, the intention is for SGAT to supersede tripEstimation and GeoLight.

# In this document we will demonstrate the threshold based variant of the tripEstimation method.

# Data
# We have a Lotek "Day Log" file for two birds deployed at Tern Island. The first is a Laysan albatross deployed between March 2009  and January 2010 and the second a Black-footed albatross deployed between March 2010 and December 2010.

# First, read the raw data from the file ensuring we have sensible interpretation of the values.

library(SGAT)
library(BAStag)
library(chron)
library(pracma)
library(raster)
library(maptools)
library(rgeos)
library(maptools)


setwd("C:/Users/Melinda/Dropbox/Chapter_03_Finalize/data/Tern_GLS_LOTEK_allrawfiles_Rdata/DayLogs")
filenames<- list.files()

load(filenames[i])
filenames[i]

deploy.time <- as.POSIXct("2004-01-11 00:00:00", tz = "HST")
recov.date <- ISOdatetime(2004, 12, 1, 0, 0, 0, tz = "HST")
endix<-which.min(abs(recov.date-as.POSIXct(as.character(rawd$datetime),"%m/%d/%Y",tz="HST")))
rawd <- rawd[1:endix,]
rawd <- rawd[as.POSIXct(as.character(rawd$date),"%m/%d/%Y",tz="HST") >= deploy.time, ]


# adjust classes (have to make just like LAT tags.)

##------------------------------ Fix lats and lons-----------------------------------------------------------------
# lons and lats are initial locations from tag algorithm given in tag datafile

rawd$lon<-as.character(rawd$lon)
rawd$lat<-as.character(rawd$lat)

lon.split<-strsplit(rawd$lon, split=" ")
lat.split<-strsplit(rawd$lat, split=" ")

for (k in 1:length(lon.split)) 
{rawd$lon[k]<-ifelse(lon.split[k][[1]][2]=="West",paste("-",lon.split[k][[1]][1],sep=""),lon.split[k][[1]][1])}
for (k in 1:length(lat.split)) 
{rawd$lat[k]<-ifelse(lat.split[k][[1]][2]=="South",paste("-",lat.split[k][[1]][1],sep=""),lat.split[k][[1]][1])}

rawd$lon<-as.numeric(rawd$lon)
rawd$lat<-as.numeric(rawd$lat)

# Record some details that we know about this tag, the deployment date will be used to ignore data prior to this time.
lon.home <- -166.29 + 360
lat.home <- 23.87

## maximum range in lon/lat (for land mask plus maximum bound of track)
lon.range <- c(140, 235)
lat.range <- c(20, 62)

# replace NA lon with previous lon
rawd$lon[which(is.na(rawd$lon))]<-rawd$lon[which(is.na(rawd$lon))-1]

#convert lons 0 to 360
rawd$lon[rawd$lon < 0] <- rawd$lon[rawd$lon < 0] + 360

plot(rawd$lon, rawd$lat)

##------------------------- FIX SUNRISE SUNSET -----------------------------------------------------------------------------

utcsunrise<-chron(times=paste(as.character(rawd$sunrise),":00",sep=""))
utcsunset<-chron(times=paste(as.character(rawd$sunset),":00",sep=""))

#change sunrise sunset minutes to hst , so that theyre like LAT2500.
#turn into numeric hst minutes (like LAT2500)

#Rise:
hst.rise<-utcsunrise-chron(times="10:00:00")
r.hrs<-hours(hst.rise)
r.min<-minutes(hst.rise)
hourminrise<-sprintf("%02d%02d", r.hrs, r.min) 
RiseNumHST<-as.numeric(hourminrise)

#Set:
hst.set<-utcsunset-chron(times="10:00:00")
s.hrs<-hours(hst.set)
s.min<-minutes(hst.set)
hourminset<-sprintf("%02d%02d", s.hrs, s.min) 
SetNumHST<-as.numeric(hourminset)


rawd$sunrise<-RiseNumHST
rawd$sunset<-SetNumHST



  ##-------------------- Create Light Data Frame ----------------------------------------------------------------------------

ldata <- data.frame(lon = rawd$lon, lat = rawd$lat, 
                    date = as.POSIXct(rawd$datetime, "%m/%d/%Y", tz="HST"), 
                    RiseMinute = rawd$sunrise, 
                    SetMinute = rawd$sunset)


## ---- sunrise sunset intervals  ----------------------------------------------------------

# Plotting the intervals of the rise and set minutes as relative to a time sequence gives a day to day profile of the light seen by the tag.
plot(RiseMinute ~ date, type = "n", ylim = range(c(RiseMinute, SetMinute)),  data = ldata, ylab = "tag day-time")
polygon(c(ldata$date, rev(ldata$date)), c(ldata$RiseMinute, rev(ldata$SetMinute)), col = "yellow")
title("daylight period (minute of the day in UTC)")

# Hampel Filter ldata$RiseMinute and ldata$SetMinute to get rid of spikes
HampelFilter <- function (x, k,t0=3){
  n <- length(x)
  y <- x
  ind <- c()
  L <- 1.4826
  for (i in (k + 1):(n - k)) {
    x0 <- median(x[(i - k):(i + k)])
    S0 <- L * median(abs(x[(i - k):(i + k)] - x0))
    if (abs(x[i] - x0) > t0 * S0) {
      y[i] <- x0
      ind <- c(ind, i)
    }
  }
  list(y = y, ind = ind)
}
# 
testrise<-HampelFilter(ldata$RiseMinute,5,t0=5) # smaller the t0, the smoother the filter
testset<-HampelFilter(ldata$SetMinute,5,t0=5)

plot(testrise$y ~ date, type = "n", ylim = range(c(testrise$y, testset$y)),  data = ldata, ylab = "tag day-time")
polygon(c(ldata$date, rev(ldata$date)), c(testrise$y, rev(testset$y)), col = "yellow")
title("HAMPLE TRAMPLED daylight period (minute of the day in UTC)")


##------ sometimes need a batch filter when sunset times are < sunrise times-------------

ix<-which(testset$y<testrise$y) # counterintuitive but b/c times are in GMT
points(ldata$date[ix],testset$y[ix],pch=21,col="black",bg="red")


newSet<-matrix(NA,length(ix),2)

for (j in 1:length(ix)) {
  
  fault=ix[j]
  
  if (fault<3){
    x00=fault
    x01=fault+10
  }else{
    unt<-ifelse(fault<11,2,15)
    x00=fault-unt
    x01=fault+unt
  }

  rng=c(x00:x01)
  
  plot(ldata$date[rng],testset$y[rng],pch=21,bg="black",ylim = range(c(testrise$y, testset$y))) 
  points(ldata$date[fault],testset$y[fault],pch=21,col="black",bg="yellow")
  
  # This new plot will show a range of sunrise or sunset times that need to be corrected. The erroneous point is in yellow. Click on a point that best estimates where the target point should be.
  p<-identify(ldata$date[rng],testset$y[rng],n=1)

  newSet[j,1]<-fault                                     
  newSet[j,2]<-testset$y[rng[p]]
  
  rm(p)
}


setixv<-as.numeric(newSet[,1])
chronx<-as.numeric(newSet[,2])

testset$y[setixv]<-chronx


plot(testrise$y~ date, type = "n", ylim = range(c(testrise$y, testset$y)),  data = ldata, ylab = "tag day-time")
polygon(c(ldata$date, rev(ldata$date)), c(testrise$y, rev(testset$y)), col = "yellow")
title("HAMPLE TRAMPLED AND MANUAL SWEEP daylight period (minute of the day in UTC)")

##---------------- sometimes needs lot of manual point selection to clean up light ---------------------

# Now do the final touches to the sunrise/sunset times
# Have to do sunrise and sunset separately (this is super clunky could that could be much improved)

coors2<-identify(ldata$date,testset$y)
# coors2<-identify(ldata$date,testrise$y)


# Manual Correction (again)
newSet2<-matrix(NA,length(coors2),2)
for (j in 1:length(coors2)) {
  fault=coors2[j]
  
  if (fault<3){
    x00=fault
    x01=fault+10
  }else{
    unt<-ifelse(fault<15,3,15)
    x00=fault-unt
    x01=fault+unt
  }  
  rng=c(x00:x01)
  
  #SET
  plot(ldata$date[rng],testset$y[rng],pch=21,bg="black",ylim = range(c(testrise$y, testset$y))) 
  points(ldata$date[fault],testset$y[fault],pch=21,col="black",bg="yellow")
  p<-identify(ldata$date[rng],testset$y[rng],n=1)
  newSet2[j,1]<-fault                                     
  newSet2[j,2]<-testset$y[rng[p]]
  
  #   RISE
#     plot(ldata$date[rng],testrise$y[rng],pch=21,bg="black",ylim = range(c(testrise$y, testset$y))) 
#     points(ldata$date[fault],testrise$y[fault],pch=21,col="black",bg="yellow")
#     p<-identify(ldata$date[rng],testrise$y[rng],n=1)
#     newSet2[j,1]<-fault                                     
#     newSet2[j,2]<-testrise$y[rng[p]]
#   
  rm(p)  
}


setixv<-as.numeric(newSet2[,1])
chronx<-as.numeric(newSet2[,2])

testset$y[setixv]<-chronx
#     testrise$y[setixv]<-chronx

plot(testset$y~ date, type = "n", ylim = range(c(testrise$y, testset$y)),  data = ldata, ylab = "tag day-time")  
polygon(c(ldata$date, rev(ldata$date)), c(testrise$y, rev(testset$y)), col = "yellow")
title("HAMPLE TRAMPLED AND MANUAL SWEEP daylight period (minute of the day in UTC)")

rm(newSet2)
rm(coors2)

##------------- Put filtered rise and set times back into original object

# first, fix times

rise.corr.str<-ifelse(nchar(as.character(testrise$y))<4, paste("0",as.character(testrise$y),sep=""), as.character(testrise$y))
set.corr.str<-ifelse(nchar(as.character(testset$y))<4, paste("0",as.character(testset$y),sep=""), as.character(testset$y))
# hst time string
rise.hms.str<-matrix(NA,length(rise.corr.str),1)
for (k in 1:length(rise.corr.str)){
  xispl<-substring(rise.corr.str[k], seq(1, nchar(rise.corr.str[k])-1, 2), seq(2, nchar(rise.corr.str[k]), 2))
  rise.hms.str[k]<-paste(xispl[1],":",xispl[2],":00",sep="")
}

set.hms.str<-matrix(NA,length(set.corr.str),1)
for (k in 1:length(set.corr.str)){
  xispl<-substring(set.corr.str[k], seq(1, nchar(set.corr.str[k])-1, 2), seq(2, nchar(set.corr.str[k]), 2))
  set.hms.str[k]<-paste(xispl[1],":",xispl[2],":00",sep="")
}



ldata$RiseMinute<-as.POSIXct(paste(ldata$date,rise.hms.str),"%Y-%m-%d %H:%M:%S",tz="HST")
ldata$SetMinute<-as.POSIXct(paste(ldata$date, set.hms.str),"%Y-%m-%d %H:%M:%S",tz="HST")


##------ save ldata ------------------------------------------------
dp2<-"/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/ldata/"
save(file=paste(dp2,'ldata_',filenames[i],sep=""),ldata)


rm(list=ls()[! ls() %in% c("filenames")])




##--------------------------------------------------------------------------------