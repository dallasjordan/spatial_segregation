library(SGAT)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# for LAT - make sure computer is in HST ..... this script only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
setwd("C:/Users/Melinda/Dropbox/Chapter_03_Finalize/data/Tern_GLS_LOTEK_allrawfiles_Rdata/DayLogs")
filenames<- list.files()

load(filenames[i])
filenames[i]


deploy.time <- as.POSIXct("2010-01-6 00:00:00", tz = "HST")
recov.date <- ISOdatetime(2010, 9,12, 0, 0, 0, tz = "HST")

endix<-which.min(abs(recov.date-as.POSIXct(rawd$datetime,"%m/%d/%Y",tz="HST")))
rawd <- rawd[27:endix,]
rawd <- rawd[as.POSIXct(as.character(rawd$date),"%m/%d/%Y",tz="HST") >= deploy.time, ]


rawd$lat<-rawd$lat*.01
rawd$lon<-rawd$lon*.01
# 

##--------- record some tag details ------------------
lon.home <- -166.29 + 360
lat.home <- 23.87
deploy.time <- as.POSIXct("2009-03-11 00:00:00", tz = "GMT")

## maximum range in lon/lat (for land mask plus maximum bound of track)
lon.range <- c(140, 235)
lat.range <- c(20, 65)

rawd$lon[rawd$lon < 0] <- rawd$lon[rawd$lon < 0] + 360
plot(rawd$lon, rawd$lat)

##----------- extract the raw sunrise and sunset values----------

ldata <- data.frame(lon = rawd$lon, lat = rawd$lat, 
                    date = as.POSIXct(as.Date(rawd$datetime, "%m/%d/%Y")), 
                    RiseMinute = rawd$sunrise, 
                    SetMinute = rawd$sunset)

## out of bounds is 65535, also a 2527 value at 2009-03-08
##ldata <- ldata[pmax(ldata$SetMinute, ldata$RiseMinute) < 2500, ]


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
testrise<-HampelFilter(ldata$RiseMinute,5,t0=5)
testset<-HampelFilter(ldata$SetMinute,5,t0=5)

# testrise<-HampelFilter(testrise$y,5,t0=2)
# testset<-HampelFilter(testset$y,5,t0=2)

plot(testrise$y ~ date, type = "n", ylim = range(c(testrise$y, testset$y)),  data = ldata, ylab = "tag day-time")
polygon(c(ldata$date, rev(ldata$date)), c(testrise$y, rev(testset$y)), col = "yellow")
title("HAMPLE TRAMPLED daylight period (minute of the day in UTC)")


##------ sometimes need a batch filter when sunset times are < sunrise times-------------

# ix<-which(testset$y<testrise$y) # counterintuitive but b/c times are in GMT
# points(ldata$date[ix],testset$y[ix],pch=21,col="black",bg="red")
# 
# 
# newSet<-matrix(NA,length(ix),2)
# 
# for (j in 1:length(ix)) {
#   
#   fault=ix[j]
#   
#   if (fault<3){
#     x00=fault
#     x01=fault+10
#   }else{
#     unt<-ifelse(fault<11,2,15)
#     x00=fault-unt
#     x01=fault+unt
#   }
# 
#   rng=c(x00:x01)
#   
#   plot(ldata$date[rng],testset$y[rng],pch=21,bg="black",ylim = range(c(testrise$y, testset$y))) 
#   points(ldata$date[fault],testset$y[fault],pch=21,col="black",bg="yellow")
#   
#   p<-identify(ldata$date[rng],testset$y[rng],n=1)
# 
#   newSet[j,1]<-fault                                     
#   newSet[j,2]<-testset$y[rng[p]]
#   
#   rm(p)
# }


# setixv<-as.numeric(newSet[,1])
# chronx<-as.numeric(newSet[,2])
# 
# testset$y[setixv]<-chronx
# 
# 
# plot(testrise$y~ date, type = "n", ylim = range(c(testrise$y, testset$y)),  data = ldata, ylab = "tag day-time")
# polygon(c(ldata$date, rev(ldata$date)), c(testrise$y, rev(testset$y)), col = "yellow")
# title("HAMPLE TRAMPLED AND MANUAL SWEEP daylight period (minute of the day in UTC)")

##---------------- sometimes needs lot of manual point selection to clean up light ---------------------

# Now do the final touches to the sunrise/sunset times
# coors2<-identify(ldata$date,testset$y)
coors2<-identify(ldata$date,testrise$y)


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
  
#   #SET
#   plot(ldata$date[rng],testset$y[rng],pch=21,bg="black",ylim = range(c(testrise$y, testset$y))) 
#   points(ldata$date[fault],testset$y[fault],pch=21,col="black",bg="yellow")
#   p<-identify(ldata$date[rng],testset$y[rng],n=1)
#   newSet2[j,1]<-fault                                     
#   newSet2[j,2]<-testset$y[rng[p]]
# # #   
#     RISE
    plot(ldata$date[rng],testrise$y[rng],pch=21,bg="black",ylim = range(c(testrise$y, testset$y))) 
    points(ldata$date[fault],testrise$y[fault],pch=21,col="black",bg="yellow")
    p<-identify(ldata$date[rng],testrise$y[rng],n=1)
    newSet2[j,1]<-fault                                     
    newSet2[j,2]<-testrise$y[rng[p]]
  
  rm(p)  
}


setixv<-as.numeric(newSet2[,1])
chronx<-as.numeric(newSet2[,2])

    testset$y[setCixv]<-chronx
    testrise$y[setixv]<-chronx
# 
plot(testset$y~ date, type = "n", ylim = range(c(testrise$y, testset$y)),  data = ldata, ylab = "tag day-time")  
polygon(c(ldata$date, rev(ldata$date)), c(testrise$y, rev(testset$y)), col = "yellow")
title("HAMPLE TRAMPLED AND MANUAL SWEEP daylight period (minute of the day in UTC)")


rm(newSet2)
rm(coors2)

##------------- Put filtered rise and set times back into original object

# first, fix times

risetmp<-ldata$date + testrise$y * 60
settmp<-ldata$date + testset$y * 60

ldata$RiseMinute<-format(risetmp,tz="HST",usetz=TRUE)
ldata$SetMinute<-format(settmp,tz="HST",usetz=TRUE)


dp2<-"/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/ldata/"
save(file=paste(dp2,'ldata_',filenames[i],sep=""),ldata)

rm(list=ls()[! ls() %in% c("filenames")])



