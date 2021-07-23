Species <- "UriLom"
ID <- "2655"

lat.calib <- 78.17
lon.calib <- 15.15

# define deployement and retrieval dates
start    <- as.Date("2012-07-11")
end      <- as.Date("2013-05-01")

wd <- "data"

raw       <- read.csv("data/RawData/UriLom/2655.csv")
raw$TimeS <- as.Date(strptime(raw$TimeS,"%d/%m/%y"))
head(raw, 2)

raw2  <- raw[raw$TimeS > start & raw$TimeS < end,] 

data(wrld_simpl)
plot(raw2$TRLon,raw2$TRLat ,type = "n", ylab="Latitude", xlab="Longitude",
     xlim = c(-180, 180), ylim = c(-60, 90), bty = "n")
plot(wrld_simpl, col = "grey90", border = "grey40", add = T)
points(raw2$TRLon,raw2$TRLat, pch=16, col="cornflowerblue", type = "o")

trn <- lotek_to_dataframe(date    = raw$TimeS,
                          sunrise = as.character(raw$Sunrise),
                          sunset  = as.character(raw$Sunset),
                          TRLat   = raw$TRLat,
                          TRLon   = raw$TRLon)
trn <- trn[!is.na(trn$tSecond),]
head(trn,2)

trn$keep <- loessFilter(trn, k = 3, plot = T)

data <- read.csv("data/RawData/UriLom/2655_Basic Log.csv")

## datetime object needs to be in POSIXct, UTC time zone and must be called 'dtime'
data$dtime     <- as.POSIXct(strptime(data[,1], format = "%H:%M:%S %d/%m/%y"),tz='UTC')
data           <- data[!is.na(data$dtime),]
act            <- data
## wet dry data column must be called 'wetdry'
act$wetdry    <- 1-act$WetDryState
act <- subset(act, select = c(dtime,wetdry))

head(act)

speed_dry  = c(17, 4, 30) 

sd <- data.frame(speed=seq(0,35,0.1),prob=dnorm(seq(0,35,0.1),speed_dry[1],speed_dry[2]))
sd$prob <- sd$prob/max(sd$prob)
sd$prob[sd$speed<speed_dry[1]] <- 1
plot(sd,type="l",xlim=c(0,31),xlab="speed [m/s]",ylab="probability",lwd=2)
abline(v=speed_dry[3],lty=2,col="orange")

speed_wet  = c(1, 1.3, 5)

sd <- data.frame(speed=seq(0,35,0.1),prob=dnorm(seq(0,35,0.1),speed_wet[1],speed_wet[2]))
sd$prob <- sd$prob/max(sd$prob)
sd$prob[sd$speed<speed_wet[1]] <- 1
plot(sd,type="l",xlim=c(0,31),xlab="speed [m/s]",ylab="probability",lwd=2)
abline(v=speed_wet[3],lty=2,col="orange")

td <- data

## only keep temperature values when the logger was immersed in salt water 
## and if the 3 previous readings have been recorded while immersed in sea water as well
td$WetDryState.before   <- c(NA,head(td$WetDryState,-1))
td$WetDryState.before.2 <- c(NA,NA,head(td$WetDryState,-2))
td$WetDryState.before.3 <- c(NA,NA,NA,head(td$WetDryState,-3))
td           <- td[td$WetDryState         ==0 & 
                     td$WetDryState.before  ==0 & 
                     td$WetDryState.before.2==0 & 
                     td$WetDryState.before.3==0,]

## determine daily SST value recorded by the logger
sst           <- sst_deduction(datetime = td$dtime, temp = td$IntTemp, temp.range = c(-2,19))
abline(h=c(-2,19),lty=2,col="orange")

tw <- twilight_error_estimation(shape = 2.49, scale = 0.94, delay = 0)

pr   <- prob_algorithm(trn                         = trn,
                       sensor                      = sst[sst$SST.remove==F,], 
                       act                         = act, 
                       tagging.date                = start, 
                       retrieval.date              = end, 
                       loess.quartile              = NULL, 
                       tagging.location            = c(lon.calib, lat.calib), 
                       particle.number             = 500, 
                       iteration.number            = 100, 
                       sunrise.sd                  = tw,
                       sunset.sd                   = tw,
                       range.solar                 = c(-7,-1), 
                       speed.wet                   = speed_wet,
                       speed.dry                   = speed_dry, 
                       sst.sd                      = 0.5, 
                       max.sst.diff                = 3, 
                       boundary.box                = c(-90,120,40,90), 
                       days.around.spring.equinox  = c(21,14),   
                       days.around.fall.equinox    = c(14,21),
                       ice.conc.cutoff             = 0.9, 
                       land.mask                   = T,   # The track is assumed not to be on land
                       med.sea                     = F,   # if T the track cannot enter the Mediterranean Sea        
                       black.sea                   = F,   # if T the track cannot enter the Black Sea 
                       baltic.sea                  = F,   # if T the track cannot enter the Baltic Sea 
                       caspian.sea                 = F,   # if T the track cannot enter the Caspian Sea 
                       east.west.comp              = F,   # if true use the east west compensation (see Biotrack manual)
                       wetdry.resolution           = 300, # in seconds, i.e. 5 minutes = 300 seconds
                       NOAA.OI.location            = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/environment_data")


plot_map(pr)

# save geographic median .csv for all points in output object
most_probable <- pr[[2]]
most_probable_dtime <- as.data.frame(date(most_probable$dtime))
most_probable_lon <- as.data.frame(most_probable$lon)
most_probable_lat <- as.data.frame(most_probable$lat)
most_probable_export <- cbind(most_probable_dtime,most_probable_lon,most_probable_lat)
colnames(most_probable_export) <- c("dtime","Longitude", "Latitude")
# If you need to save this:
wd <- getwd()
write.csv(most_probable_export,file="/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/for_sarah.csv", row.names = FALSE) 
