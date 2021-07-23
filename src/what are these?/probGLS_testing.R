install.packages("remotes")
remotes::install_github("benjamin-merkel/probGLS", force= TRUE)
library(GeoLight)
library(raster)
library(ncdf4)
library(stringr)
library(probGLS)
library(dplyr)

######

library(SGAT)
library(BAStag)
library(chron)
library(pracma)
library(maptools)
library(rgeos)
library(maptools)
library(zoo)

#=================================================================================

# lotek_to_dataframe to convert your calculated riseset into a format usable by probGLS
# necessary formats are as follows:
 # date	column "TimeS" in daylog file as as.Date() [DONE]
 # sunrise column "Sunrise" in daylog file as as.character()
 # sunset column "Sunset" in daylog file as as.character()
 # lon/lat is in numeric [DONE]


###
#bc_bird110801 <- read.csv("toy.csv", stringsAsFactors = F)
#Fill empty cells with NA
#bc_bird110801 <- as.data.frame(apply(bc_bird110801, 2, function(x) gsub("^$|^ $", NA, x)))
# Filter out rows containing all NAs
#bc_bird110801 <- bc_bird110801 %>% filter_all(any_vars(!is.na(.)))


#######################################################################################################################

#Removing NAs

bc_bird110801 <- bc_bird110801[!is.na(bc_bird110801$sunrise), ]
#Fill empty cells with NA
bc_bird110801 <- as.data.frame(apply(bc_bird110801, 2, function(x) gsub("^$|^ $", NA, x)))
# Filter out rows containing all NAs
bc_bird110801 <- bc_bird110801 %>% filter_all(any_vars(!is.na(.)))

#######################################################################################################################

#Formatting for lotek_to_dataframe

sunrise <- as.character(sapply(strsplit(as.character(bc_bird110801$sunrise), split = " "), "[", 2))
sunrise.1 <- sapply(strsplit(sunrise, split = ":"), "[", 1)
sunrise.2 <- sapply(strsplit(sunrise, split = ":"), "[", 2)
sunrise <- paste(sunrise.1, sunrise.2, sep=":")

sunset <- as.character(sapply(strsplit(as.character(bc_bird110801$sunset), split = " "), "[", 2))
sunset.1 <- sapply(strsplit(sunset, split = ":"), "[", 1)
sunset.2 <- sapply(strsplit(sunset, split = ":"), "[", 2)
sunset <- paste(sunset.1, sunset.2, sep=":")

date <- bc_bird110801$date
TRLat <- bc_bird110801$lat
TRLon <- bc_bird110801$lon
test <- lotek_to_dataframe(date = date, sunrise = sunrise, sunset = sunset,
                           TRLon = TRLon, TRLat = TRLat)


start <- as.POSIXct("2007-12-28 00:00", tz="Pacific/Midway")
end   <- as.POSIXct("2008-12-04 08:55", tz="Pacific/Midway")
pr   <- prob_algorithm(trn                         = test, 
                       sensor                      = NULL, 
                       act                         = NULL, 
                       tagging.date                = start, 
                       retrieval.date              = end, 
                       loess.quartile              = NULL, 
                       tagging.location            = c(-177.35,28.2), 
                       particle.number             = 1000, 
                       iteration.number            = 100,
                       sunrise.sd                  = tw,
                       sunset.sd                   = tw,
                       range.solar                 = c(-7,-1),
                       boundary.box                = c(-220,-140,-90,0),
                       days.around.spring.equinox  = c(0,0), 
                       days.around.fall.equinox    = c(0,0),
                       speed.dry                   = c(12,6,45),
                       speed.wet                   = c(1,1.3,5), 
                       sst.sd                      = 0.5,       
                       max.sst.diff                = 3,          
                       east.west.comp              = T,
                       land.mask                   = T, 
                       ice.conc.cutoff             = 1, 
                       wetdry.resolution           = 1,
                       NOAA.OI.location            = "folder with environmental data and land mask")

xy <- coordinates(pr$`most probable track`)
xy

