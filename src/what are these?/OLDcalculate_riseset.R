# Estimating sunrise/sunset from estimated lat/longs on Lotek tags

############################################################################################################

install.packages("suncalc")
install.packages("remotes")
remotes::install_github("benjamin-merkel/probGLS", force= TRUE)
#library(GeoLight)
#library(SGAT)
library(chron)
library(dplyr)
library(suncalc)

setwd("/Users/dallasjordan/Desktop/Stony Brook/SoMAS/Thesis/R/Spatial segregation/Spatial segregation")
filenames<- list.files()

############################################################################################################

#Load in trimmed data set for first tag, will have to repeat manually for each tag. 

bc_bird110801 <- read.csv("1108_01.csv", stringsAsFactors = F)
bc_bird110801 <- select(bc_bird110801, -Bird.ID, -SST)
bc_bird110801 <- bc_bird110801[,c(1,3,2)]
bc_bird110801$Date.and.Time <- strptime(bc_bird110801$Date.and.Time, format="%m/%d/%y", tz="Pacific/Midway")
bc_bird110801$Date.and.Time <- as.Date(bc_bird110801$Date.and.Time, format= "%F", tz="Pacific/Midway")
#bc_bird110801$Date.and.Time <- format(bc_bird110801$Date.and.Time, "%F")
bc_bird110801 <- bc_bird110801 %>% rename(date = Date.and.Time,lat = Raw.Latitude,lon = Raw.Longitude)
rise_set <- getSunlightTimes(data=bc_bird110801, keep= c("sunrise", "sunset","nauticalDawn"), tz="Pacific/Midway")
bc_bird110801$sunrise <-rise_set[,4]
bc_bird110801$sunset <- rise_set[,5]
#bc_bird110801$sunrise <-strptime(bc_bird110801$sunrise, format="%H:%M", tz="Pacific/Midway")

############################################################################################################

#Test probGLS. Should have required arguments for lotek_to_dataframe, which you save as an object, then use as an argument in prob_algorithm

#bc_bird110801$Date.and.Time <- as.Date(as.character(bc_bird110801$Date.and.Time), format = c("%m/%d/%y %H:%M"))
#bc_bird110801$Raw.Longitude <- as.numeric(as.character(bc_bird110801$Raw.Longitude))
#bc_bird110801$Raw.Latitude <- as.numeric(as.character(bc_bird110801$Raw.Latitude))
#bc_bird110801$Sunrise <- as.character(bc_bird110801$Sunrise)
#bc_bird110801$Sunset <- as.character(bc_bird110801$Sunset)
#test <- lotek_to_dataframe(bc_bird110801$Date.and.Time,bc_bird110801$Sunrise, bc_bird110801$Sunset, bc_bird110801$Raw.Longitude, bc_bird110801$Raw.Latitude)
#test


############################################################################################################

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

#Set:
sst.set<-utcsunset-chron(times="10:00:00")
s.hrs<-hours(hst.set)
s.min<-minutes(hst.set)
hourminset<-sprintf("%02d%02d", s.hrs, s.min) 
SetNumHST<-as.numeric(hourminset)