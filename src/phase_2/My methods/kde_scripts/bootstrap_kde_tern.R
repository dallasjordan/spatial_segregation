library(dplyr)
library(maptools)
library(rgdal)
library(GeoLocTools)
setupGeolocation()
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)


convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
years <- c("2008","2009","2010","2011","2012")



# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  colnames(data_loop)= c("dtime","x","y")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","x","y")]


# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  colnames(data_loop)= c("dtime","x","y")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","x","y")]


# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)

######### RESAMPLING WITHOUT REPLACEMENT

data.spL <- split(data.sp, data.sp$id)

resample_2008 <- rbind(data.spL[["LAAL_2008"]],data.spL[["BFAL_2008"]])
resample_2009 <- rbind(data.spL[["LAAL_2009"]],data.spL[["BFAL_2009"]])
resample_2010 <- rbind(data.spL[["LAAL_2010"]],data.spL[["BFAL_2010"]])
resample_2011 <- rbind(data.spL[["LAAL_2011"]],data.spL[["BFAL_2011"]])
resample_2012 <- rbind(data.spL[["LAAL_2012"]],data.spL[["BFAL_2012"]])

year_im_on <- 2012
test_stat <- 0.367115062446283 # manually putting this in because its easier then writing code to get it from other script
iter <- 1000

significance_tally <- 0
iter_tally <- 0
results_storage <- vector()
for (i in 1:iter){
  resample_this <- resample_2012
  resample_this$x <- sample(resample_this$x,length(resample_this$x), replace=F)
  resample_this$y <- sample(resample_this$y,length(resample_this$y), replace=F)
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  coordinates(resample_this) <- c("x","y")
  proj4string(resample_this) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
  resample_this <- spTransform(resample_this,CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  kernel.boot <- kernelUD(resample_this, same4all = T)
  
  results<-kerneloverlaphr(kernel.boot, method="UDOI", percent=90)
  results_storage[i] <- results[1,2]
  if (results[1,2]<test_stat){
    print(results)
    print("+1!")
    significance_tally <- significance_tally+1
  }
}

p_value <- significance_tally/iter
mean_value <- mean(results_storage)
sd_value <- sd(results_storage)

year_im_on
p_value
mean_value
sd_value

