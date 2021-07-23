library(maptools)
library(rgdal)
library(GeoLocTools)
setupGeolocation()
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)

# February 2021 tasks
  # redo KDE here with grid = 300,000 in the correct projection (EPSG 3832?) that uses m as unit
  # get the area of each for sample size analysis
  # on your call with Julia, she showed you how to measure the grid cell size, view the calculated 
  # grid cell size ('Source' tab under properties for density layer), and then how to change it manually. 
  # Need to re-run the tool in arcGIS in order to have the correct resolution

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
years <- c("2008","2009","2010","2011","2012")

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  colnames(data_loop)= c("dtime","x","y")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","x","y")]


# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  colnames(data_loop)= c("dtime","x","y")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","x","y")]


# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)
data.plot <- data.sp # for later
coordinates(data.sp) <- c("x","y")
proj4string(data.sp) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs") # placeholder
pdc_mercator_proj<-"+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
data.sp <- spTransform(data.sp,pdc_mercator_proj)

#reference bandwidth 
href <- (0.5*(sd(data.sp@coords[1:(nrow(data.sp@coords))])+sd(data.sp@coords[(nrow(data.sp@coords)+1):(nrow(data.sp@coords)*2)])))*(nrow(data.sp@coords)^(-1/6))

kernel.ref <- kernelUD(data.sp, same4all = T, extent = 1, grid=50)
image(kernel.ref)

# data.kernel.poly <- getverticeshr(kernel.ref, percent = 80, unin = "m", unout = "km2") 
# print(data.kernel.poly)
# 
# data(world2HiresMapEnv)
# plot(data.plot$x,data.plot$y ,type = "p", ylab="Latitude", xlab="Longitude",
#      xlim= c(120,240), ylim = c(-10, 80), bty = "n", cex=0.1)
# maps::map('world2Hires', xlim= c(120,240), ylim = c(-10, 75),  col = "grey", add = T)
# 
# data.kernel.poly <- spTransform(data.kernel.poly,CRS("+proj=longlat +zone=2 +ellps=WGS84 +datum=WGS84 +no_defs"))
# data.sp <- spTransform(data.sp,pdc_mercator_proj$proj4string)
# for (i in 1:length(data.kernel.poly@polygons)){
#   for (j in 1:length(data.kernel.poly@polygons[[i]]@Polygons)){
#     data.kernel.poly@polygons[[i]]@Polygons[[j]]@coords <- convlon(data.kernel.poly@polygons[[i]]@Polygons[[j]]@coords)
#   }
# }
# data.sp@coords <- convlon(data.sp@coords)
# 
# sp::plot(data.kernel.poly, col = data.kernel.poly@data$id, add=T)
# sp::plot(data.sp, col = "red", pch = 16, cex=0.1, add=T)

results<-kerneloverlaphr(kernel.ref, method="UDOI", percent=90)
results
View(results)
write.csv(results,"/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/LAALvsBFAL_Midway_by_year.csv")
