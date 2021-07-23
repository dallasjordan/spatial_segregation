library(maptools)
library(rgdal)
library(GeoLocTools)
setupGeolocation()
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)


setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/2008")
convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}

files <- list.files(getwd())
test <- read.csv(files[1], 
                    stringsAsFactors = FALSE) 

colnames(test)= c("id","dtime","x","y")
test <- test[!is.na(test$x) & !is.na(test$y),]

test.sp <- test[,c("id","x","y")]
coordinates(test.sp) <- c("x","y")
proj4string(test.sp) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
# alternate method - SpatialPointsDataFrame()

test.sp

test.sp <- spTransform(test.sp,CRS("+proj=utm +zone=2 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
#reference bandwidth 
href <- (0.5*(sd(test.sp@coords[1:(nrow(test.sp@coords))])+sd(test.sp@coords[(nrow(test.sp@coords)+1):(nrow(test.sp@coords)*2)])))*(nrow(test.sp@coords)^(-1/6))


kernel.ref <- kernelUD(test.sp, h=href)
image(kernel.ref)

test.kernel.poly <- getverticeshr(kernel.ref, percent = 95, unin = "m", unout = "km2") 
print(test.kernel.poly)

color <- rep("white", nrow(test.sp@data))
data(world2HiresMapEnv)
plot(test$x,test$y ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,240), ylim = c(-10, 80), bty = "n")
map('world2Hires', xlim= c(120,240), ylim = c(-10, 75),  col = "grey", add = T)

test.kernel.poly <- spTransform(test.kernel.poly,CRS("+proj=longlat +zone=2 +ellps=WGS84 +datum=WGS84 +no_defs"))
test.sp <- spTransform(test.sp,CRS("+proj=longlat +zone=2 +ellps=WGS84 +datum=WGS84 +no_defs"))
test.kernel.poly@polygons[[1]]@Polygons[[1]]@coords <- convlon(test.kernel.poly@polygons[[1]]@Polygons[[1]]@coords)
test.sp@coords <- convlon(test.sp@coords)

sp::plot(test.kernel.poly, col = test.kernel.poly@data$id, add=T)
plot(test.sp, col = "red", pch = 16, cex=0.5, add=T)



#######################################
# Popular, alternate method of producing h - least square cross validation
test <- read.csv(files[1], 
                 stringsAsFactors = FALSE) 
colnames(test)= c("id","dtime","x","y")
test <- test[!is.na(test$x) & !is.na(test$y),]
test.sp <- test[,c("id","x","y")]
coordinates(test.sp) <- c("x","y")
proj4string(test.sp) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
test.sp <- spTransform(test.sp,CRS("+proj=utm +zone=2 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

kernel.LSCV<- kernelUD(test.sp, h="LSCV")
plotLSCV(kernel.LSCV)
image(kernel.LSCV)


test.kernel.poly <- getverticeshr(kernel.LSCV, percent = 95, unin = "m", unout = "km2") 
print(test.kernel.poly)

color <- rep("white", nrow(test.sp@data))
color[(test.sp@data$id == "1108_02")] <- "red"
color[(test.sp@data$id == "1108_03")] <- "green"
color[(test.sp@data$id == "1108_04")] <- "blue"
color[(test.sp@data$id == "1108_05")] <- "cyan"

data(world2HiresMapEnv)
plot(test$x,test$y ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,240), ylim = c(-10, 80), bty = "n")
map('world2Hires', xlim= c(120,240), ylim = c(-10, 75),  col = "grey", add = T)

test.kernel.poly <- spTransform(test.kernel.poly,CRS("+proj=longlat +zone=2 +ellps=WGS84 +datum=WGS84 +no_defs"))
test.sp <- spTransform(test.sp,CRS("+proj=longlat +zone=2 +ellps=WGS84 +datum=WGS84 +no_defs"))
test.kernel.poly@polygons[[1]]@Polygons[[1]]@coords <- convlon(test.kernel.poly@polygons[[1]]@Polygons[[1]]@coords)
test.sp@coords <- convlon(test.sp@coords)

sp::plot(test.kernel.poly, col = test.kernel.poly@data$id, add=T)
plot(test.sp, col = "red", pch = 16, cex=0.5, add=T)

# No dip/convergence is common when using the LSCV method of choosing h. 
# In addition, when we look at the heat map of this kernel (above), the 
# home ranges are heavily fragmented into many “islands.” In cases with 
# infrequent relocation data (eg. every few days or less), I would not 
# recommend this method of choosing h. If you have GPS collar data that 
# collected locations very frequently, the LSCV method of selecting h may 
# be the most appropriate. 




