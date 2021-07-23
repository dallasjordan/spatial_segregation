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
  print(files)
  print(nrow(data_loop))
  colnames(data_loop)= c("dtime","x","y","delete1","delete2")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","x","y")]


# Make BFAL data object
convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
years <- c("2008","2009","2010","2011","2012")
SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  print(files)
  print(nrow(data_loop))
  colnames(data_loop)= c("dtime","x","y","delete1","delete2")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","x","y")]


# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)
data.plot <- data.sp # for later
coordinates(data.sp) <- c("x","y")
proj4string(data.sp) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")

data.sp <- spTransform(data.sp,CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
#reference bandwidth 
href <- (0.5*(sd(data.sp@coords[1:(nrow(data.sp@coords))])+sd(data.sp@coords[(nrow(data.sp@coords)+1):(nrow(data.sp@coords)*2)])))*(nrow(data.sp@coords)^(-1/6))


kernel.ref <- kernelUD(data.sp, same4all = T)
image(kernel.ref)

data.kernel.poly <- getverticeshr(kernel.ref, percent = 90, unin = "m", unout = "km2") 
print(data.kernel.poly)

data(world2HiresMapEnv)
plot(data.plot$x,data.plot$y ,type = "p", ylab="Latitude", xlab="Longitude",
     xlim= c(120,240), ylim = c(0, 80), bty = "n", cex=0.1)
maps::map('world2Hires', xlim= c(120,260), ylim = c(0, 75),  col = "grey", add = T)

data.kernel.poly <- spTransform(data.kernel.poly,CRS("+proj=longlat +zone=2 +ellps=WGS84 +datum=WGS84 +no_defs"))
data.sp <- spTransform(data.sp,CRS("+proj=longlat +zone=2 +ellps=WGS84 +datum=WGS84 +no_defs"))
for (i in 1:length(data.kernel.poly@polygons)){
  for (j in 1:length(data.kernel.poly@polygons[[i]]@Polygons)){
    data.kernel.poly@polygons[[i]]@Polygons[[j]]@coords <- convlon(data.kernel.poly@polygons[[i]]@Polygons[[j]]@coords)
  }
}
data.sp@coords <- convlon(data.sp@coords)

sp::plot(data.kernel.poly, col = data.kernel.poly@data$id, add=T)
sp::plot(data.sp, col = "red", pch = 16, cex=0.1, add=T)

results<-kerneloverlaphr(kernel.ref, method="UDOI", percent=90)
View(results)
write.csv(results,"/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/LAALvsBFAL_Tern_by_year.csv")
