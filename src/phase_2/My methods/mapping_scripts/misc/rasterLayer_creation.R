###### This script generates rasters from the KDEs made in my KDE_ scripts
# March 18 2021
# Not sure if this script is useful anymore, June 28 2021. Workflow has been adjusted to generate rasters in master_script, then plotted in overallLAALBFAL_mapping

########## EXAMPLE USED FOR TESTING

if (!require("rspatial")) devtools::install_github('rspatial/rspatial')
install.packages('rspatial')
library(rspatial)
city <- sp_data('city')
crime <- sp_data('crime')

plot(city, col='light blue')
points(crime, col='red', cex=.5, pch='+')

CityArea <- raster::area(city)
dens <- nrow(xy) / CityArea
r <- raster(city)
res(r) <- 1000
r

r <- rasterize(city, r)
plot(r)
quads <- as(r, 'SpatialPolygons')
plot(quads, add=TRUE)
points(crime, col='red', cex=.5)

nc <- rasterize(coordinates(crime), r, fun='count', background=0)
plot(nc)
plot(city, add=TRUE)


r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = "+proj=laea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

xlim <- range(LAAL[,2])
ylim <- range(LAAL[,3])


#####
# ACTUAL START 

library(sf)
library(sp)
library(rspatial)
library(raster)
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures")
load('BFALdata_tern.Rdata')
data_tern<-BFAL

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures")
load('BFALdata_midway.Rdata')
data_midway<-BFAL

data <- rbind(data_tern,data_midway)# load point data
BFAL <- data
BFAL <- BFAL[,2:3]

load("BFAL_combined_KDE.Rdata")
r <- raster(fhat)
crs(r) <- "+proj=laea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
# resolution is already what you need it to be because you set that when you created fhat
r
raster::plot(r)

# now take an spdf of the LAAL point and use r <- rasterize(spdf, r)
# spdf <- st_as_sf(LAAL, coords=c("x","y"), crs = "+proj=longlat +datum=WGS84")
# spdf <- st_transform(spdf, "+proj=laea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
# st_drop_geometry(spdf)


coordinates(BFAL) <- c("x","y")
proj4string(BFAL) <- "+proj=longlat +datum=WGS84"
BFAL <- spTransform(BFAL, CRS("+proj=laea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
r <- rasterize(BFAL,r, fun='count', background=0)
raster::plot(r)
quads <- as(r, 'SpatialPolygons')
save(r, file="BFAL_combined_pointdens_raster.Rdata")

r_df <- as.data.frame(r, xy = TRUE) 

str(r_df)
# plot(quads, add=TRUE)


############## export raster for Arc
extent(velocity) <- extent(r) # assign it the same extent as the ith sst raster "x"
projection(velocity) <- CRS("+proj=longlat +datum=WGS84")
writeRaster(r, filename = "Rast_combined_BFAL.tif", options=c('TFW=YES'), overwrite=TRUE)

############## export spPolyDF for Arc
rgdal::writeOGR(spPolyDF, shapefilewritepath, "pws_KDE" ,driver="ESRI Shapefile", verbose=TRUE, overwrite=TRUE)



