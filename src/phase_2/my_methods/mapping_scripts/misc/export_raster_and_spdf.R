# Export rasters and SpatialPolygonsDataFrames for plotting in Arcmap
# as of April 1 2021, need to modify line 16 (as(results$___))
# This is the script you use to get rasters/contours in Arc! as of April 9 2021
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
library(sf)

names(midLAAL_v_ternLAAL)
r<- raster(as(midLAAL_v_ternLAAL$ternLAAL,"SpatialPixelsDataFrame"))
writeRaster(r, filename = "Rast_LAAL_tern_normalized.tif", options=c('TFW=YES'), overwrite=TRUE)


x<-getverticeshr(midLAAL_v_ternLAAL$midLAAL, percent=95)
y<-getverticeshr(midLAAL_v_ternLAAL$ternLAAL, percent=95)
x1<-getverticeshr(midLAAL_v_ternLAAL$midLAAL, percent=50)
y1<-getverticeshr(midLAAL_v_ternLAAL$ternLAAL, percent=50)
rgdal::writeOGR(y1, getwd(), "LAALtern50normalized" ,driver="ESRI Shapefile", verbose=TRUE, overwrite=TRUE) 
