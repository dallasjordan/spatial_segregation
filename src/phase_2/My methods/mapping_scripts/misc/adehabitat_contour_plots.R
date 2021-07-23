x <- seq(-50, 150, by=1.) # resolution is the pixel size you desire 
y <- seq(-50, 150, by=1.)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)


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

calculate_sp_obj_extent <- function(sp_obj){
  x_range <- sp_obj@bbox[1,2] - sp_obj@bbox[1,1]
  y_range <- sp_obj@bbox[2,2] - sp_obj@bbox[2,1]
  x_range_km <- x_range/1000
  y_range_km <- y_range/1000
  print1<- paste0("the x-range in km is ",x_range_km)
  print2<- paste0("the y-range in km is ",y_range_km)
  grid_calc <- x_range_km/300
  print3<- paste0("for a grid size of 300km^2 use grid parameter ",grid_calc)
  print(print1)
  print(print2)
  print(print3)
  return(grid_calc)
}

pdc_mercator_proj<-"+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
laea <- "+proj=laea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# load aggregated data

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_midway_withTrackID.Rdata")
BFALmid <- BFAL
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_tern_withTrackID.Rdata")
BFALtern <- BFAL

BFALmid$id <- "BFALmid"
BFALtern$id <- "BFALtern"
both <- rbind(BFALmid,BFALtern)

##### Need to calculate a single grid for both:
both_1<-both[,1:3]
sp::coordinates(both_1) <- c("x", "y")

proj4string(both_1) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs") # placeholder
both_1 <- spTransform(both_1,laea)
grid_input <- calculate_sp_obj_extent(both_1)
results <- kernelUD(both_1, grid=grid_input,same4all=T,extent=0.15)
image(results)
vud_points<-getvolumeUD(results)
image(vud_points)

# 8. Get contour
levels <- c(50, 95)
list <- vector(mode="list", length = 2)

list[[1]] <- as.image.SpatialGridDataFrame(vud_points[[1]])
list[[2]] <- as.image.SpatialGridDataFrame(vud_points[[2]])


# 9. Plot
par(mfrow = c(2, 1))
image(vud_points[[1]])
contour(list[[1]], add=TRUE, levels=levels)
image(vud_points[[2]])
contour(list[[2]], add=TRUE, levels=levels)

# 10. Get vertices (Will be an Error)
vkde_points <- getverticeshr(kud_points, percent = 75,
                             unin = 'm', unout='m2')
plot(vkde_points)








