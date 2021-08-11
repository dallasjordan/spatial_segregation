# Script for additional analyses, per Scott's suggestions on improvements to manuscript
  # KDE Area
  # Centroid
  # Boxplots of lat/lon by species and colony
  # Oceanographic analysis
# Dallas Jordan Aug 7 2021

# reference...not really for this script but helpful
# https://stackoverflow.com/questions/62488162/use-of-tilde-and-period-in-r

# Setup -------------------------------------------------------------------

sf::sf_use_s2(FALSE)

# Load in KDEs and lat/lon data

  # loading lat/lon files
    # for mac

    load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_midway_withTrackID.Rdata")
    LAALmid <- LAAL
    LAALmid$id <- paste0("lm",LAALmid$track)
    load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_tern_withTrackID.Rdata")
    LAALtern <- LAAL
    LAALtern$id <- paste0("lt",LAALtern$track)
    
    LAAL <- rbind(LAALmid, LAALtern)
    
    load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_midway_withTrackID.Rdata")
    BFALmid <- BFAL
    BFALmid$id <- paste0("bm",BFALmid$track)
    load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_tern_withTrackID.Rdata")
    BFALtern <- BFAL
    BFALtern$id <- paste0("bt",BFALtern$track)
    
    BFAL <- rbind(BFALmid, BFALtern)
    
    # comparisons: allLAAL v allBFAL
    #              ternLAAL v midwayLAAL
    #              ternBFAL v midwayBFAL
    #              ternLAAL v ternBFAL
    #              midwayLAAL v midwayBFAL
    all_data <- rbind(LAAL,BFAL)
    lm <- all_data[grep("lm", all_data$id), ]
    lt <- all_data[grep("lt", all_data$id), ]
    bm <- all_data[grep("bm", all_data$id), ]
    bt <- all_data[grep("bt", all_data$id), ]
    
    # for pc
    setwd("E:/project_data/spatial_segregation/data")
    
    load("LAALdata_midway_withTrackID.Rdata")
    LAALmid <- LAAL
    LAALmid$id <- paste0("lm",LAALmid$track)
    load("LAALdata_tern_withTrackID.Rdata")
    LAALtern <- LAAL
    LAALtern$id <- paste0("lt",LAALtern$track)
    
    LAAL <- rbind(LAALmid, LAALtern)
    
    load("BFALdata_midway_withTrackID.Rdata")
    BFALmid <- BFAL
    BFALmid$id <- paste0("bm",BFALmid$track)
    load("BFALdata_tern_withTrackID.Rdata")
    BFALtern <- BFAL
    BFALtern$id <- paste0("bt",BFALtern$track)
    
    BFAL <- rbind(BFALmid, BFALtern)
    
    # comparisons: allLAAL v allBFAL
    #              ternLAAL v midwayLAAL
    #              ternBFAL v midwayBFAL
    #              ternLAAL v ternBFAL
    #              midwayLAAL v midwayBFAL
    all_data <- rbind(LAAL,BFAL)
    lm <- all_data[grep("lm", all_data$id), ]
    lt <- all_data[grep("lt", all_data$id), ]
    bm <- all_data[grep("bm", all_data$id), ]
    bt <- all_data[grep("bt", all_data$id), ]
    
  # Load in KDE contours to calculate area
    
    # for pc - You can load in the UDs and make contours or load in the contours
    load("E:/project_data/spatial_segregation/data/final_ud/allLAAL.Rdata")
    load("E:/project_data/spatial_segregation/data/final_ud/allBFAL.Rdata")
    load("E:/project_data/spatial_segregation/data/final_ud/midLAAL.Rdata")
    load("E:/project_data/spatial_segregation/data/final_ud/ternLAAL.Rdata")
    load("E:/project_data/spatial_segregation/data/final_ud/midBFAL.Rdata")
    load("E:/project_data/spatial_segregation/data/final_ud/ternBFAL.Rdata")
    
    # for pc - loading in the contours
    setwd("E:/project_data/spatial_segregation/figures/individual/midLAAL/master_script_contours/")
    load("vert95_midLAAL.Rdata")
    load("vert50_midLAAL.Rdata")
    setwd("E:/project_data/spatial_segregation/figures/individual/midBFAL/master_script_contours/")
    load("vert95_midBFAL.Rdata")
    load("vert50_midBFAL.Rdata")
    setwd("E:/project_data/spatial_segregation/figures/individual/ternLAAL/master_script_contours/")
    load("vert95_ternLAAL.Rdata")
    load("vert50_ternLAAL.Rdata")
    setwd("E:/project_data/spatial_segregation/figures/individual/ternBFAL/master_script_contours/")
    load("vert95_ternBFAL.Rdata")
    load("vert50_ternBFAL.Rdata")
    setwd("E:/project_data/spatial_segregation/figures/allLAAL_allBFAL/master_script_contours/")
    load("vert95_allLAAL.Rdata")
    load("vert50_allLAAL.Rdata")
    load("vert95_allBFAL.Rdata")
    load("vert50_allBFAL.Rdata")


  # SOME EXPLANATION: 
    # TWO WAYS TO GET AREA -> 
      # 1. getverticeshr(estUD object) -> has $area column of the SPDF
      # 2. kernel.area(estUD object, percent= ___ )
    
    # Note that the home-range sizes returned by this function are slightly
    # different from the home-range size stored in the SpatialPolygonsDataFrame
    # returned by the function getverticeshr. Indeed, while the former measures
    # the area covered by the rasterized home range (area covered by the set of pixels
    # of the grid included in the home range), the latter measures the area of the
    # vector home range (with smoother contour). However, note that the difference
    # between the two estimates decrease as the resolution of the grid becomes finer
    
    # We can see that 
    
    # this is where vert95 came from, in master script: 
    plot(getverticeshr(midLAAL,percent=95))
    # see, same thing: 
    plot(vert95_midLAAL)
    
    # compare the two methods to get area: 
    kernel.area(midLAAL, percent=95)
    vert95_midLAAL$area
    
    # slightly different, but we are rolling with vert95_midLAAL$area because it is the area with the smooth contours, should be a more conservative estimate? 
    


# Calculate area ----------------------------------------------------------

  # redoing getvertices hr to get the units I want
    vert95_midLAAL <- getverticeshr(midLAAL,percent=95, unin = "m", unout = "km2")
    vert50_midLAAL <- getverticeshr(midLAAL,percent=50, unin = "m", unout = "km2")
    vert95_midBFAL <- getverticeshr(midBFAL,percent=95, unin = "m", unout = "km2")
    vert50_midBFAL <- getverticeshr(midBFAL,percent=50, unin = "m", unout = "km2")
    vert95_ternLAAL <- getverticeshr(ternLAAL,percent=95, unin = "m", unout = "km2")
    vert50_ternLAAL <- getverticeshr(ternLAAL,percent=50, unin = "m", unout = "km2")
    vert95_ternBFAL <- getverticeshr(ternBFAL,percent=95, unin = "m", unout = "km2")
    vert50_ternBFAL <- getverticeshr(ternBFAL,percent=50, unin = "m", unout = "km2")
    vert95_allLAAL <- getverticeshr(allLAAL,percent=95, unin = "m", unout = "km2")
    vert50_allLAAL <- getverticeshr(allLAAL,percent=50, unin = "m", unout = "km2")
    vert95_allBFAL <- getverticeshr(allBFAL,percent=95, unin = "m", unout = "km2")
    vert50_allBFAL <- getverticeshr(allBFAL,percent=50, unin = "m", unout = "km2")

    # just copying over into Excel...
    vert95_midLAAL$area
    vert50_midLAAL$area
    vert95_midBFAL$area
    vert50_midBFAL$area
    vert95_ternLAAL$area
    vert50_ternLAAL$area
    vert95_ternBFAL$area
    vert50_ternBFAL$area
    vert95_allLAAL$area
    vert50_allLAAL$area
    vert95_allBFAL$area
    vert50_allBFAL$area
    

# Calculate centroid ------------------------------------------------------

library(geosphere)
# use function centroid()
    
    # Compute the centroid of longitude/latitude polygons. Unlike other functions 
    # in this package, there is no spherical trigonomery involved in the implementation 
    # of this function. Instead, the function projects the polygon to the (conformal) 
    # Mercator coordinate reference system, computes the centroid, and then inversely 
    # projects it to longitude and latitude. This approach fails for polygons that 
    # include one of the poles (and is rather biased for anything close to the poles).
    # The function should work for polygons that cross the -180/180 meridian (date line).

# allLAAL 
calc_allLAAL <- LAAL[,c(2,3)]
c_allLAAL<-centroid(calc_allLAAL)
plot(calc_allLAAL, pch = 16, cex = 0.05, col="grey")
points(x=c_allLAAL[1], y=c_allLAAL[2], pch = 18, cex=4, col= "red")
c_allLAAL

# allBFAL
calc_allBFAL <- BFAL[,c(2,3)]
c_allBFAL<-centroid(calc_allBFAL)
plot(calc_allBFAL, pch = 16, cex = 0.05, col="grey")
points(x=c_allBFAL[1], y=c_allBFAL[2], pch = 18, cex=4, col= "red")
c_allBFAL

      # testing because this looks funky
      test2 <- cbind(BFAL$x, BFAL$y)
      test2[,1] <- 1
      plot(test2)
      mean(test2[,2])
      points(1,46.60866, pch = 18, cex = 4, col = "red")
      
# midLAAL
calc_midLAAL <- lm[,c(2,3)]
c_midLAAL<-centroid(calc_midLAAL)
plot(calc_midLAAL, pch = 16, cex = 0.05, col="grey")
points(x=c_midLAAL[1], y=c_midLAAL[2], pch = 18, cex=4, col= "red")
c_midLAAL

# midBFAL
calc_midBFAL <- bm[,c(2,3)]
c_midBFAL<-centroid(calc_midBFAL)
plot(calc_midBFAL, pch = 16, cex = 0.05, col="grey")
points(x=c_midBFAL[1], y=c_midBFAL[2], pch = 18, cex=4, col= "red")
c_midBFAL

# ternLAAL
calc_ternLAAL <- lt[,c(2,3)]
c_ternLAAL<-centroid(calc_ternLAAL)
plot(calc_ternLAAL, pch = 16, cex = 0.05, col="grey")
points(x=c_ternLAAL[1], y=c_ternLAAL[2], pch = 18, cex=4, col= "red")
c_ternLAAL

# ternBFAL
calc_ternBFAL <- bt[,c(2,3)]
c_ternBFAL<-centroid(calc_ternBFAL)
plot(calc_ternBFAL, pch = 16, cex = 0.05, col="grey")
points(x=c_ternBFAL[1], y=c_ternBFAL[2], pch = 18, cex=4, col= "red")
c_ternBFAL


# Calculate centroid based on KDE/weighted means --------------------------
  
# https://stackoverflow.com/questions/23613655/calculating-weighted-polygon-centroids-in-r
library(plyr)
library(raster)

# load in KDEs

setwd("E:/project_data/spatial_segregation/data/final_ud")
path <- getwd()
file_list <- list.files(path=path)

for (i in 1:length(file_list)){
  load(file_list[i])
}

#midLAAL
midLAAL_rast <- raster(as(midLAAL, "SpatialPixelsDataFrame"))
plot(midLAAL_rast)
cell_num <- seq(1:893)
ud_value <- midLAAL_rast@data@values
cell_value <- cbind(cell_num, ud_value)
xy_of_cells <- xyFromCell(midLAAL_rast,cell_num)
cell_value_xy <- cbind(cell_value, xy_of_cells)

weighted_x <- stats::weighted.mean(cell_value_xy[,3],cell_value_xy[,2])
weighted_y <- stats::weighted.mean(cell_value_xy[,4],cell_value_xy[,2])

proj <- data.frame(x=weighted_x,y=weighted_y) 
coordinates(proj) <- ~x+y 
class(proj)
proj4string(proj) <- CRS("+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 
ll <- spTransform(proj,CRS("+proj=longlat +datum=WGS84"))

print(ll)
c_midLAAL <- c(ll@coords[1],ll@coords[2])

pr1 <- projectRaster(midLAAL_rast, crs="+proj=longlat +datum=WGS84")
plot(pr1)
points(x=ll@coords[1],y=ll@coords[2])




# midBFAL
midBFAL_rast <- raster(as(midBFAL, "SpatialPixelsDataFrame"))
plot(midBFAL_rast)
cell_num <- seq(1:893)
ud_value <- midBFAL_rast@data@values
cell_value <- cbind(cell_num, ud_value)
xy_of_cells <- xyFromCell(midBFAL_rast,cell_num)
cell_value_xy <- cbind(cell_value, xy_of_cells)

weighted_x <- stats::weighted.mean(cell_value_xy[,3],cell_value_xy[,2])
weighted_y <- stats::weighted.mean(cell_value_xy[,4],cell_value_xy[,2])

proj <- data.frame(x=weighted_x,y=weighted_y) 
coordinates(proj) <- ~x+y 
class(proj)
proj4string(proj) <- CRS("+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 
ll <- spTransform(proj,CRS("+proj=longlat +datum=WGS84"))

print(ll)
c_midBFAL <- c(ll@coords[1],ll@coords[2])

pr1 <- projectRaster(midBFAL_rast, crs="+proj=longlat +datum=WGS84")
plot(pr1)
points(x=ll@coords[1],y=ll@coords[2])





#ternLAAL
ternLAAL_rast <- raster(as(ternLAAL, "SpatialPixelsDataFrame"))
plot(ternLAAL_rast)
cell_num <- seq(1:893)
ud_value <- ternLAAL_rast@data@values
cell_value <- cbind(cell_num, ud_value)
xy_of_cells <- xyFromCell(ternLAAL_rast,cell_num)
cell_value_xy <- cbind(cell_value, xy_of_cells)

weighted_x <- stats::weighted.mean(cell_value_xy[,3],cell_value_xy[,2])
weighted_y <- stats::weighted.mean(cell_value_xy[,4],cell_value_xy[,2])

proj <- data.frame(x=weighted_x,y=weighted_y) 
coordinates(proj) <- ~x+y 
class(proj)
proj4string(proj) <- CRS("+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 
ll <- spTransform(proj,CRS("+proj=longlat +datum=WGS84"))

print(ll)
c_ternLAAL <- c(ll@coords[1],ll@coords[2])

pr1 <- projectRaster(ternLAAL_rast, crs="+proj=longlat +datum=WGS84")
plot(pr1)
points(x=ll@coords[1],y=ll@coords[2])




# ternBFAL
ternBFAL_rast <- raster(as(ternBFAL, "SpatialPixelsDataFrame"))
plot(ternBFAL_rast)
cell_num <- seq(1:893)
ud_value <- ternBFAL_rast@data@values
cell_value <- cbind(cell_num, ud_value)
xy_of_cells <- xyFromCell(ternBFAL_rast,cell_num)
cell_value_xy <- cbind(cell_value, xy_of_cells)

weighted_x <- stats::weighted.mean(cell_value_xy[,3],cell_value_xy[,2])
weighted_y <- stats::weighted.mean(cell_value_xy[,4],cell_value_xy[,2])

proj <- data.frame(x=weighted_x,y=weighted_y) 
coordinates(proj) <- ~x+y 
class(proj)
proj4string(proj) <- CRS("+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 
ll <- spTransform(proj,CRS("+proj=longlat +datum=WGS84"))

print(ll)
c_ternBFAL <- c(ll@coords[1],ll@coords[2])

print(ll)

pr1 <- projectRaster(ternBFAL_rast, crs="+proj=longlat +datum=WGS84")
plot(pr1)
points(x=ll@coords[1],y=ll@coords[2])




#allLAAL
allLAAL_rast <- raster(as(allLAAL, "SpatialPixelsDataFrame"))
plot(allLAAL_rast)
cell_num <- seq(1:893)
ud_value <- allLAAL_rast@data@values
cell_value <- cbind(cell_num, ud_value)
xy_of_cells <- xyFromCell(allLAAL_rast,cell_num)
cell_value_xy <- cbind(cell_value, xy_of_cells)

weighted_x <- stats::weighted.mean(cell_value_xy[,3],cell_value_xy[,2])
weighted_y <- stats::weighted.mean(cell_value_xy[,4],cell_value_xy[,2])

proj <- data.frame(x=weighted_x,y=weighted_y) 
coordinates(proj) <- ~x+y 
class(proj)
proj4string(proj) <- CRS("+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 
ll <- spTransform(proj,CRS("+proj=longlat +datum=WGS84"))

print(ll)
c_allLAAL <- c(ll@coords[1],ll@coords[2])

pr1 <- projectRaster(allLAAL_rast, crs="+proj=longlat +datum=WGS84")
plot(pr1)
points(x=ll@coords[1],y=ll@coords[2])




# allBFAL
allBFAL_rast <- raster(as(allBFAL, "SpatialPixelsDataFrame"))
plot(allBFAL_rast)
cell_num <- seq(1:893)
ud_value <- allBFAL_rast@data@values
cell_value <- cbind(cell_num, ud_value)
xy_of_cells <- xyFromCell(allBFAL_rast,cell_num)
cell_value_xy <- cbind(cell_value, xy_of_cells)

weighted_x <- stats::weighted.mean(cell_value_xy[,3],cell_value_xy[,2])
weighted_y <- stats::weighted.mean(cell_value_xy[,4],cell_value_xy[,2])

proj <- data.frame(x=weighted_x,y=weighted_y) 
coordinates(proj) <- ~x+y 
class(proj)
proj4string(proj) <- CRS("+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 
ll <- spTransform(proj,CRS("+proj=longlat +datum=WGS84"))

print(ll)
c_allBFAL <- c(ll@coords[1],ll@coords[2])

pr1 <- projectRaster(allBFAL_rast, crs="+proj=longlat +datum=WGS84")
plot(pr1)
points(x=ll@coords[1],y=ll@coords[2])

# distance and bearing to colonies
  # shortest path on ellipsoid (geodesic)
library(geosphere)

# midway <- c(463347.57, 3120213.67)
# tern <- c(369236.99, 2640401.49)

midway <- c(182.63, 28.2)
tern <- c(193.716, 23.87)

distGeo(midway, c_midLAAL)/1000
bearing(midway, c_midLAAL) # output is in degrees

distGeo(midway, c_midBFAL)/1000
bearing(midway, c_midBFAL) # output is in degrees

distGeo(tern, c_ternLAAL)/1000
bearing(tern, c_ternLAAL) # output is in degrees

distGeo(tern, c_ternBFAL)/1000
bearing(tern, c_ternBFAL) # output is in degrees




# Lat/lon boxplot ---------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(gridExtra)
library(wesanderson)
library(lattice)
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggpubr)


lm_box <- lm
lm_box$id <- "midLAAL"

bm_box <- bm
bm_box$id <- "midBFAL"

lt_box <- lt
lt_box$id <- "ternLAAL"

bt_box <- bt
bt_box$id <- "ternBFAL"

al_box <- rbind(lm, lt)
al_box$id <- "allLAAL"

ab_box <- rbind(bm,bt)
ab_box$id <- "allBFAL"

all_box <- rbind(lm_box, bm_box, lt_box, bt_box, al_box, ab_box)
all_box <- all_box[,1:3]

pal <- wes_palette("Zissou1", 4, type="discrete")
lon_box_lm <- ggplot(all_box, aes(id, x, fill=factor(id))) +
  geom_boxplot() +
  theme_classic() +
  labs(title ="Longitude by class", x="class", y="Longitude") + 
  theme(legend.position = "none")+
  scale_fill_manual(name="class", values=c(pal[1],pal[2], pal[3], pal[4], "wheat1", "wheat3"))
lon_box_lm

lat_box_lm <- ggplot(all_box, aes(id, y, fill=factor(id))) +
  geom_boxplot() +
  theme_classic() +
  labs(title ="Latitude by class", x="class", y="Latitude") + 
  theme(legend.position = "none")+
  scale_fill_manual(name="class", values=c(pal[1],pal[2], pal[3], pal[4], "wheat1", "wheat3"))
lat_box_lm
