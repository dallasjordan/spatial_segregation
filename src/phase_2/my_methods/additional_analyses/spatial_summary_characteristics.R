# Calculate centroid and distance/angle to colony of polygons
# calculate lat/lon ranges
# Dec 3 2021
# Dallas Jordan


# Setup -------------------------------------------------------------------
install.packages("remotes")
remotes::install_github("RodrigoAgronomia/PAR")

library(dplyr)
library(ggplot2)
library(geosphere)
library(sf)
library(pracma)
library(lwgeom)

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

mid_data <- rbind(lm,bm)
tern_data <-rbind(lt,bt)

midway <- data.frame(longitude = -177.3761,
                     latitude = 28.2101)
tern <- data.frame(longitude = -166.284,
                   latitude = 23.870)

midway <- st_as_sf(midway, coords = c("longitude", "latitude"), 
                   crs = 4326, agr = "constant")
tern <- st_as_sf(tern, coords = c("longitude", "latitude"), 
                 crs = 4326, agr = "constant")


# by class

lm_mean_x <- mean(lm$x)
lm_mean_y <- mean(lm$y)
plot(lm$x,lm$y, cex=0.1)
points(lm_mean_x,lm_mean_y, cex=2, col="red")

lt_mean_x <- mean(lt$x)
lt_mean_y <- mean(lt$y)
plot(lt$x,lt$y, cex=0.1)
points(lt_mean_x,lt_mean_y, cex=2, col="red")

bm_mean_x <- mean(bm$x)
bm_mean_y <- mean(bm$y)
plot(bm$x,bm$y, cex=0.1)
points(bm_mean_x,bm_mean_y, cex=2, col="red")

bt_mean_x <- mean(bt$x)
bt_mean_y <- mean(bt$y)
plot(bt$x,bt$y, cex=0.1)
points(bt_mean_x,bt_mean_y, cex=2, col="red")

# stopping here because Scott recommended weighted center of polygons. So, need to load in polygons and calculate centroid of those

# Centroid by polygons ----------------------------------------------------

# load in KDE polygons

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/final_ud/")
path <- getwd()
file_list <- list.files(path=path)

for (i in 1:length(file_list)){
  load(file_list[i])
}

midLAAL.poly95 <- getverticeshr(midLAAL, percent = 95, unin = "m", unout = "km2") 
midLAAL.poly50 <- getverticeshr(midLAAL, percent = 50, unin = "m", unout = "km2") 
midBFAL.poly95 <- getverticeshr(midBFAL, percent = 95, unin = "m", unout = "km2") 
midBFAL.poly50 <- getverticeshr(midBFAL, percent = 50, unin = "m", unout = "km2") 
ternLAAL.poly95 <- getverticeshr(ternLAAL, percent = 95, unin = "m", unout = "km2") 
ternLAAL.poly50 <- getverticeshr(ternLAAL, percent = 50, unin = "m", unout = "km2") 
ternBFAL.poly95 <- getverticeshr(ternBFAL, percent = 95, unin = "m", unout = "km2") 
ternBFAL.poly50 <- getverticeshr(ternBFAL, percent = 50, unin = "m", unout = "km2") 

test.df <- as.data.frame(midLAAL.poly95)

st_area(midLAAL.poly95.sf)

midLAAL.poly95.sf <- midLAAL.poly95 %>% st_as_sf()
mlp95centroid <- st_centroid(midLAAL.poly95.sf)
plot(midLAAL.poly95.sf$geometry)
plot(mlp95centroid$geometry, add=T)
mlp95centroid <- mlp95centroid %>% st_transform(4326)
mlp95centroid
st_distance(mlp95centroid,midway)
p <- st_sfc(st_point(c(-177.3761,28.2101)), st_point(c(172.1167,48.32789)), crs=4326)
rads<- st_geod_azimuth(p)
angle <- rad2deg(rads)
angle



midLAAL.poly50.sf <- midLAAL.poly50 %>% st_as_sf()
mlp50centroid <- st_centroid(midLAAL.poly50.sf)
plot(midLAAL.poly50.sf$geometry)
plot(mlp50centroid$geometry, add=T)
mlp50centroid <- mlp50centroid %>% st_transform(4326)
mlp50centroid
st_distance(mlp50centroid,midway)
p <- st_sfc(st_point(c(-177.3761,28.2101)), st_point(c(166.0309,47.71621)), crs=4326)
rads<- st_geod_azimuth(p)
angle <- rad2deg(rads)
angle

midBFAL.poly95.sf <- midBFAL.poly95 %>% st_as_sf()
mbp95centroid <- st_centroid(midBFAL.poly95.sf)
plot(midBFAL.poly95.sf$geometry)
plot(mbp95centroid$geometry, add=T)
mbp95centroid <- mbp95centroid %>% st_transform(4326)
mbp95centroid
st_distance(mbp95centroid,midway)
p <- st_sfc(st_point(c(-177.3761,28.2101)), st_point(c(-177.309,45.83291)), crs=4326)
rads<- st_geod_azimuth(p)
angle <- rad2deg(rads)
angle

midBFAL.poly50.sf <- midBFAL.poly50 %>% st_as_sf()
mbp50centroid <- st_centroid(midBFAL.poly50.sf)
plot(midBFAL.poly50.sf$geometry)
plot(mbp50centroid$geometry, add=T)
mbp50centroid <- mbp50centroid %>% st_transform(4326)
mbp50centroid
st_distance(mbp50centroid,midway)
p <- st_sfc(st_point(c(-177.3761,28.2101)), st_point(c(-177.8553,49.81271)), crs=4326)
rads<- st_geod_azimuth(p)
angle <- rad2deg(rads)
angle

ternLAAL.poly95.sf <- ternLAAL.poly95 %>% st_as_sf()
tlp95centroid <- st_centroid(ternLAAL.poly95.sf)
plot(ternLAAL.poly95.sf$geometry)
plot(tlp95centroid$geometry, add=T)
tlp95centroid <- tlp95centroid %>% st_transform(4326)
tlp95centroid
st_distance(tlp95centroid,tern)
p <- st_sfc(st_point(c(-166.284,23.87)), st_point(c(-179.3864,44.23886)), crs=4326)
rads<- st_geod_azimuth(p)
angle <- rad2deg(rads)
angle

ternLAAL.poly50.sf <- ternLAAL.poly50 %>% st_as_sf()
tlp50centroid <- st_centroid(ternLAAL.poly50.sf)
plot(ternLAAL.poly50.sf$geometry)
plot(tlp50centroid$geometry, add=T)
tlp50centroid <- tlp50centroid %>% st_transform(4326)
tlp50centroid
st_distance(tlp50centroid,tern)
p <- st_sfc(st_point(c(-166.284,23.87)), st_point(c(-178.1476,48.05546)), crs=4326)
rads<- st_geod_azimuth(p)
angle <- rad2deg(rads)
angle

ternBFAL.poly95.sf <- ternBFAL.poly95 %>% st_as_sf()
tbp95centroid <- st_centroid(ternBFAL.poly95.sf)
plot(ternBFAL.poly95.sf$geometry)
plot(tbp95centroid$geometry, add=T)
tbp95centroid <- tbp95centroid %>% st_transform(4326)
tbp95centroid
st_distance(tbp95centroid,tern)
p <- st_sfc(st_point(c(-166.284,23.87)), st_point(c(-159.8199,43.8779)), crs=4326)
rads<- st_geod_azimuth(p)
angle <- rad2deg(rads)
angle

ternBFAL.poly50.sf <- ternBFAL.poly50 %>% st_as_sf()
tbp50centroid <- st_centroid(ternBFAL.poly50.sf)
plot(ternBFAL.poly50.sf$geometry)
plot(tbp50centroid$geometry, add=T)
tbp50centroid <- tbp50centroid %>% st_transform(4326)
tbp50centroid
st_distance(tbp50centroid,tern)
p <- st_sfc(st_point(c(-166.284,23.87)), st_point(c(-137.8018,47.4593)), crs=4326)
rads<- st_geod_azimuth(p)
angle <- rad2deg(rads)
angle

# Calculate lat/lon ranges ------------------------------------------------
#longitude
min(lm$x)
min(lt$x)
min(bm$x)
min(bt$x)

max(lm$x)
max(lt$x)
max(bm$x)
max(bt$x)

diff_lm_x <- max(lm$x)-min(lm$x)
diff_lt_x <- max(lt$x)-min(lt$x)
diff_bm_x <- max(bm$x)-min(bm$x)
diff_bt_x <- max(bt$x)-min(bt$x)

diff_lm_x
diff_lt_x
diff_bm_x
diff_bt_x

# latitude 
min(lm$y)
min(lt$y)
min(bm$y)
min(bt$y)

max(lm$y)
max(lt$y)
max(bm$y)
max(bt$y)

diff_lm_y <- max(lm$y)-min(lm$y)
diff_lt_y <- max(lt$y)-min(lt$y)
diff_bm_y <- max(bm$y)-min(bm$y)
diff_bt_y <- max(bt$y)-min(bt$y)

diff_lm_y
diff_lt_y
diff_bm_y
diff_bt_y


