# Calculate centroid and distance/angle to colony of polygons
# calculate lat/lon ranges
# Dec 3 2021
# Dallas Jordan


# Setup -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(geosphere)
library(sf)

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

mid_data <- rbind(lm,bm)
tern_data <-rbind(lt,bt)

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

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/final_ud/")
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


