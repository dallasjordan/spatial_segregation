# Sharing code with Julia to calc KDEs and overlap
# Sept 29th, using code to prep before presenting at lab meeting on inital findings

# README JULIA:
  # Three steps and their line numbers: 
  #   1. Load in data and calculate grid extent to get 300km resolution of KDEs - Line 68
  #   2. Calculate KDEs for each individual track using the extent you calculated. You will have 3000 KDEs. - Line 84
  #   3. Average the KDEs for each "class" (e.g. Bird spp. A) to generate 1 final KDE - Line 96
  #   4. Create estUDm class objects for each comparison - Line 204
  #   5. Calculate test stats - Line 260

# This script is most up to date and changes randomization procedure to permute cell values rather than the labels
# This script: 
# 1. Loads post-processed data, converts to SpatialPointsDataFrame
# 2. Creates KDEs using the adehabitatHR package for each individual, saves as "results"
# 3. Does not normalize these for the relevant group (there are 5 total, LAAL/BFAL for each island and combined) - 
#    we are dealing with probability densities, don't need to be normalized
# 4. Averages each group to output a single raster. 
# 5. Assigns these estUD objects to coded names for overlap comparisons, e.g. lm_averaged_estUD = LAALmidway
# 6. Makes estUDm objects for each comparison (combines estUD objects pairwise according to the comparison of interest)
# 7. Calculates overlap of 95 UD, 50 UD, 95 PHR, 50 PHR, 95 BA, and 50 BA for each comparison of original data
# 8. Randomization procedure to test for significance

# Setup -------------------------------------------------------------------

# packages
library(dplyr)
library(maptools)
library(rgdal)
# library(GeoLocTools)
# setupGeolocation()
# Not available for R>=4.0
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)
library(sf)
library(ggplot2)

# functions and objects
calculate_sp_obj_extent <- function(sp_obj, extent_parameter){
  x_range <- sp_obj@bbox[1,2] - sp_obj@bbox[1,1]
  x_min <- sp_obj@bbox[1,1]-(extent_parameter*x_range)
  x_max <- sp_obj@bbox[1,2]+(extent_parameter*x_range)
  ade_x_range <- x_max-x_min
  y_range <- sp_obj@bbox[2,2] - sp_obj@bbox[2,1]
  y_min <- sp_obj@bbox[2,1]-(extent_parameter*y_range)
  y_max <- sp_obj@bbox[2,2]+(extent_parameter*y_range)
  ade_y_range <- y_max-y_min
  x_range_km <- ade_x_range/1000
  y_range_km <- ade_y_range/1000
  print1<- paste0("the x-range in km is ",x_range_km)
  print2<- paste0("the y-range in km is ",y_range_km)
  grid_calc <- x_range_km/300
  print3<- paste0("for a grid size of 300km x 300km use grid parameter ",grid_calc)
  print(print1)
  print(print2)
  print(print3)
  return(grid_calc)
}
pdc_mercator_proj<-"+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# lcea <- "+proj=cea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" # what I used to use - is the 20, 180 necessary?
lcea <- "+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs" 
years <- c("2008","2009","2010","2011","2012")


# Load data ---------------------------------------------------------------

# JULIA LOADS IN DATA HERE

# Calculate grid value ----------------------------------------------------

# MAKE SP OBJECT so you can assign a CRS - did this before I learned about sf
sp::coordinates(all_data_1) <- c("x", "y")
proj4string(all_data_1) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # placeholder
all_data_1 <- spTransform(all_data_1,lcea)
# Next line is critical for 300km resolution of KDEs
grid_input <- calculate_sp_obj_extent(all_data_1,0.1) 

# Generate individual KDE -------------------------------------------------

# kernelUD output currently is probability density: absolute density (events per unit area) divided by total number of events
results <- kernelUD(all_data_1, grid=grid_input,same4all=T,extent=0.1,h=150000)# grid-input guarantees AT LEAST 300km x 300km - in actuality, it is very slightly more than 300km x 300km.
image(results[[1]])

# Alter lines 87-90 with the results of "names" to see individual tracks. You might not want to do this because you have thousands of tracks
names <- names(results)
lm_kde <- results[grep("lm",names)]
lt_kde <- results[grep("lt",names)]
bm_kde <- results[grep("bm",names)]
bt_kde <- results[grep("bt",names)]
all_LAAL <- append(lm_kde, lt_kde)
all_BFAL <- append(bm_kde, bt_kde)

# Averaging -> new KDEs, account for track length ------------------------

# "0" values are preserved here, not converted to NA
# lt_holder length will change here if your cell size changes (calc_sp_obj_extent)

# average lm #
lm_holder <- numeric(length = 893)
for (i in 1:length(lm_kde)) {
  lm_kde[[i]]@data$ud[is.na(lm_kde[[i]]@data$ud)] <- 0
  add <- lm_kde[[i]]@data$ud
  lm_holder <- lm_holder+add
}
lm_holder <- lm_holder/length(lm_kde)
lm_holder

##### modify existing estUD object with averaged values, then rename
lm_averaged_estUD <- lm_kde[[1]]
lm_averaged_estUD@data$ud <- lm_holder
lm_averaged_estUD@data$ud[is.na(lm_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(lm_kde[[1]])
image(lm_averaged_estUD) 

# average lt #
lt_holder <- numeric(length = 893)
for (i in 1:length(lt_kde)) {
  lt_kde[[i]]@data$ud[is.na(lt_kde[[i]]@data$ud)] <- 0
  add <- lt_kde[[i]]@data$ud
  lt_holder <- lt_holder+add
}
lt_holder <- lt_holder/length(lt_kde)
lt_holder

##### modify existing estUD object with averaged values, then rename
lt_averaged_estUD <- lt_kde[[1]]
lt_averaged_estUD@data$ud <- lt_holder
lt_averaged_estUD@data$ud[is.na(lt_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(lt_kde[[1]])
image(lt_averaged_estUD)         

# average bm #
bm_holder <- numeric(length = 893)
for (i in 1:length(bm_kde)) {
  bm_kde[[i]]@data$ud[is.na(bm_kde[[i]]@data$ud)] <- 0
  add <- bm_kde[[i]]@data$ud
  bm_holder <- bm_holder+add
}
bm_holder <- bm_holder/length(bm_kde)
bm_holder

##### modify existing estUD object with averaged values, then rename
bm_averaged_estUD <- bm_kde[[1]]
bm_averaged_estUD@data$ud <- bm_holder
bm_averaged_estUD@data$ud[is.na(bm_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(bm_kde[[1]])
image(bm_averaged_estUD)

# average bt #
bt_holder <- numeric(length = 893)
for (i in 1:length(bt_kde)) {
  bt_kde[[i]]@data$ud[is.na(bt_kde[[i]]@data$ud)] <- 0
  add <- bt_kde[[i]]@data$ud
  bt_holder <- bt_holder+add
}
bt_holder <- bt_holder/length(bt_kde)
bt_holder

##### modify existing estUD object with averaged values, then rename
bt_averaged_estUD <- bt_kde[[1]]
bt_averaged_estUD@data$ud <- bt_holder
bt_averaged_estUD@data$ud[is.na(bt_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(bt_kde[[1]])
image(bt_averaged_estUD)

# average all_LAAL #
all_LAAL_holder <- numeric(length = 893)
for (i in 1:length(all_LAAL)) {
  all_LAAL[[i]]@data$ud[is.na(all_LAAL[[i]]@data$ud)] <- 0
  add <- all_LAAL[[i]]@data$ud
  all_LAAL_holder <- all_LAAL_holder+add
}
all_LAAL_holder <- all_LAAL_holder/length(all_LAAL)
all_LAAL_holder

##### modify existing estUD object with averaged values, then rename
all_LAAL_averaged_estUD <- all_LAAL[[1]]
all_LAAL_averaged_estUD@data$ud <- all_LAAL_holder
all_LAAL_averaged_estUD@data$ud[is.na(all_LAAL_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(all_LAAL[[1]])
image(all_LAAL_averaged_estUD)

# average all_BFAL #
all_BFAL_holder <- numeric(length = 893)
for (i in 1:length(all_BFAL)) {
  all_BFAL[[i]]@data$ud[is.na(all_BFAL[[i]]@data$ud)] <- 0
  add <- all_BFAL[[i]]@data$ud
  all_BFAL_holder <- all_BFAL_holder+add
}
all_BFAL_holder <- all_BFAL_holder/length(all_BFAL)
all_BFAL_holder

##### modify existing estUD object with averaged values, then rename
all_BFAL_averaged_estUD <- all_BFAL[[1]]
all_BFAL_averaged_estUD@data$ud <- all_BFAL_holder
all_BFAL_averaged_estUD@data$ud[is.na(all_BFAL_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(all_BFAL[[1]])
image(all_BFAL_averaged_estUD)    


# Generate class KDE ------------------------------------------------------

allLAAL <- all_LAAL_averaged_estUD
allBFAL <- all_BFAL_averaged_estUD
midLAAL <- lm_averaged_estUD
midBFAL <- bm_averaged_estUD
ternLAAL <- lt_averaged_estUD
ternBFAL <- bt_averaged_estUD

allLAAL_v_allBFAL <- list(allLAAL,allBFAL)
class(allLAAL_v_allBFAL)<-"estUDm"
names(allLAAL_v_allBFAL)<-c("allLAAL","allBFAL")

# LAAL between islands
midLAAL_v_ternLAAL <- list(midLAAL,ternLAAL)
class(midLAAL_v_ternLAAL)<-"estUDm"
names(midLAAL_v_ternLAAL)<-c("midLAAL","ternLAAL")

# BFAL between islands
midBFAL_v_ternBFAL <- list(midBFAL,ternBFAL)
class(midBFAL_v_ternBFAL)<-"estUDm"
names(midBFAL_v_ternBFAL)<-c("midBFAL","ternBFAL")

# within Midway
midLAAL_v_midBFAL <- list(midLAAL,midBFAL)
class(midLAAL_v_midBFAL)<-"estUDm"
names(midLAAL_v_midBFAL)<-c("midLAAL","midBFAL")

# within Tern
ternLAAL_v_ternBFAL <- list(ternLAAL,ternBFAL)
class(ternLAAL_v_ternBFAL)<-"estUDm"
names(ternLAAL_v_ternBFAL)<-c("ternLAAL","ternBFAL")

# save data classes; now you don't need to re-run all of the above
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/permutation_tests/data_to_load/"
save(allLAAL_v_allBFAL, file=paste0(path,"allLAAL_v_allBFAL.Rdata"))
save(midLAAL_v_ternLAAL, file=paste0(path,"midLAAL_v_ternLAAL.Rdata"))
save(midBFAL_v_ternBFAL, file=paste0(path,"midBFAL_v_ternBFAL.Rdata"))
save(midLAAL_v_midBFAL, file=paste0(path,"midLAAL_v_midBFAL.Rdata"))
save(ternLAAL_v_ternBFAL, file=paste0(path,"ternLAAL_v_ternBFAL.Rdata"))


# Overlap calculations ----------------------------------------------------

### load in files ###

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/permutation_tests/data_to_load/")
path <- getwd()
file_list <- list.files(path=path)

for (i in 1:length(file_list)){
  load(file_list[i])
}

### done loading ###

# all LAAL vs all BFAL
alab_UDOI_95 <- kerneloverlaphr(allLAAL_v_allBFAL, method="UDOI", percent=95, conditional=T)
alab_UDOI_50 <- kerneloverlaphr(allLAAL_v_allBFAL, method="UDOI", percent=50, conditional=T)
alab_PHR_95 <- kerneloverlaphr(allLAAL_v_allBFAL, method="PHR", percent=95, conditional=T)
alab_PHR_50 <- kerneloverlaphr(allLAAL_v_allBFAL, method="PHR", percent=50, conditional=T)
alab_95_BA <- kerneloverlaphr(allLAAL_v_allBFAL, method="BA", percent=95, conditional = T) 
alab_50_BA <- kerneloverlaphr(allLAAL_v_allBFAL, method="BA", percent=50, conditional = T)
alab_UDOI_95_test_stat <- alab_UDOI_95[1,2]
alab_UDOI_50_test_stat <- alab_UDOI_50[1,2]
alab_PHR_95_test_stat <- alab_PHR_95[1,2]
alab_PHR_50_test_stat <- alab_PHR_50[1,2]
alab_95_BA_test_stat <- alab_95_BA[1,2]
alab_50_BA_test_stat <- alab_50_BA[1,2]

# LAAL between islands
mltl_UDOI_95 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="UDOI", percent=95, conditional=T)
mltl_UDOI_50 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="UDOI", percent=50, conditional=T)
mltl_PHR_95 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="PHR", percent=95, conditional=T)
mltl_PHR_50 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="PHR", percent=50, conditional=T)
mltl_95_BA <- kerneloverlaphr(midLAAL_v_ternLAAL, method="BA", percent=95, conditional = T) 
mltl_50_BA <- kerneloverlaphr(midLAAL_v_ternLAAL, method="BA", percent=50, conditional = T)
mltl_UDOI_95_test_stat <- mltl_UDOI_95[1,2]
mltl_UDOI_50_test_stat <- mltl_UDOI_50[1,2]
mltl_PHR_95_test_stat <- mltl_PHR_95[1,2]
mltl_PHR_50_test_stat <- mltl_PHR_50[1,2]
mltl_95_BA_test_stat <- mltl_95_BA[1,2]
mltl_50_BA_test_stat <- mltl_50_BA[1,2]

# BFAL between islands
mbtb_UDOI_95 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="UDOI", percent=95, conditional=T)
mbtb_UDOI_50 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="UDOI", percent=50, conditional=T)
mbtb_PHR_95 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="PHR", percent=95, conditional=T)
mbtb_PHR_50 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="PHR", percent=50, conditional=T)
mbtb_95_BA <- kerneloverlaphr(midBFAL_v_ternBFAL, method="BA", percent=95, conditional = T) 
mbtb_50_BA <- kerneloverlaphr(midBFAL_v_ternBFAL, method="BA", percent=50, conditional = T)
mbtb_UDOI_95_test_stat <- mbtb_UDOI_95[1,2]
mbtb_UDOI_50_test_stat <- mbtb_UDOI_50[1,2]
mbtb_PHR_95_test_stat <- mbtb_PHR_95[1,2]
mbtb_PHR_50_test_stat <- mbtb_PHR_50[1,2]
mbtb_95_BA_test_stat <- mbtb_95_BA[1,2]
mbtb_50_BA_test_stat <- mbtb_50_BA[1,2]

# within Midway
mlmb_UDOI_95 <- kerneloverlaphr(midLAAL_v_midBFAL, method="UDOI", percent=95, conditional=T)
mlmb_UDOI_50 <- kerneloverlaphr(midLAAL_v_midBFAL, method="UDOI", percent=50, conditional=T)
mlmb_PHR_95 <- kerneloverlaphr(midLAAL_v_midBFAL, method="PHR", percent=95, conditional=T)
mlmb_PHR_50 <- kerneloverlaphr(midLAAL_v_midBFAL, method="PHR", percent=50, conditional=T)
mlmb_95_BA <- kerneloverlaphr(midLAAL_v_midBFAL, method="BA", percent=95, conditional = T) 
mlmb_50_BA <- kerneloverlaphr(midLAAL_v_midBFAL, method="BA", percent=50, conditional = T)
mlmb_UDOI_95_test_stat <- mlmb_UDOI_95[1,2]
mlmb_UDOI_50_test_stat <- mlmb_UDOI_50[1,2]
mlmb_PHR_95_test_stat <- mlmb_PHR_95[1,2]
mlmb_PHR_50_test_stat <- mlmb_PHR_50[1,2]
mlmb_95_BA_test_stat <- mlmb_95_BA[1,2]
mlmb_50_BA_test_stat <- mlmb_50_BA[1,2]

# within Tern
tltb_UDOI_95 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="UDOI", percent=95, conditional=T)
tltb_UDOI_50 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="UDOI", percent=50, conditional=T)
tltb_PHR_95 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="PHR", percent=95, conditional=T)
tltb_PHR_50 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="PHR", percent=50, conditional=T)
tltb_95_BA <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="BA", percent=95, conditional = T) 
tltb_50_BA <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="BA", percent=50, conditional = T)
tltb_UDOI_95_test_stat <- tltb_UDOI_95[1,2]
tltb_UDOI_50_test_stat <- tltb_UDOI_50[1,2]
tltb_PHR_95_test_stat <- tltb_PHR_95[1,2]
tltb_PHR_50_test_stat <- tltb_PHR_50[1,2]
tltb_95_BA_test_stat <- tltb_95_BA[1,2]
tltb_50_BA_test_stat <- tltb_50_BA[1,2]

path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/permutation_tests/test_stats/"
setwd(path)
save(list = ls(all.names = TRUE), file = "all_test_stats.RData")