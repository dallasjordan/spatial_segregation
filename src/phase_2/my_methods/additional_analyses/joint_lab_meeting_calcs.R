# Permutation calcs for UDOI for joint lab meeting
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
wrap360 = function(lon) {lon360<-ifelse(lon<0,lon+360,lon);return(lon360)}

# Load data ---------------------------------------------------------------

load_wd <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/julia_simulations/"
# just manual load those in the finder in rStudio

# change this next line for each iteration
all_data_1 <- sim_points_separated
all_data_1 <- all_data_1 %>% select("Animal_ID","lon","lat")
# permutations
all_data <- sim_points_separated
all_data <- all_data %>% select("Animal_ID","lon","lat", "species")


sp::coordinates(all_data_1) <- c("lon", "lat")

# Generate individual KDE -------------------------------------------------

# kernelUD output currently is probability density: absolute density (events per unit area) divided by total number of events
results <- kernelUD(all_data_1,same4all=T)# grid-input guarantees AT LEAST 300km x 300km - in actuality, it is very slightly more than 300km x 300km.
image(results[[1]])

# Alter lines 87-90 with the results of "names" to see individual tracks. You might not want to do this because you have thousands of tracks
names <- names(results)
a_kde <- results[grep("A",names)]
b_kde <- results[grep("B",names)]

# Averaging -> new KDEs, account for track length ------------------------

# "0" values are preserved here, not converted to NA
# lt_holder length will change here if your cell size changes (calc_sp_obj_extent)

# average a #
a_holder <- numeric(length=length(a_kde$A1@data$ud))
for (i in 1:length(a_kde)) {
  a_kde[[i]]@data$ud[is.na(a_kde[[i]]@data$ud)] <- 0
  add <- a_kde[[i]]@data$ud
  a_holder <- a_holder+add
}
a_holder <- a_holder/length(a_kde)
a_holder

##### modify existing estUD object with averaged values, then rename
a_averaged_estUD <- a_kde[[1]]
a_averaged_estUD@data$ud <- a_holder
a_averaged_estUD@data$ud[is.na(a_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(a_kde[[1]])
image(a_averaged_estUD) 



# average b #
b_holder <- numeric(length=length(b_kde$B1@data$ud))
for (i in 1:length(b_kde)) {
  b_kde[[i]]@data$ud[is.na(b_kde[[i]]@data$ud)] <- 0
  add <- b_kde[[i]]@data$ud
  b_holder <- b_holder+add
}
b_holder <- b_holder/length(b_kde)
b_holder

##### modify existing estUD object with averaged values, then rename
b_averaged_estUD <- b_kde[[1]]
b_averaged_estUD@data$ud <- b_holder
b_averaged_estUD@data$ud[is.na(b_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(b_kde[[1]])
image(b_averaged_estUD)        

# Generate class KDE ------------------------------------------------------

a <- a_averaged_estUD
b <- b_averaged_estUD

a_v_b_noOL <- list(a,b)
class(a_v_b_noOL)<-"estUDm"
names(a_v_b_noOL)<-c("a","b")

# save data classes; now you don't need to re-run all of the above
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/julia_simulations/"
save(a_v_b_noOL, file=paste0(path,"a_v_b_noOL.Rdata"))

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
a_v_b_noOL_UDOI_95 <- kerneloverlaphr(a_v_b_noOL, method="UDOI", percent=95, conditional=T)
a_v_b_noOL_UDOI_50 <- kerneloverlaphr(a_v_b_noOL, method="UDOI", percent=50, conditional=T)
a_v_b_noOL_PHR_95 <- kerneloverlaphr(a_v_b_noOL, method="PHR", percent=95, conditional=T)
a_v_b_noOL_PHR_50 <- kerneloverlaphr(a_v_b_noOL, method="PHR", percent=50, conditional=T)
a_v_b_noOL_95_BA <- kerneloverlaphr(a_v_b_noOL, method="BA", percent=95, conditional = T) 
a_v_b_noOL_50_BA <- kerneloverlaphr(a_v_b_noOL, method="BA", percent=50, conditional = T)
a_v_b_noOL_UDOI_95_test_stat <- a_v_b_noOL_UDOI_95[1,2]
a_v_b_noOL_UDOI_50_test_stat <- a_v_b_noOL_UDOI_50[1,2]
a_v_b_noOL_PHR_95_test_stat <- a_v_b_noOL_PHR_95[1,2]
a_v_b_noOL_PHR_50_test_stat <- a_v_b_noOL_PHR_50[1,2]
a_v_b_noOL_95_BA_test_stat <- a_v_b_noOL_95_BA[1,2]
a_v_b_noOL_50_BA_test_stat <- a_v_b_noOL_50_BA[1,2]

# permutations

# Permute lm lt -----------------------------------------------------------

# prep and bookkeeping

# have to run results up above!
results_ud <- results

resample_this <- all_data
resample_this_species <- split(resample_this, resample_this$species)

resample_this_a_track <- split(resample_this_species[[1]],resample_this_species[[1]]$Animal_ID)
resample_this_b_track <- split(resample_this_species[[2]],resample_this_species[[2]]$Animal_ID)
resample_all_tracks <- c(resample_this_a_track,resample_this_b_track)
available_a_tracks <- unique(resample_this_species[[1]]$Animal_ID)
available_b_tracks <- unique(resample_this_species[[2]]$Animal_ID)
n_aat <- length(available_a_tracks) # you are calculating these to preserve sample size in randomizations
n_abt <- length(available_b_tracks)
n_all <- n_aat+n_abt
counts <- resample_this %>% count(species,Animal_ID)

# 95th/50th isopleth and BA
iter <- 1000
significance_tallyUDOI_95 <- 0
significance_tallyUDOI_50 <- 0
significance_tallyBA_95 <- 0
significance_tallyBA_50 <- 0
iter_tally <- 0
resultsUDOI_95_storage <- vector(length=iter, mode="numeric")
resultsUDOI_50_storage <- vector(length=iter, mode="numeric")
resultsBA_95_storage <- vector(length=iter, mode="numeric")
resultsBA_50_storage <- vector(length=iter, mode="numeric")

for (i in 1:iter){
  pool <- 1:n_all
  a_nums <- sample(pool, n_aat, replace=F)
  b_nums <- setdiff(pool,a_nums)
  a_iter <- list()
  b_iter <- list()
  
  for (j in 1:length(a_nums)){
    a_loop <- results_ud[[a_nums[j]]]
    a_iter <- append(a_iter,a_loop)
  }
  
  for (w in 1:length(b_nums)){
    b_loop <- results_ud[[b_nums[w]]]
    b_iter <- append(b_iter,b_loop)
  }
  
  # average a_iter #
  a_iter_holder <- numeric(length = length(a_kde$A1@data$ud))
  for (d in 1:length(a_iter)) {
    a_iter[[d]]@data$ud[is.na(a_iter[[d]]@data$ud)] <- 0
    add <- a_iter[[d]]@data$ud
    a_iter_holder <- a_iter_holder+add
  }
  a_iter_holder <- a_iter_holder/length(a_iter)
  a_iter_holder
  
  ##### modify existing estUD object with averaged values, then rename
  a_iter_avg <- a_iter[[1]]
  a_iter_avg@data$ud <- a_iter_holder
  a_iter_avg@data$ud[is.na(a_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
  image(a_iter[[1]])
  image(a_iter_avg) 
  
  # average b_iter #
  b_iter_holder <- numeric(length = length(a_kde$A1@data$ud))
  for (f in 1:length(b_iter)) {
    b_iter[[f]]@data$ud[is.na(b_iter[[f]]@data$ud)] <- 0
    add <- b_iter[[f]]@data$ud
    b_iter_holder <- b_iter_holder+add
  }
  b_iter_holder <- b_iter_holder/length(b_iter)
  b_iter_holder
  
  ##### modify existing estUD object with averaged values, then rename
  b_iter_avg <- b_iter[[1]]
  b_iter_avg@data$ud <- b_iter_holder
  b_iter_avg@data$ud[is.na(b_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
  image(b_iter[[1]])
  image(b_iter_avg) 
  
  # Create estUDm 
  ia<-a_iter_avg
  ib<-b_iter_avg
  ia_v_ib <- list(ia,ib)
  class(ia_v_ib)<-"estUDm"
  names(ia_v_ib)<-c("ia","ib")
  image(ia_v_ib)
  
  iter_ab_UDOI_95 <- kerneloverlaphr(ia_v_ib, method="UDOI", percent=95, conditional=T)
  iter_ab_UDOI_50 <- kerneloverlaphr(ia_v_ib, method="UDOI", percent=50, conditional=T)
  iter_ab_BA_95 <- kerneloverlaphr(ia_v_ib, method="BA", percent=95, conditional=T) 
  iter_ab_BA_50 <- kerneloverlaphr(ia_v_ib, method="BA", percent=50, conditional=T)
  
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  # proportion of randomized overlaps that are less than observed
  # so if all randomized overlaps are greater than the observed, p=0, there is significant segregation
  # if the observed overlap was less than all randomizations, then P ≤ 0.001, reject the null hypothesis that 
  # there is no difference in the spatial distributions of the two groups
  
  resultsUDOI_95_storage[i] <- iter_ab_UDOI_95[1,2]
  print(resultsUDOI_95_storage[i])
  if (iter_ab_UDOI_95[1,2]<a_v_b_noOL_UDOI_95_test_stat){
    print(iter_ab_UDOI_95)
    print("+1!")
    significance_tallyUDOI_95 <- significance_tallyUDOI_95+1
  }
  
  resultsUDOI_50_storage[i] <- iter_ab_UDOI_50[1,2]
  print(resultsUDOI_50_storage[i])
  if (iter_ab_UDOI_50[1,2]<a_v_b_noOL_UDOI_50_test_stat){
    print(iter_ab_UDOI_50)
    print("+1!")
    significance_tallyUDOI_50 <- significance_tallyUDOI_50+1
  }
  
  resultsBA_95_storage[i] <- iter_ab_BA_95[1,2]
  print(resultsBA_95_storage[i])
  if (iter_ab_BA_95[1,2]<a_v_b_noOL_95_BA_test_stat){
    print(iter_ab_BA_95)
    print("+1!")
    significance_tallyBA_95 <- significance_tallyBA_95+1
  }
  
  resultsBA_50_storage[i] <- iter_ab_BA_50[1,2]
  print(resultsBA_50_storage[i])
  if (iter_ab_BA_50[1,2]<a_v_b_noOL_50_BA_test_stat){
    print(iter_ab_BA_50)
    print("+1!")
    significance_tallyBA_50 <- significance_tallyBA_50+1
  }
}

p_valueUDOI_95 <- significance_tallyUDOI_95/iter
mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)

p_valueUDOI_95 # 0
mean_valueUDOI_95 # 1.286
sd_valueUDOI_95 # 0.055

p_valueUDOI_50 <- significance_tallyUDOI_50/iter
mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)

p_valueUDOI_50 # 0
mean_valueUDOI_50 # 0.149
sd_valueUDOI_50 # 0.025

p_valueBA_95 <- significance_tallyBA_95/iter
mean_valueBA_95 <- mean(resultsBA_95_storage)
sd_valueBA_95 <- sd(resultsBA_95_storage)

p_valueBA_95 # 0
mean_valueBA_95 # 0.887
sd_valueBA_95 # 0.013

p_valueBA_50 <- significance_tallyBA_50/iter
mean_valueBA_50 <- mean(resultsBA_50_storage)
sd_valueBA_50 <- sd(resultsBA_50_storage)

p_valueBA_50 # 0
mean_valueBA_50 # 0.374
sd_valueBA_50 # 0.033

UDOI95 <- as.data.frame(resultsUDOI_95_storage)
ggplot(data=UDOI95, aes(x=resultsUDOI_95_storage))+geom_histogram()+geom_vline(aes(xintercept=a_v_b_noOL_UDOI_95_test_stat, color="red"))

UDOI50 <- as.data.frame(resultsUDOI_50_storage)
ggplot(data=UDOI50, aes(x=resultsUDOI_50_storage))+geom_histogram()+geom_vline(aes(xintercept=a_v_b_noOL_UDOI_50_test_stat, color="red"))
                                                                               