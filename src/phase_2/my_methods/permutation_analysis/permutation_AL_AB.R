# Permutation test with permutation of track labels - allLAAL and allBFAL
# This is an updated version of the first attempt you had at randomization tests - in this 
# version, we are following Clay et al. 2016 and permuting track labels (preserving sample size
# of original groups) and then creating kernels -> calculating overlap, comparing to observed

# What you will need
# sample size of groups
# the tracks 
# test stats

# Pseudocode
# save sample size of groups
# load in tracks
# load in test stats
# take the tracks -> randomly assign into a group for a specific comparison
# preserve sample size!
# create kernelUD
# calculate overlap
# compare to test stat


# Setup -------------------------------------------------------------------

library(adehabitatHR)
library(dplyr)
library(sp)

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

#load tracks in
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

al <- rbind(lm,lt)
ab <- rbind(bm,bt)



# grid parameter calc
all_data_1<-all_data[,c(1,3,4)]
sp::coordinates(all_data_1) <- c("x", "y")
proj4string(all_data_1) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # placeholder
all_data_1 <- spTransform(all_data_1,lcea)
grid_input <- calculate_sp_obj_extent(all_data_1,0.1)

# all kernels
ak <- rbind(al,ab)
ak <- ak[,c(1,3,4)]
sp::coordinates(ak) <- c("x", "y")
proj4string(ak) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # placeholder
ak <- spTransform(ak,lcea)
grid_input <- calculate_sp_obj_extent(ak,0.1)

a_ud <-  kernelUD(ak, grid=grid_input,same4all=T,extent=0.1,h=150000)
image(a_ud)

path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/permutation_tests/test_stats/"
setwd(path)
load("all_test_stats.RData") # generated from master_script

# Permute lm bm -----------------------------------------------------------

# prep and bookkeeping
add_al <-al
add_ab <-ab
add_al$species <- "LAAL"
add_ab$species <- "BFAL"
resample_this <- rbind(add_al,add_ab)
resample_this_species <- split(resample_this, resample_this$species)
resample_this_BFAL_track <- split(resample_this_species[[1]],resample_this_species[[1]]$id)
resample_this_LAAL_track <- split(resample_this_species[[2]],resample_this_species[[2]]$id)
resample_all_tracks <- c(resample_this_BFAL_track,resample_this_LAAL_track)
available_BFAL_tracks <- unique(resample_this_species[[1]]$id)
available_LAAL_tracks <- unique(resample_this_species[[2]]$id)
n_abt <- length(available_BFAL_tracks) # you are calculating these to preserve sample size in randomizations
n_alt <- length(available_LAAL_tracks)
n_all <- n_abt+n_alt
counts <- resample_this %>% count(id,track)

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
  BFAL_nums <- sample(pool, n_abt, replace=F)
  LAAL_nums <- setdiff(pool,BFAL_nums)
  BFAL_iter <- list()
  LAAL_iter <- list()
  
  for (j in 1:length(BFAL_nums)){
    BFAL_loop <- a_ud[[BFAL_nums[j]]]
    BFAL_iter <- append(BFAL_iter,BFAL_loop)
  }
  
  for (w in 1:length(LAAL_nums)){
    LAAL_loop <- a_ud[[LAAL_nums[w]]]
    LAAL_iter <- append(LAAL_iter,LAAL_loop)
  }
  
  # average mid_iter #
  BFAL_iter_holder <- numeric(length = 893)
  for (d in 1:length(BFAL_iter)) {
    BFAL_iter[[d]]@data$ud[is.na(BFAL_iter[[d]]@data$ud)] <- 0
    add <- BFAL_iter[[d]]@data$ud
    BFAL_iter_holder <- BFAL_iter_holder+add
  }
  BFAL_iter_holder <- BFAL_iter_holder/length(BFAL_iter)
  BFAL_iter_holder
  
  ##### modify existing estUD object with averaged values, then rename
  BFAL_iter_avg <- BFAL_iter[[1]]
  BFAL_iter_avg@data$ud <- BFAL_iter_holder
  BFAL_iter_avg@data$ud[is.na(BFAL_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
  image(BFAL_iter[[1]])
  image(BFAL_iter_avg) 
  
  # average tern_iter #
  LAAL_iter_holder <- numeric(length = 893)
  for (f in 1:length(LAAL_iter)) {
    LAAL_iter[[f]]@data$ud[is.na(LAAL_iter[[f]]@data$ud)] <- 0
    add <- LAAL_iter[[f]]@data$ud
    LAAL_iter_holder <- LAAL_iter_holder+add
  }
  LAAL_iter_holder <- LAAL_iter_holder/length(LAAL_iter)
  LAAL_iter_holder
  
  ##### modify existing estUD object with averaged values, then rename
  LAAL_iter_avg <- LAAL_iter[[1]]
  LAAL_iter_avg@data$ud <- LAAL_iter_holder
  LAAL_iter_avg@data$ud[is.na(LAAL_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
  image(LAAL_iter[[1]])
  image(LAAL_iter_avg) 
  
  # Create estUDm 
  iBFAL<-BFAL_iter_avg
  iLAAL<-LAAL_iter_avg
  iLAAL_v_iBFAL <- list(iLAAL,iBFAL)
  class(iLAAL_v_iBFAL)<-"estUDm"
  names(iLAAL_v_iBFAL)<-c("iLAAL","iBFAL")
  image(iLAAL_v_iBFAL)
  
  iter_alab_UDOI_95 <- kerneloverlaphr(iLAAL_v_iBFAL, method="UDOI", percent=95, conditional=T)
  iter_alab_UDOI_50 <- kerneloverlaphr(iLAAL_v_iBFAL, method="UDOI", percent=50, conditional=T)
  iter_alab_BA_95 <- kerneloverlaphr(iLAAL_v_iBFAL, method="BA", percent=95, conditional=T) 
  iter_alab_BA_50 <- kerneloverlaphr(iLAAL_v_iBFAL, method="BA", percent=50, conditional=T)
  
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  # proportion of randomized overlaps that are less than observed
  # so if all randomized overlaps are greater than the observed, p=0, there is significant segregation
  # if the observed overlap was less than all randomizations, then P ≤ 0.001, reject the null hypothesis that 
  # there is no difference in the spatial distributions of the two groups
  
  resultsUDOI_95_storage[i] <- iter_alab_UDOI_95[1,2]
  print(resultsUDOI_95_storage[i])
  if (iter_alab_UDOI_95[1,2]<alab_UDOI_95_test_stat){
    print(iter_alab_UDOI_95)
    print("+1!")
    significance_tallyUDOI_95 <- significance_tallyUDOI_95+1
  }
  
  resultsUDOI_50_storage[i] <- iter_alab_UDOI_50[1,2]
  print(resultsUDOI_50_storage[i])
  if (iter_alab_UDOI_50[1,2]<alab_UDOI_50_test_stat){
    print(iter_alab_UDOI_50)
    print("+1!")
    significance_tallyUDOI_50 <- significance_tallyUDOI_50+1
  }
  
  resultsBA_95_storage[i] <- iter_alab_BA_95[1,2]
  print(resultsBA_95_storage[i])
  if (iter_alab_BA_95[1,2]<alab_95_BA_test_stat){
    print(iter_alab_BA_95)
    print("+1!")
    significance_tallyBA_95 <- significance_tallyBA_95+1
  }
  
  resultsBA_50_storage[i] <- iter_alab_BA_50[1,2]
  print(resultsBA_50_storage[i])
  if (iter_alab_BA_50[1,2]<alab_50_BA_test_stat){
    print(iter_alab_BA_50)
    print("+1!")
    significance_tallyBA_50 <- significance_tallyBA_50+1
  }
}

# values up to date as of dec 19 2021

p_valueUDOI_95 <- significance_tallyUDOI_95/iter
mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)

p_valueUDOI_95 # <0.001
mean_valueUDOI_95 # 1.191
sd_valueUDOI_95 # 0.400

p_valueUDOI_50 <- significance_tallyUDOI_50/iter
mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)

p_valueUDOI_50 # <0.001
mean_valueUDOI_50 # 0.157
sd_valueUDOI_50 # 0.200

p_valueBA_95 <- significance_tallyBA_95/iter
mean_valueBA_95 <- mean(resultsBA_95_storage)
sd_valueBA_95 <- sd(resultsBA_95_storage)

p_valueBA_95 # <0.001
mean_valueBA_95 # 0.887
sd_valueBA_95 # 0.012

p_valueBA_50 <- significance_tallyBA_50/iter
mean_valueBA_50 <- mean(resultsBA_50_storage)
sd_valueBA_50 <- sd(resultsBA_50_storage)

p_valueBA_50 # 0
mean_valueBA_50 # 0.386
sd_valueBA_50 # 0.024









