# Permutation test with permutation of track labels - midwayLAAL and ternLAAL
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
# for mac
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

# for pc
load("G:/Other computers/My MacBook Pro/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

# grid parameter calc
all_data_1<-all_data[,c(1,3,4)]
sp::coordinates(all_data_1) <- c("x", "y")
proj4string(all_data_1) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # placeholder
all_data_1 <- spTransform(all_data_1,lcea)
grid_input <- calculate_sp_obj_extent(all_data_1,0.1)

# LAAL kernels
lk <- rbind(lm,lt)
lk <- lk[,c(1,3,4)]
sp::coordinates(lk) <- c("x", "y")
proj4string(lk) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # placeholder
lk <- spTransform(lk,lcea)
grid_input <- calculate_sp_obj_extent(lk,0.1)

a_ud <-  kernelUD(lk, grid=grid_input,same4all=T,extent=0.1,h=150000)
image(a_ud)

# for mac
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/permutation_tests/test_stats/"
setwd(path)
load("all_test_stats.RData") # generated from master_script

# for pc
path <- "G:/Other computers/My MacBook Pro/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/permutation_tests/test_stats/"
setwd(path)
load("all_test_stats.RData") 

# Permute lm lt -----------------------------------------------------------

# prep and bookkeeping
add_lm <-lm
add_lt <-lt
add_lm$island <- "LAALmid"
add_lt$island <- "LAALtern"
resample_this <- rbind(add_lm,add_lt)
resample_this_island <- split(resample_this, resample_this$island)
resample_this_mid_track <- split(resample_this_island[[1]],resample_this_island[[1]]$track)
resample_this_tern_track <- split(resample_this_island[[2]],resample_this_island[[2]]$track)
resample_all_tracks <- c(resample_this_mid_track,resample_this_tern_track)
available_mid_tracks <- unique(resample_this_island[[1]]$track)
available_tern_tracks <- unique(resample_this_island[[2]]$track)
n_amt <- length(available_mid_tracks) # you are calculating these to preserve sample size in randomizations
n_att <- length(available_tern_tracks)
n_all <- n_amt+n_att
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
  mid_nums <- sample(pool, n_amt, replace=F)
  tern_nums <- setdiff(pool,mid_nums)
  mid_iter <- list()
  tern_iter <- list()
  
  for (j in 1:length(mid_nums)){
    mid_loop <- a_ud[[mid_nums[j]]]
    mid_iter <- append(mid_iter,mid_loop)
  }

  for (w in 1:length(tern_nums)){
    tern_loop <- a_ud[[tern_nums[w]]]
    tern_iter <- append(tern_iter,tern_loop)
  }
    
      # average mid_iter #
      mid_iter_holder <- numeric(length = 640)
      for (d in 1:length(mid_iter)) {
        mid_iter[[d]]@data$ud[is.na(mid_iter[[d]]@data$ud)] <- 0
        add <- mid_iter[[d]]@data$ud
        mid_iter_holder <- mid_iter_holder+add
      }
      mid_iter_holder <- mid_iter_holder/length(mid_iter)
      mid_iter_holder
      
      ##### modify existing estUD object with averaged values, then rename
      mid_iter_avg <- mid_iter[[1]]
      mid_iter_avg@data$ud <- mid_iter_holder
      mid_iter_avg@data$ud[is.na(mid_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
      image(mid_iter[[1]])
      image(mid_iter_avg) 
      
      # average tern_iter #
      tern_iter_holder <- numeric(length = 640)
      for (f in 1:length(tern_iter)) {
        tern_iter[[f]]@data$ud[is.na(tern_iter[[f]]@data$ud)] <- 0
        add <- tern_iter[[f]]@data$ud
        tern_iter_holder <- tern_iter_holder+add
      }
      tern_iter_holder <- tern_iter_holder/length(tern_iter)
      tern_iter_holder
      
      ##### modify existing estUD object with averaged values, then rename
      tern_iter_avg <- tern_iter[[1]]
      tern_iter_avg@data$ud <- tern_iter_holder
      tern_iter_avg@data$ud[is.na(tern_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
      image(tern_iter[[1]])
      image(tern_iter_avg) 
  
  # Create estUDm 
  imLAAL<-mid_iter_avg
  itLAAL<-tern_iter_avg
  imLAAL_v_itLAAL <- list(imLAAL,itLAAL)
  class(imLAAL_v_itLAAL)<-"estUDm"
  names(imLAAL_v_itLAAL)<-c("imLAAL","itLAAL")
  image(imLAAL_v_itLAAL)
  
  iter_mltl_UDOI_95 <- kerneloverlaphr(imLAAL_v_itLAAL, method="UDOI", percent=95, conditional=T)
  iter_mltl_UDOI_50 <- kerneloverlaphr(imLAAL_v_itLAAL, method="UDOI", percent=50, conditional=T)
  iter_mltl_BA_95 <- kerneloverlaphr(imLAAL_v_itLAAL, method="BA", percent=95, conditional=T) 
  iter_mltl_BA_50 <- kerneloverlaphr(imLAAL_v_itLAAL, method="BA", percent=50, conditional=T)
  
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  # proportion of randomized overlaps that are less than observed
  # so if all randomized overlaps are greater than the observed, p=0, there is significant segregation
  # if the observed overlap was less than all randomizations, then P ≤ 0.001, reject the null hypothesis that 
  # there is no difference in the spatial distributions of the two groups
  
  resultsUDOI_95_storage[i] <- iter_mltl_UDOI_95[1,2]
  print(resultsUDOI_95_storage[i])
  if (iter_mltl_UDOI_95[1,2]<mltl_UDOI_95_test_stat){
    print(iter_mltl_UDOI_95)
    print("+1!")
    significance_tallyUDOI_95 <- significance_tallyUDOI_95+1
  }
  
  resultsUDOI_50_storage[i] <- iter_mltl_UDOI_50[1,2]
  print(resultsUDOI_50_storage[i])
  if (iter_mltl_UDOI_50[1,2]<mltl_UDOI_50_test_stat){
    print(iter_mltl_UDOI_50)
    print("+1!")
    significance_tallyUDOI_50 <- significance_tallyUDOI_50+1
  }
  
  resultsBA_95_storage[i] <- iter_mltl_BA_95[1,2]
  print(resultsBA_95_storage[i])
  if (iter_mltl_BA_95[1,2]<mltl_95_BA_test_stat){
    print(iter_mltl_BA_95)
    print("+1!")
    significance_tallyBA_95 <- significance_tallyBA_95+1
  }
  
  resultsBA_50_storage[i] <- iter_mltl_BA_50[1,2]
  print(resultsBA_50_storage[i])
  if (iter_mltl_BA_50[1,2]<mltl_50_BA_test_stat){
    print(iter_mltl_BA_50)
    print("+1!")
    significance_tallyBA_50 <- significance_tallyBA_50+1
  }
}

# values updated as of Nov 23 2021
p_valueUDOI_95 <- significance_tallyUDOI_95/iter
mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)

p_valueUDOI_95 # 0
mean_valueUDOI_95 # 1.290
sd_valueUDOI_95 # 0.054

p_valueUDOI_50 <- significance_tallyUDOI_50/iter
mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)

p_valueUDOI_50 # 0.001
mean_valueUDOI_50 # 0.148
sd_valueUDOI_50 # 0.025

p_valueBA_95 <- significance_tallyBA_95/iter
mean_valueBA_95 <- mean(resultsBA_95_storage)
sd_valueBA_95 <- sd(resultsBA_95_storage)

p_valueBA_95 # 0
mean_valueBA_95 # 0.887
sd_valueBA_95 # 0.012

p_valueBA_50 <- significance_tallyBA_50/iter
mean_valueBA_50 <- mean(resultsBA_50_storage)
sd_valueBA_50 <- sd(resultsBA_50_storage)

p_valueBA_50 # 0.001
mean_valueBA_50 # 0.373
sd_valueBA_50 # 0.032









