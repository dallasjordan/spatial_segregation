# Test script for permutation methods - 
  # All comparisons follow same pattern, so for testing purposes, just have the comparison of LAAL at Tern and LAAL at Midway

#########################
######## SETUP ##########
#########################

# packages
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

#########################
###### LOAD DATA ########
#########################

# Load LAALmidway and LAALtern postprocessed data

load("INSERT DOWNLOAD FOLDER PATH/LAALdata_midway_withTrackID.Rdata")
LAALmid <- LAAL
LAALmid$id <- paste0("lm",LAALmid$track)
load("/INSERT DOWNLOADS FOLDER PATH/LAALdata_tern_withTrackID.Rdata")
LAALtern <- LAAL
LAALtern$id <- paste0("lt",LAALtern$track)

LAAL <- rbind(LAALmid, LAALtern)

all_data <- LAAL

lm <- all_data[grep("lm", all_data$id), ] # object for just Midway LAAL
lt <- all_data[grep("lt", all_data$id), ] # object for just Tern LAAL

# Load list of calculated KDE, object is called "results"

load("/INSERT DOWNLOADS FOLDER PATH/calculated_KDEs.Rdata")

###########################
### OVERLAP CALCULATIONS ##
###########################

# These are the originally calculated test statistics - there are several steps to get to these values, and they are 
# present in my script "master_script". Here, Just need to assign the values. 

# LAAL between islands
mltl95_test_stat <- 0.9385553
mltl50_test_stat <- 0.09385553 
mltlBA_test_stat <- 0.8509857

######################################################
######################################################

### RANDOMIZATION TESTS ##

######################################################
######################################################


###########################
#### midLAAL v ternLAAL ###
###########################

##### SETUP #####
# Lines 72 - 99 are getting things organized and ready - original dataframes are duplicated and subsetted
# in order to be able to resample entire tracks (permute labels)

    # holder objects (do not want to modify lm and lt)
    add_lm <-lm
    add_lt <-lt
    
    # add columns to make splitting resampling dataframe easier
    add_lm$island <- "LAALmid" 
    add_lt$island <- "LAALtern"
    
    # creating another duplicate dataframe to resample from (the reason this is here is because 'lt and lm' 
    # need to be preserved for other comparisons, not included in this script)
    resample_this <- rbind(add_lm,add_lt)
    
    # split resampling dataframe in islands, then tracks -> necessary in order to permute labels of tracks
    resample_this_island <- split(resample_this, resample_this$island)
    resample_this_mid_track <- split(resample_this_island[[1]],resample_this_island[[1]]$track)
    resample_this_tern_track <- split(resample_this_island[[2]],resample_this_island[[2]]$track)
    
    # combined list of individual tracks
    resample_all_tracks <- c(resample_this_mid_track,resample_this_tern_track)
    
    # tallying counts of each type of track, used to make sure all tracks were preserved between
    # original dataframe and the creation of resampling dataframes. Variables used in initial resampling 
    # for each loop iteration.
    available_mid_tracks <- unique(resample_this_island[[1]]$track)
    available_tern_tracks <- unique(resample_this_island[[2]]$track)
    n_amt <- length(available_mid_tracks) # you are calculating these to preserve sample size in randomizations
    n_att <- length(available_tern_tracks)
    n_all <- n_amt+n_att # n_all is number of available tracks to randomize with
    counts <- resample_this %>% count(id,track)

###### BEGIN PERMUTATION TESTS ######
    
# 95th/50th isopleth and BA
iter <- 1000 # one thousand iterations takes a long time; you can test with a smaller number here. 
significance_tally95 <- 0
significance_tally50 <- 0
significance_tallyBA <- 0
iter_tally <- 0
results95_storage <- vector()
results50_storage <- vector()
resultsBA_storage <- vector()
for (i in 1:iter){
  pool <- 1:n_all # reminder, n_all is the number of total tracks that you are randomizing
  mid_nums <- sample(pool, n_amt, replace=F) # designate which tracks are going to be assigned "midway" from the available tracks
  tern_nums <- setdiff(pool,mid_nums) # designate which tracks are going to be assigned "tern" from the available tracks
  mid_iter <- NULL
  tern_iter <- NULL
  for (j in 1:length(mid_nums)){
    mid_loop <- resample_all_tracks[[mid_nums[j]]] # here, this step pulls the tracks assigned to Midway from the list of all tracks
    mid_iter <- rbind(mid_iter,mid_loop) # The pulled tracks are now mid_iter (some of these tracks were originally Tern tracks)
  }
  mid_tracks <- unique(mid_iter$id) # can call "mid_tracks" to see that what is Midway for this randomization includes both originally tern and midway tracks
  mid_iter_ud <- results[mid_tracks] # IMPORTANT! grabbing the utilization distributions (KDEs) from the KDE object. We are doing this because we have to normalize the tracks within this procedure.
  for (j in 1:length(tern_nums)){
    tern_loop <- resample_all_tracks[[tern_nums[j]]] # doing what we did above, for tern
    tern_iter <- rbind(tern_iter,tern_loop)
  }
  tern_tracks <- unique(tern_iter$id) # doing what we did above, for tern 
  tern_iter_ud <- results[tern_tracks]
  
  ##### normalize Midway group ##### 
  m_iter_holder <- numeric(length = 570) # All KDEs were calculated on the same grid (at least 300km), which left 570 cells in the underlying raster. This value is hard coded here.
  for (k in 1:length(mid_iter_ud)) { 
    add <- mid_iter_ud[[k]]@data$ud # Description of lines 125 - 129: this portion is the normalization of the density value. It is adding respective cell values for every KDE then dividing by the number of KDEs to have one normalized KDE
    m_iter_holder <- m_iter_holder+add # Line 125-129 need examining! High potential something is incorrect. However,normalization will not drastically change the overall problem we are trying to address (still worth examining though)
  }
  m_iter_holder <- m_iter_holder/length(mid_iter_ud)
  m_iter_holder
  
      ##### modify existing estUD object with averaged values, then rename
      m_iter_normalized_estUD <- mid_iter_ud[[1]] # Replicating one KDE to serve as a template raster
      m_iter_normalized_estUD@data$ud <- m_iter_holder # Overwriting template with the normalized UD
      image(mid_iter_ud[[1]]) # prints so you can see it as it loops
      image(m_iter_normalized_estUD) # prints so you can see it as it loops
  
  ##### normalize Tern group ##### 
  t_iter_holder <- numeric(length = 570) # Lines 138 - 150 are the same as 123-135, just for the tracks assigned randomly as "Tern"
  for (m in 1:length(tern_iter_ud)) {
    add <- tern_iter_ud[[m]]@data$ud
    t_iter_holder <- t_iter_holder+add
  }
  t_iter_holder <- t_iter_holder/length(tern_iter_ud)
  t_iter_holder
  
      ##### modify existing estUD object with averaged values, then rename
      t_iter_normalized_estUD <- tern_iter_ud[[1]]
      t_iter_normalized_estUD@data$ud <- t_iter_holder
      image(tern_iter_ud[[1]])
      image(t_iter_normalized_estUD) 
  
  ##### Create estUDm ##### 
  imLAAL<-m_iter_normalized_estUD # Lines 153-157 are just adjusting classes to get them to work with esoteric adehabitat functions
  itLAAL<-t_iter_normalized_estUD
  imLAAL_v_itLAAL <- list(imLAAL,itLAAL)
  class(imLAAL_v_itLAAL)<-"estUDm"
  names(imLAAL_v_itLAAL)<-c("imLAAL","itLAAL")
  
  iter_mltl_95 <- kerneloverlaphr(imLAAL_v_itLAAL, method="UDOI", percent=95) # The actual overlap calculation for 95th percentile 
  iter_mltl_50 <- kerneloverlaphr(imLAAL_v_itLAAL, method="UDOI", percent=50) # The actual overlap calculation for 50th percentile
  iter_mltl_BA <- kerneloverlaphr(imLAAL_v_itLAAL, method="BA") # The actual overlap calculation for BA
  
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  # Storing the results of the overlap calculation - fairly self explanatory, it "dings" whenever there is a randomized overlap 
  # that is less than the originally calculated value. Out of a 1000 iterations, the number of times that the randomized value
  # is less than the original values is the signficance. IMPORTANT: Hence, if randomly generated overlap values being lower than the original 
  # data overlap is rare, we have a low p-value -> significant segregation. 
  results95_storage[i] <- iter_mltl_95[1,2]
  print(results95_storage[i])
  if (iter_mltl_95[1,2]<mltl95_test_stat){
    print(iter_mltl_95)
    print("+1!")
    significance_tally95 <- significance_tally95+1
  }
  
  results50_storage[i] <- iter_mltl_50[1,2]
  print(results50_storage[i])
  if (iter_mltl_50[1,2]<mltl50_test_stat){
    print(iter_mltl_50)
    print("+1!")
    significance_tally50 <- significance_tally50+1
  }
  
  resultsBA_storage[i] <- iter_mltl_BA[1,2]
  print(resultsBA_storage[i])
  if (iter_mltl_BA[1,2]<mltlBA_test_stat){
    print(iter_mltl_BA)
    print("+1!")
    significance_tallyBA <- significance_tallyBA+1
  }
}

# display of final values - the commented numbers are what my last randomization tests calculated. 

p_value95 <- significance_tally95/iter
mean_value95 <- mean(results95_storage)
sd_value95 <- sd(results95_storage)

p_value95 # 0
mean_value95 # 1.446012
sd_value95 # 0.05678921

p_value50 <- significance_tally50/iter
mean_value50 <- mean(results50_storage)
sd_value50 <- sd(results50_storage)

p_value50 # 0
mean_value50 # 0.2648847
sd_value50 # 0.02871768

p_valueBA <- significance_tallyBA/iter
mean_valueBA <- mean(resultsBA_storage)
sd_valueBA <- sd(resultsBA_storage)

p_valueBA # 0
mean_valueBA # 0.9740942
sd_valueBA # 0.007144034











