# *New* Master KDE Script
# This script: 
  # 1. Loads post-processed data
  # 2. Creates KDEs using the adehabitatHR package for each individual, saves as "results"
  # 3. Normalizes these for the relevant group (there are 5 total, LAAL/BFAL for each island and combined)
  # 4. Assigns these estUD objects to coded names for overlap comparisons, e.g. lm_normalized_estUD = LAALmidway
  # 5. Makes estUDm objects for each comparison (combines estUD objects pairwise according to the comparison of interest)
  # 6. Calculates overlap of 95 UD, 50 UD, and BA for each comparison of original data
  # 7. Randomization procedure


# Dallas Jordan April 1 2021

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

# functions and objects
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
years <- c("2008","2009","2010","2011","2012")

#########################
###### LOAD DATA ########
#########################

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

# comparisons: allLAAL v all BFAL
#              ternLAAL v midwayLAAL
#              ternBFAL v midwayBFAL
#              ternLAAL v ternBFAL
#              midwayLAAL v midwayBFAL
all_data <- rbind(LAAL,BFAL)
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

#########################
#### CALC GRID VALUE ####
#########################

all_data_1<-all_data[,1:3]
sp::coordinates(all_data_1) <- c("x", "y")
proj4string(all_data_1) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs") # placeholder
all_data_1 <- spTransform(all_data_1,laea)
grid_input <- calculate_sp_obj_extent(all_data_1)
# https://trac.faunalia.it/animove/wiki/AnimoveHowto#Whatunitisthegridcellsizein

###########################
# GENERATE INDIVIDUAL KDE #
###########################

# kernelUD output currently is probability density: absolute density divided by total number of events
results <- kernelUD(all_data_1, grid=grid_input,same4all=T,extent=0.1) # grid-input guarantees AT LEAST 300km x 300km - in actuality, it is more than 300km x 300km.
  # image(results)

###############################################
# NORMALIZATION THROUGH AVERAGING -> NEW KDEs #
###############################################

# average all individuals in an estUDm, save as estUD
names <- names(results)
lm_kde <- results[grep("lm",names)]
lt_kde <- results[grep("lt",names)]
bm_kde <- results[grep("bm",names)]
bt_kde <- results[grep("bt",names)]

##### normalize lm ##### 
    lm_holder <- numeric(length = 570)
    for (i in 1:length(lm_kde)) {
      add <- lm_kde[[i]]@data$ud
      lm_holder <- lm_holder+add
    }
    lm_holder <- lm_holder/length(lm_kde)
    lm_holder
    
    ##### modify existing estUD object with averaged values, then rename
        lm_normalized_estUD <- lm_kde[[1]]
        lm_normalized_estUD@data$ud <- lm_holder
        image(lm_kde[[1]])
        image(lm_normalized_estUD) 

##### normalize lt ##### 
    lt_holder <- numeric(length = 570)
    for (i in 1:length(lt_kde)) {
      add <- lt_kde[[i]]@data$ud
      lt_holder <- lt_holder+add
    }
    lt_holder <- lt_holder/length(lt_kde)
    lt_holder
        
    ##### modify existing estUD object with averaged values, then rename
        lt_normalized_estUD <- lt_kde[[1]]
        lt_normalized_estUD@data$ud <- lt_holder
        image(lt_kde[[1]])
        image(lt_normalized_estUD)         
    
##### normalize bm ##### 
    bm_holder <- numeric(length = 570)
    for (i in 1:length(bm_kde)) {
      add <- bm_kde[[i]]@data$ud
      bm_holder <- bm_holder+add
    }
    bm_holder <- bm_holder/length(bm_kde)
    bm_holder
        
    ##### modify existing estUD object with averaged values, then rename
        bm_normalized_estUD <- bm_kde[[1]]
        bm_normalized_estUD@data$ud <- bm_holder
        image(bm_kde[[1]])
        image(bm_normalized_estUD)
    
##### normalize bt ##### 
    bt_holder <- numeric(length = 570)
    for (i in 1:length(bt_kde)) {
      add <- bt_kde[[i]]@data$ud
      bt_holder <- bt_holder+add
    }
    bt_holder <- bt_holder/length(bt_kde)
    bt_holder
        
    ##### modify existing estUD object with averaged values, then rename
        bt_normalized_estUD <- bt_kde[[1]]
        bt_normalized_estUD@data$ud <- bt_holder
        image(bt_kde[[1]])
        image(bt_normalized_estUD)

###########################
### GENERATE CLASS KDE ####
###########################
# udspdf <- estUDm2spixdf(lm_kde)

midLAAL <- lm_normalized_estUD
midBFAL <- bm_normalized_estUD
ternLAAL <- lt_normalized_estUD
ternBFAL <- bt_normalized_estUD

all_classes <- list(midLAAL,midBFAL,ternLAAL,ternBFAL)
class(all_classes)<-"estUDm"
names(all_classes)<-c("midLAAL","midBFAL","ternLAAL","ternBFAL")

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

###########################
### OVERLAP CALCULATIONS ##
###########################

# LAAL between islands
mltl_95 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="UDOI", percent=95, conditional = T)
mltl_50 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="UDOI", percent=50)
mltl_BA <- kerneloverlaphr(midLAAL_v_ternLAAL, method="BA") # 0.8509857
mltl95_test_stat <- 0.9385553 # manually putting this in because its easier then calling from matrix
mltl50_test_stat <- 0.09385553 # manually putting this in because its easier then calling from matrix
mltlBA_test_stat <- 0.8509857

# BFAL between islands
mbtb_95 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="UDOI", percent=95)
mbtb_50 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="UDOI", percent=50)
mbtb_BA <- kerneloverlaphr(midBFAL_v_ternBFAL, method="BA") # 0.7314891
mbtb95_test_stat <- 0.5460556 # manually putting this in because its easier then calling from matrix
mbtb50_test_stat <- 0.02642205 # manually putting this in because its easier then calling from matrix
mbtbBA_test_stat <- 0.7314891

# within Midway
mlmb_95 <- kerneloverlaphr(midLAAL_v_midBFAL, method="UDOI", percent=95)
mlmb_50 <- kerneloverlaphr(midLAAL_v_midBFAL, method="UDOI", percent=50)
mlmb_BA <- kerneloverlaphr(midLAAL_v_midBFAL, method="BA") # 0.8230154 
mlmb95_test_stat <- 0.7865138 # manually putting this in because its easier then calling from matrix
mlmb50_test_stat <-	0.0453758	# manually putting this in because its easier then calling from matrix
mlmbBA_test_stat <- 0.8230154 

# within Tern
tltb_95 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="UDOI", percent=95)
tltb_50 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="UDOI", percent=50)
tltb_BA <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="BA") # 0.5907059
tltb95_test_stat <- 0.2868636  # manually putting this in because its easier then calling from matrix
tltb50_test_stat <-	0.002957356# manually putting this in because its easier then calling from matrix
tltbBA_test_stat <- 0.5907059


######################################################
######################################################

### RANDOMIZATION TESTS ##

######################################################
######################################################

###########################
#### midLAAL v ternLAAL ###
###########################

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
  iter <- 100
  significance_tally95 <- 0
  significance_tally50 <- 0
  significance_tallyBA <- 0
  iter_tally <- 0
  results95_storage <- vector()
  results50_storage <- vector()
  resultsBA_storage <- vector()
  for (i in 1:iter){
    pool <- 1:n_all
    mid_nums <- sample(pool, n_amt, replace=F)
    tern_nums <- setdiff(pool,mid_nums)
    mid_iter <- NULL
    tern_iter <- NULL
    for (j in 1:length(mid_nums)){
      mid_loop <- resample_all_tracks[[mid_nums[j]]]
      mid_iter <- rbind(mid_iter,mid_loop)
    }
    mid_tracks <- unique(mid_iter$id)
    mid_iter_ud <- results[mid_tracks]
    for (j in 1:length(tern_nums)){
      tern_loop <- resample_all_tracks[[tern_nums[j]]]
      tern_iter <- rbind(tern_iter,tern_loop)
    }
    tern_tracks <- unique(tern_iter$id)
    tern_iter_ud <- results[tern_tracks]
    ##### normalize ##### 
            m_iter_holder <- numeric(length = 570) # if grid parameter gets changed, gotta change this
            for (k in 1:length(mid_iter_ud)) {
              add <- mid_iter_ud[[k]]@data$ud
              m_iter_holder <- m_iter_holder+add
            }
            m_iter_holder <- m_iter_holder/length(mid_iter_ud)
            m_iter_holder
            
                ##### modify existing estUD object with averaged values, then rename
                m_iter_normalized_estUD <- mid_iter_ud[[1]]
                m_iter_normalized_estUD@data$ud <- m_iter_holder
                image(mid_iter_ud[[1]])
                image(m_iter_normalized_estUD) 
            
            t_iter_holder <- numeric(length = 570) # if grid parameter gets changed, gotta change this
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
    imLAAL<-m_iter_normalized_estUD
    itLAAL<-t_iter_normalized_estUD
    imLAAL_v_itLAAL <- list(imLAAL,itLAAL)
    class(imLAAL_v_itLAAL)<-"estUDm"
    names(imLAAL_v_itLAAL)<-c("imLAAL","itLAAL")
    
    iter_mltl_95 <- kerneloverlaphr(imLAAL_v_itLAAL, method="UDOI", percent=95)
    iter_mltl_50 <- kerneloverlaphr(imLAAL_v_itLAAL, method="UDOI", percent=50)
    iter_mltl_BA <- kerneloverlaphr(imLAAL_v_itLAAL, method="BA") 
    
    iter_tally <- iter_tally+1
    print(iter_tally)
    
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
  
  
###########################
#### midBFAL v ternBFAL ###
###########################