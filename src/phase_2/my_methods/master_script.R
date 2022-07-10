# Master KDE Script
# Dallas Jordan 
# Originated: Dec 18 2021
# Last edited: June 7 2021

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

####### IMPORTANT NOTES, in no particular order. Especially important are notes concerning bandwidth, 'h'. 

# no need for asymptote analysis because of Gutowsky paper that Melinda found

# Normalization not needed when dealing with probability density UDs

# Brownian bridge UDs are newer, but they make assumptions about staying at foraging spots that 
# aren't necessarily true here

# smoothing parameter 'h' is in the units of your relocations - here, I have projected into meters, so smoothing parameter should
# be in meters. The UD is estimated at the center of each pixel of a grid - adehabitatHR vignette. I set 150km 'h', because I am assuming that each kernel is smoothed from the center
# of a grid, so 150km in each direction = 300km. From the vignette - the size and resolution of the grid does not have a large effect on the estimates

# Concerning the value of h/bandwidth/search window - Low values of h give nearby locations the greatest influence on the shape of the kernel, revealing small-scale detail, while large 
# values allow more influence from distant locations, which reveals the outlying shape of the distribution

# LSCV only works when you have frequent and NONAUTOCORRELATED (Row, 2006) location data and there aren't many 'islands' of relocation. It
# is calculated by minimizing the error by comparing the prediction from all data points to the data minus each point.
# More about better smoothing parameters: https://vita.had.co.nz/papers/density-estimation.pdf

# href tends to estimate a larger home range than other methods (oversmoothing). However, we aren't interested
# in the size of the home range here; rather we are interested in the relative use. So to meet our needs, we just need 
# to pick a consistent bandwidth that will highlight different densities of use across space. We pick a consistent bandwidth
# because it will assume that space use is uniform (no area of the ocean gets a bigger search window to boost its density)

# href bandwidth for my dataset is about 690km - this is pretty big. Chapman pg. 44 - honestly, subjective choice can 
# be best. "A natural method for choosing the smoothing parameter is to plot out several curves and choose the estimate that is most
# in accordance with one's prior ideas about the density". 

# Furthermore on this point, "The reference bandwidth supposes that
# the UD is a bivariate normal distribution, which is disputable in most ecological
# studies. When the animals uses several centers of activity, the reference smoothing parameter is often too large, which results 
# into a strong oversmoothing of the data (the estimated UD predicts the frequent presence of the animal in areas
# which are not actually used)" - adehabitatHR vignette

# In general, home range analyses are rarely documented with enough detail to reproduce them or to allow for comparative studies. 
# At a minimum, h should be specified and it should be mentioned what method was used to calculate it. 

# Setup -------------------------------------------------------------------

# packages
library(dplyr)
library(maptools)
library(rgdal) # probably not needed, this is outdated by stars, terra, and sf
# library(GeoLocTools)
# setupGeolocation()
# Not available for R>=4.0, though you can get around that if you need
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
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/")

load("./LAALdata_midway_withTrackID.Rdata")
LAALmid <- LAAL
LAALmid$id <- paste0("lm",LAALmid$track)
load("./LAALdata_tern_withTrackID.Rdata")
LAALtern <- LAAL
LAALtern$id <- paste0("lt",LAALtern$track)

LAAL <- rbind(LAALmid, LAALtern)

load("./BFALdata_midway_withTrackID.Rdata")
BFALmid <- BFAL
BFALmid$id <- paste0("bm",BFALmid$track)
load("./BFALdata_tern_withTrackID.Rdata")
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

# Calculate grid value ----------------------------------------------------

all_data_1<-all_data[,1:3]
sp::coordinates(all_data_1) <- c("x", "y")
proj4string(all_data_1) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # placeholder
all_data_1 <- spTransform(all_data_1,lcea)
grid_input <- calculate_sp_obj_extent(all_data_1,0.1)

# Effect of grid size, however, is negligible...
# http://r-sig-geo.2731867.n2.nabble.com/adehabitatHR-kerneloverlap-and-kernelUD-etc-td7586144.html

# https://trac.faunalia.it/animove/wiki/AnimoveHowto#Whatunitisthegridcellsizein
# some more good notes on parameters: 
  # https://lists.faunalia.it/pipermail/animov/2006-May/000137.html


# Generate individual KDE -------------------------------------------------

        #### FOR TESTING PURPOSES: COMBINE ALL SPP INTO 1 KDE, ELIMINATE THE AVERAGING STEP LATER: 
        # load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_midway_withTrackID.Rdata")
        # LAALmid <- LAAL
        # LAALmid$id <- paste0("lm")
        # load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_tern_withTrackID.Rdata")
        # LAALtern <- LAAL
        # LAALtern$id <- paste0("lt")
        # 
        # LAAL <- rbind(LAALmid, LAALtern)
        # 
        # load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_midway_withTrackID.Rdata")
        # BFALmid <- BFAL
        # BFALmid$id <- paste0("bm")
        # load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_tern_withTrackID.Rdata")
        # BFALtern <- BFAL
        # BFALtern$id <- paste0("bt")
        # 
        # BFAL <- rbind(BFALmid, BFALtern)
        # 
        # all_data <- rbind(LAAL,BFAL)
        # lm <- all_data[grep("lm", all_data$id), ]
        # lt <- all_data[grep("lt", all_data$id), ]
        # bm <- all_data[grep("bm", all_data$id), ]
        # bt <- all_data[grep("bt", all_data$id), ]
        
        
        # CALC GRID VALUE #
        # all_data_1<-all_data[,1:3]
        # sp::coordinates(all_data_1) <- c("x", "y")
        # proj4string(all_data_1) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs") # placeholder
        # all_data_1 <- spTransform(all_data_1,lcea)
        # grid_input <- calculate_sp_obj_extent(all_data_1)
        # 
        # bound <- structure(list(x = barrier$x, y = barrier$y, .Names = c("x", "y")))
        # 
        # bound <- cbind(barrier$x,barrier$y)
        # Slo1 <- Line(bound)
        # Sli1 <- Lines(list(Slo1), ID="frontier1")
        # barrier <- SpatialLines(list(Sli1))
        # test <- kernelUD(all_data_1, grid=grid_input,same4all=T,extent=0.1, h=150000, boundary = barrier)
        # turning angles > pi/2
        
# kernelUD output currently is probability density: absolute density (events per unit area) divided by total number of events
results <- kernelUD(all_data_1, grid=grid_input,same4all=T,extent=0.1,h=150000)# grid-input guarantees AT LEAST 300km x 300km - in actuality, it is very slightly more than 300km x 300km.
image(results[[1]])
# Very very important: 
# https://meridian.allenpress.com/copeia/article/2006/4/797/260240
# Large H is ok; pick what is biologically relevant
# More discussion of h and grid parameter: 
# https://lists.faunalia.it/pipermail/animov/2006-May/000137.html
# GRID JUST CHANGES HOW THE RASTER LOOKS - THE SEARCH WINDOW/BANDWIDTH/K IS WHAT CONTROLS BIOLOGICAL RELEVANCE

names <- names(results)
lm_kde <- results[grep("lm",names)]
lt_kde <- results[grep("lt",names)]
bm_kde <- results[grep("bm",names)]
bt_kde <- results[grep("bt",names)]
all_LAAL <- append(lm_kde, lt_kde)
all_BFAL <- append(bm_kde, bt_kde)

# reminder on raster functions: not used becuase lines 176-242 are functionally equivalent
# ud1 <- stack(lapply(lm_kde, raster))
# plot(raster(lm_kde[[1]]))


# Normalization -> values range between 0 and 1 ---------------------------

# Normalize each raster, then average those rasters in the next step. 
# Normalizing to account for sample size differences in data points. 
# Unless there are two identical max values in a UD, there will be only 1 
# entry that equals 1 after normalizing. 

# THIS DOES NOT NEED TO BE DONE FOR PROBABILITY DENSITIES, WHICH ARE ALREADY NORMALIZED

# normalize lm #
# for (i in 1:length(lm_kde)){
#   print(length(lm_kde[[i]]$ud[lm_kde[[i]]$ud!=0])) # how many entries don't equal 0? 
#   lm_holder <- lm_kde[[i]]
#   lm_holder$ud[lm_holder$ud==0] <- NA
#   min_value <- min(lm_holder$ud[!is.na(lm_holder$ud)])
#   max_value <- max(lm_holder$ud[!is.na(lm_holder$ud)])
#   lm_holder$ud <- sapply(lm_holder$ud, function(x) (x-min_value)/(max_value-min_value))
#   lm_kde[[i]]$ud <- lm_holder$ud 
#   print(length(lm_kde[[1]]$ud[lm_kde[[1]]$ud!=0]))
#   print(i)
# }

# normalize lt #
# for (i in 1:length(lt_kde)){
#   lt_holder <- lt_kde[[i]]
#   lt_holder$ud[lt_holder$ud==0] <- NA
#   min_value <- min(lt_holder$ud[!is.na(lt_holder$ud)])
#   max_value <- max(lt_holder$ud[!is.na(lt_holder$ud)])
#   lt_holder$ud <- sapply(lt_holder$ud, function(x) (x-min_value)/(max_value-min_value))
#   lt_kde[[i]]$ud <- lt_holder$ud 
# }

# normalize bm #
# for (i in 1:length(bm_kde)){
#   bm_holder <- bm_kde[[i]]
#   bm_holder$ud[bm_holder$ud==0] <- NA
#   min_value <- min(bm_holder$ud[!is.na(bm_holder$ud)])
#   max_value <- max(bm_holder$ud[!is.na(bm_holder$ud)])
#   bm_holder$ud <- sapply(bm_holder$ud, function(x) (x-min_value)/(max_value-min_value))
#   bm_kde[[i]]$ud <- bm_holder$ud 
# }

# normalize bt #
# for (i in 1:length(bt_kde)){
#   bt_holder <- bt_kde[[i]]
#   bt_holder$ud[bt_holder$ud==0] <- NA
#   min_value <- min(bt_holder$ud[!is.na(bt_holder$ud)])
#   max_value <- max(bt_holder$ud[!is.na(bt_holder$ud)])
#   bt_holder$ud <- sapply(bt_holder$ud, function(x) (x-min_value)/(max_value-min_value))
#   bt_kde[[i]]$ud <- bt_holder$ud 
# }


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
    
# reminder for converting into a spatial pixel dataframe if needed: 
  # udspdf <- estUDm2spixdf(lm_kde)

allLAAL <- all_LAAL_averaged_estUD
allBFAL <- all_BFAL_averaged_estUD
midLAAL <- lm_averaged_estUD
midBFAL <- bm_averaged_estUD
ternLAAL <- lt_averaged_estUD
ternBFAL <- bt_averaged_estUD

path <- paste0(getwd(),"/final_push/final_ud/")
# save(allLAAL, file=paste0(path,"allLAAL.Rdata"))
# save(allBFAL, file=paste0(path,"allBFAL.Rdata"))
# save(midLAAL, file=paste0(path,"midLAAL.Rdata"))
# save(midBFAL, file=paste0(path,"midBFAL.Rdata"))
# save(ternLAAL, file=paste0(path,"ternLAAL.Rdata"))
# save(ternBFAL, file=paste0(path,"ternBFAL.Rdata"))

### load in files ###

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/final_ud/")
path <- getwd()
file_list <- list.files(path=path)

for (i in 1:length(file_list)){
  load(file_list[i])
}

### done loading ###

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
# save(allLAAL_v_allBFAL, file=paste0(path,"allLAAL_v_allBFAL.Rdata"))
# save(midLAAL_v_ternLAAL, file=paste0(path,"midLAAL_v_ternLAAL.Rdata"))
# save(midBFAL_v_ternBFAL, file=paste0(path,"midBFAL_v_ternBFAL.Rdata"))
# save(midLAAL_v_midBFAL, file=paste0(path,"midLAAL_v_midBFAL.Rdata"))
# save(ternLAAL_v_ternBFAL, file=paste0(path,"ternLAAL_v_ternBFAL.Rdata"))


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

# Rasters/Contours and GGplot visualization -------------------------------

### allLAAL v allBFAL ###

#95% contours and 100% UD RASTERS
# LAAL 
image(allLAAL)
plot(getverticeshr(allLAAL,percent=95), add=T)
vert95_allLAAL <- getverticeshr(allLAAL,percent=95)
gg95_allLAAL <- fortify(vert95_allLAAL)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_contours/vert95_allLAAL.Rdata"
save(vert95_allLAAL,file=path)

# for plotting UD raster in overallLAALBFAL_mapping: 
allLAAL.ud.vol <- getvolumeUD(allLAAL, standardize=TRUE)
plot(allLAAL.ud.vol)
allLAAL.ud.vol.raster <- raster(allLAAL.ud.vol)
plot(allLAAL.ud.vol.raster)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_rasters/allLAAL_ud_vol_rast.Rdata"
save(allLAAL.ud.vol.raster,file=path)

# BFAL 
image(allBFAL)
plot(getverticeshr(allBFAL,percent=95), add=T)
vert95_allBFAL <- getverticeshr(allBFAL,percent=95)
gg95_allBFAL <- fortify(vert95_allBFAL)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_contours/vert95_allBFAL.Rdata"
save(vert95_allBFAL,file=path)

# for plotting UD raster in overallLAALBFAL_mapping: 
allBFAL.ud.vol <- getvolumeUD(allBFAL, standardize=TRUE)
plot(allBFAL.ud.vol)
allBFAL.ud.vol.raster <- raster(allBFAL.ud.vol)
plot(allBFAL.ud.vol.raster)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_rasters/allBFAL_ud_vol_rast.Rdata"
save(allBFAL.ud.vol.raster,file=path)


# compare contours at 95UD: 
ggplot()+
  geom_polygon(data=gg95_allLAAL,aes(x=long,y=lat,group=group), color="red",fill=NA)+
  geom_polygon(data=gg95_allBFAL,aes(x=long,y=lat,group=group), color="blue",fill=NA)+
  ggtitle("95UD allLAAL (red), 95UD allBFAL (blue) - \n averaged and not normalized") +
  annotate("label", x = 3800000, y = 6000000, label = "UDOI = 0.432, BA = 0.574", cex=2) +
  coord_equal()


#50
# LAAL 
image(allLAAL)
plot(getverticeshr(allLAAL,percent=50), add=T)
vert50_allLAAL <- getverticeshr(allLAAL,percent=50)
gg50_allLAAL <- fortify(vert50_allLAAL)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_contours/vert50_allLAAL.Rdata"
save(vert50_allLAAL,file=path)

# BFAL
image(allBFAL)
plot(getverticeshr(allBFAL,percent=50), add=T)
vert50_allBFAL <- getverticeshr(allBFAL,percent=50)
gg50_allBFAL <- fortify(vert50_allBFAL)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_contours/vert50_allBFAL.Rdata"
save(vert50_allBFAL,file=path)

# compare contours at 50UD
ggplot()+
  geom_polygon(data=gg50_allLAAL,aes(x=long,y=lat,group=group), color="red",fill=NA)+
  geom_polygon(data=gg50_allBFAL,aes(x=long,y=lat,group=group), color="blue",fill=NA)+
  ggtitle("50UD allLAAL (red), 50UD allBFAL (blue) - \n averaged and not normalized") +
  annotate("label", x = 3800000, y = 6000000, label = "UDOI = 0.0005, BA = 0.024", cex=2) +
  coord_equal()

#10
# LAAL 
image(allLAAL)
plot(getverticeshr(allLAAL,percent=10), add=T)
vert10_allLAAL <- getverticeshr(allLAAL,percent=10)
gg10_allLAAL <- fortify(vert10_allLAAL)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_contours/vert10_allLAAL.Rdata"
save(vert10_allLAAL,file=path)

# BFAL
image(allBFAL)
plot(getverticeshr(allBFAL,percent=10), add=T)
vert10_allBFAL <- getverticeshr(allBFAL,percent=10)
gg10_allBFAL <- fortify(vert10_allBFAL)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_contours/vert10_allBFAL.Rdata"
save(vert10_allBFAL,file=path)

### midLAAL ###

#95-50-10% contours
image(midLAAL)
vert95_midLAAL <- getverticeshr(midLAAL,percent=95)
vert50_midLAAL <- getverticeshr(midLAAL,percent=50)
vert10_midLAAL <- getverticeshr(midLAAL,percent=10)
plot(getverticeshr(midLAAL,percent=95), add=T)
plot(getverticeshr(midLAAL,percent=50), add=T)
plot(getverticeshr(midLAAL,percent=10), add=T)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/midLAAL/master_script_contours/"
save(vert95_midLAAL,file=paste0(path,"vert95_midLAAL.Rdata"))
save(vert50_midLAAL,file=paste0(path,"vert50_midLAAL.Rdata"))
save(vert10_midLAAL,file=paste0(path,"vert10_midLAAL.Rdata"))

# for plotting UD raster in midwayLAAL_mapping: 
midLAAL.ud.vol <- getvolumeUD(midLAAL, standardize=TRUE)
midLAAL.ud.vol.raster <- raster(midLAAL.ud.vol)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/midLAAL/master_script_rasters/"
save(midLAAL.ud.vol.raster,file=paste0(path,"midLAAL_ud_vol_rast.Rdata"))

### midBFAL ###

#95-50-10% contours
image(midBFAL)
vert95_midBFAL <- getverticeshr(midBFAL,percent=95)
vert50_midBFAL <- getverticeshr(midBFAL,percent=50)
vert10_midBFAL <- getverticeshr(midBFAL,percent=10)
plot(getverticeshr(midBFAL,percent=95), add=T)
plot(getverticeshr(midBFAL,percent=50), add=T)
plot(getverticeshr(midBFAL,percent=10), add=T)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/midBFAL/master_script_contours/"
save(vert95_midBFAL,file=paste0(path,"vert95_midBFAL.Rdata"))
save(vert50_midBFAL,file=paste0(path,"vert50_midBFAL.Rdata"))
save(vert10_midBFAL,file=paste0(path,"vert10_midBFAL.Rdata"))

# for plotting UD raster in midwayLAAL_mapping: 
midBFAL.ud.vol <- getvolumeUD(midBFAL, standardize=TRUE)
midBFAL.ud.vol.raster <- raster(midBFAL.ud.vol)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/midBFAL/master_script_rasters/"
save(midBFAL.ud.vol.raster,file=paste0(path,"midBFAL_ud_vol_rast.Rdata"))

### ternLAAL ###

#95-50-10% contours
image(ternLAAL)
vert95_ternLAAL <- getverticeshr(ternLAAL,percent=95)
vert50_ternLAAL <- getverticeshr(ternLAAL,percent=50)
vert10_ternLAAL <- getverticeshr(ternLAAL,percent=10)
plot(getverticeshr(ternLAAL,percent=95), add=T)
plot(getverticeshr(ternLAAL,percent=50), add=T)
plot(getverticeshr(ternLAAL,percent=10), add=T)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternLAAL/master_script_contours/"
save(vert95_ternLAAL,file=paste0(path,"vert95_ternLAAL.Rdata"))
save(vert50_ternLAAL,file=paste0(path,"vert50_ternLAAL.Rdata"))
save(vert10_ternLAAL,file=paste0(path,"vert10_ternLAAL.Rdata"))

# for plotting UD raster in midwayLAAL_mapping: 
ternLAAL.ud.vol <- getvolumeUD(ternLAAL, standardize=TRUE)
ternLAAL.ud.vol.raster <- raster(ternLAAL.ud.vol)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternLAAL/master_script_rasters/"
save(ternLAAL.ud.vol.raster,file=paste0(path,"ternLAAL_ud_vol_rast.Rdata"))

### ternBFAL ###

#95-50-10% contours
image(ternBFAL)
vert95_ternBFAL <- getverticeshr(ternBFAL,percent=95)
vert50_ternBFAL <- getverticeshr(ternBFAL,percent=50)
vert10_ternBFAL <- getverticeshr(ternBFAL,percent=10)
plot(getverticeshr(ternBFAL,percent=95), add=T)
plot(getverticeshr(ternBFAL,percent=50), add=T)
plot(getverticeshr(ternBFAL,percent=10), add=T)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternBFAL/master_script_contours/"
save(vert95_ternBFAL,file=paste0(path,"vert95_ternBFAL.Rdata"))
save(vert50_ternBFAL,file=paste0(path,"vert50_ternBFAL.Rdata"))
save(vert10_ternBFAL,file=paste0(path,"vert10_ternBFAL.Rdata"))

# for plotting UD raster in midwayLAAL_mapping: 
ternBFAL.ud.vol <- getvolumeUD(ternBFAL, standardize=TRUE)
ternBFAL.ud.vol.raster <- raster(ternBFAL.ud.vol)
path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternBFAL/master_script_rasters/"
save(ternBFAL.ud.vol.raster,file=paste0(path,"ternBFAL_ud_vol_rast.Rdata"))

#>>> DONE WITH MAPPING CONTOURS AND RASTERS <<< 

##### Home range in raster mode
# GUIDE: https://ecosystems.psu.edu/research/labs/walter-lab/manual/home-range-estimation/link-to-pdf
# file:///Users/dallasjordan/Downloads/Chapter04_2016.pdf starting bottom of page 80
#the value of a pixel is equal to the percentage of the smallest home range containing this pixel:

volumeLAAL <- getvolumeUD(all_LAAL_averaged_estUD)
image(volumeLAAL)
plot(getverticeshr(allLAAL,percent=95), add=T)

volumeBFAL <- getvolumeUD(all_BFAL_averaged_estUD)
image(volumeBFAL)
plot(getverticeshr(allBFAL,percent=95), add=T)

allLAAL.ud.vol.raster <- raster(volumeLAAL)
plot.new()
breaks <- c(0, 50, 80, 90, 95)
plot(volumeLAAL, col=alpha(viridis(4, direction=-1),0.4), breaks=breaks,interpolate=TRUE,legend.shrink=0.80,legend.args=list(text="UD by
Volume (%)",side=4, font=2, line=2.5, cex=0.8))
plot(volumeBFAL, col=alpha(plasma(4),0.4), breaks=breaks,interpolate=TRUE,
     main="Kernel Density Estimation, Plug-in Bandwidth",xlab="Coords X",
     ylab="Coords Y", legend.shrink=0.80,legend.args=list(text="UD by
Volume (%)",side=4, font=2, line=2.5, cex=0.8))


allLAAL.50vol <- getverticeshr(allLAAL, percent = 50,ida = NULL, unin = "m",
                              unout = "ha", standardize=TRUE)
plot(allLAAL.50vol,add=T)

allBFAL.50vol <- getverticeshr(allBFAL, percent = 50,ida = NULL, unin = "m",
                               unout = "ha", standardize=TRUE)
plot(allBFAL.50vol,add=T)

### midLAAL v ternLAAL 

#95 
image(midLAAL)
plot(getverticeshr(midLAAL,percent=95), add=T)
vert95_midLAAL <- getverticeshr(midLAAL,percent=95)
gg95_midLAAL <- fortify(vert95_midLAAL)


image(ternLAAL)
plot(getverticeshr(ternLAAL,percent=95), add=T)
vert95_ternLAAL <- getverticeshr(ternLAAL,percent=95)
gg95_ternLAAL <- fortify(vert95_ternLAAL)
ggplot()+
  geom_polygon(data=gg95_midLAAL,aes(x=long,y=lat,group=group), color="red",fill=NA)+
  geom_polygon(data=gg95_ternLAAL,aes(x=long,y=lat,group=group), color="blue",fill=NA)+
  ggtitle("95UD midLAAL (red), 95UD ternLAAL (blue) - \n averaged and not normalized") +
  annotate("label", x = 3800000, y = 6000000, label = "UDOI = 0.779, BA = 0.772", cex=2) +
  coord_equal()

#50 

image(midLAAL)
plot(getverticeshr(midLAAL,percent=50), add=T)
vert50_midLAAL <- getverticeshr(midLAAL,percent=50)
gg50_midLAAL <- fortify(vert50_midLAAL)

image(ternLAAL)
plot(getverticeshr(ternLAAL,percent=50), add=T)
vert50_ternLAAL <- getverticeshr(ternLAAL,percent=50)
gg50_ternLAAL <- fortify(vert50_ternLAAL)
ggplot()+
  geom_polygon(data=gg50_midLAAL,aes(x=long,y=lat,group=group), color="red",fill=NA)+
  geom_polygon(data=gg50_ternLAAL,aes(x=long,y=lat,group=group), color="blue",fill=NA)+
  ggtitle("50UD midLAAL (red), 50UD ternLAAL (blue) - \n averaged and not normalized") +
  annotate("label", x = 1200000, y = 5400000, label = "UDOI = 0.114, HR= 0.45") +
  coord_cartesian(xlim = c(-4e+06, 2e+06), ylim = c(4000000, 5500000))

# Randomization tests -----------------------------------------------------

path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/permutation_tests/test_stats/"
setwd(path)
load("all_test_stats.RData")

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/permutation_tests/data_to_load/")
path <- getwd()
file_list <- list.files(path=path)

for (i in 1:length(file_list)){
  load(file_list[i])
}

#### midLAAL v ternLAAL ###
  
  # 95th/50th isopleth and BA
  iter <- 1000
  significance_tallyUDOI_95 <- 0
  significance_tallyUDOI_50 <- 0
  significance_tallyBA_95 <- 0
  significance_tallyBA_50 <- 0
  iter_tally <- 0
  resultsUDOI_95_storage <- vector()
  resultsUDOI_50_storage <- vector()
  resultsBA_95_storage <- vector()
  resultsBA_50_storage <- vector()
  
  # Permute the cell values rather than the track labels (permuting spatially instead of by labels)
      # take normalized and averaged estUD
      # random sample order of grid cells
      # recalculate overlap metrics of these randomized UDs
      # repeat! simple. 
  
  # New take: 
      # see permutation_XX_XX.R scripts for track-ID permutation
  
  # Begin permutation test
  
  for (i in 1:iter){
    ##### Randomize order of normalized/averaged UD grid cells
    # This is doing it over the whole grid - is that right?
    midLAAL_iter <- midLAAL_v_ternLAAL[[1]]
    midLAAL_randomized_indices <- sample(nrow(midLAAL_v_ternLAAL[[1]]@data))
    midLAAL_iter@data$ud <- midLAAL_iter@data[midLAAL_randomized_indices,]
    
    ternLAAL_iter <- midLAAL_v_ternLAAL[[2]]
    ternLAAL_randomized_indices <- sample(nrow(midLAAL_v_ternLAAL[[2]]@data))
    ternLAAL_iter@data$ud <- ternLAAL_iter@data[ternLAAL_randomized_indices,]
    
    ##### Create estUDm ##### 
    imLAAL<-midLAAL_iter
    itLAAL<-ternLAAL_iter
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
  
  
  p_valueUDOI_95 <- significance_tallyUDOI_95/iter
  mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
  sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)
  
  p_valueUDOI_95 # 1
  mean_valueUDOI_95 # 0.023
  sd_valueUDOI_95 # 0.011
  
  p_valueUDOI_50 <- significance_tallyUDOI_50/iter
  mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
  sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)
  
  p_valueUDOI_50 # 1
  mean_valueUDOI_50 # 0.0005
  sd_valueUDOI_50 # 0.0008
  
  p_valueBA_95 <- significance_tallyBA_95/iter
  mean_valueBA_95 <- mean(resultsBA_95_storage)
  sd_valueBA_95 <- sd(resultsBA_95_storage)
  
  p_valueBA_95 # 1
  mean_valueBA_95 # 0.122
  sd_valueBA_95 # 0.028
  
  p_valueBA_50 <- significance_tallyBA_50/iter
  mean_valueBA_50 <- mean(resultsBA_50_storage)
  sd_valueBA_50 <- sd(resultsBA_50_storage)
  
  p_valueBA_50 # 1
  mean_valueBA_50 # 0.015
  sd_valueBA_50 # 0.016
  
  # book keeping - track counts, do not need to run when testing # 
  add_lm <-lm
  add_lt <-lt
  add_lm$island <- "LAALmid"
  add_lt$island <- "LAALtern"
  combined_lm_lt <- rbind(add_lm,add_lt)
  by_island <- split(combined_lm_lt, combined_lm_lt$island)
  midLAAL_tracks <- split(by_island[[1]],by_island[[1]]$track)
  ternLAAL_tracks <- split(by_island[[2]],by_island[[2]]$track)
  all_tracks <- c(midLAAL_tracks,ternLAAL_tracks)
  available_mid_tracks <- unique(by_island[[1]]$track)
  available_tern_tracks <- unique(by_island[[2]]$track)
  n_amt <- length(available_mid_tracks) # you are calculating these to preserve sample size in randomizations
  n_att <- length(available_tern_tracks)
  n_all <- n_amt+n_att
  counts <- combined_lm_lt %>% count(id,track)
  
  
  
  
  
#### midBFAL v ternBFAL ###
  
  # 95th/50th isopleth and BA
  iter <- 1000
  significance_tallyUDOI_95 <- 0
  significance_tallyUDOI_50 <- 0
  significance_tallyBA_95 <- 0
  significance_tallyBA_50 <- 0
  iter_tally <- 0
  resultsUDOI_95_storage <- vector()
  resultsUDOI_50_storage <- vector()
  resultsBA_95_storage <- vector()
  resultsBA_50_storage <- vector()
  
  for (i in 1:iter){
    ##### Randomize order of normalized/averaged UD grid cells
    # This is doing it over the whole grid - is that right?
    midBFAL_iter <- midBFAL_v_ternBFAL[[1]]
    midBFAL_randomized_indices <- sample(nrow(midBFAL_v_ternBFAL[[1]]@data))
    midBFAL_iter@data$ud <- midBFAL_iter@data[midBFAL_randomized_indices,]
    
    ternBFAL_iter <- midBFAL_v_ternBFAL[[2]]
    ternBFAL_randomized_indices <- sample(nrow(midBFAL_v_ternBFAL[[2]]@data))
    ternBFAL_iter@data$ud <- ternBFAL_iter@data[ternBFAL_randomized_indices,]
    
    ##### Create estUDm ##### 
    imBFAL<-midBFAL_iter
    itBFAL<-ternBFAL_iter
    imBFAL_v_itBFAL <- list(imBFAL,itBFAL)
    class(imBFAL_v_itBFAL)<-"estUDm"
    names(imBFAL_v_itBFAL)<-c("imBFAL","itBFAL")
    
    iter_mbtb_UDOI_95 <- kerneloverlaphr(imBFAL_v_itBFAL, method="UDOI", percent=95, conditional=T)
    iter_mbtb_UDOI_50 <- kerneloverlaphr(imBFAL_v_itBFAL, method="UDOI", percent=50, conditional=T)
    iter_mbtb_BA_95 <- kerneloverlaphr(imBFAL_v_itBFAL, method="BA", percent=95, conditional=T) 
    iter_mbtb_BA_50 <- kerneloverlaphr(imBFAL_v_itBFAL, method="BA", percent=50, conditional=T)
    
    
    iter_tally <- iter_tally+1
    print(iter_tally)
    
    resultsUDOI_95_storage[i] <- iter_mbtb_UDOI_95[1,2]
    print(resultsUDOI_95_storage[i])
    if (iter_mbtb_UDOI_95[1,2]<mbtb_UDOI_95_test_stat){
      print(iter_mbtb_UDOI_95)
      print("+1!")
      significance_tallyUDOI_95 <- significance_tallyUDOI_95+1
    }
    
    resultsUDOI_50_storage[i] <- iter_mbtb_UDOI_50[1,2]
    print(resultsUDOI_50_storage[i])
    if (iter_mbtb_UDOI_50[1,2]<mbtb_UDOI_50_test_stat){
      print(iter_mbtb_UDOI_50)
      print("+1!")
      significance_tallyUDOI_50 <- significance_tallyUDOI_50+1
    }
    
    resultsBA_95_storage[i] <- iter_mbtb_BA_95[1,2]
    print(resultsBA_95_storage[i])
    if (iter_mbtb_BA_95[1,2]<mbtb_95_BA_test_stat){
      print(iter_mbtb_BA_95)
      print("+1!")
      significance_tallyBA_95 <- significance_tallyBA_95+1
    }
    
    resultsBA_50_storage[i] <- iter_mbtb_BA_50[1,2]
    print(resultsBA_50_storage[i])
    if (iter_mbtb_BA_50[1,2]<mbtb_50_BA_test_stat){
      print(iter_mbtb_BA_50)
      print("+1!")
      significance_tallyBA_50 <- significance_tallyBA_50+1
    }
  }
  
  
  p_valueUDOI_95 <- significance_tallyUDOI_95/iter
  mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
  sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)
  
  p_valueUDOI_95 # 1
  mean_valueUDOI_95 # 0.049
  sd_valueUDOI_95 # 0.017
  
  p_valueUDOI_50 <- significance_tallyUDOI_50/iter
  mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
  sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)
  
  p_valueUDOI_50 # 0.913
  mean_valueUDOI_50 # 0.001
  sd_valueUDOI_50 # 0.001
  
  p_valueBA_95 <- significance_tallyBA_95/iter
  mean_valueBA_95 <- mean(resultsBA_95_storage)
  sd_valueBA_95 <- sd(resultsBA_95_storage)
  
  p_valueBA_95 # 1
  mean_valueBA_95 # 0.183
  sd_valueBA_95 # 0.026
  
  p_valueBA_50 <- significance_tallyBA_50/iter
  mean_valueBA_50 <- mean(resultsBA_50_storage)
  sd_valueBA_50 <- sd(resultsBA_50_storage)
  
  p_valueBA_50 # 0.921
  mean_valueBA_50 # 0.022
  sd_valueBA_50 # 0.016
  

  #### midLAAL v midBFAL ###
  
  # 95th/50th isopleth and BA
  iter <- 1000
  significance_tallyUDOI_95 <- 0
  significance_tallyUDOI_50 <- 0
  significance_tallyBA_95 <- 0
  significance_tallyBA_50 <- 0
  iter_tally <- 0
  resultsUDOI_95_storage <- vector()
  resultsUDOI_50_storage <- vector()
  resultsBA_95_storage <- vector()
  resultsBA_50_storage <- vector()
  
  for (i in 1:iter){
    ##### Randomize order of normalized/averaged UD grid cells
    # This is doing it over the whole grid - is that right?
    midLAAL_iter <- midLAAL_v_midBFAL[[1]]
    midLAAL_randomized_indices <- sample(nrow(midLAAL_v_midBFAL[[1]]@data))
    midLAAL_iter@data$ud <- midLAAL_iter@data[midLAAL_randomized_indices,]
    
    midBFAL_iter <- midLAAL_v_midBFAL[[2]]
    midBFAL_randomized_indices <- sample(nrow(midLAAL_v_midBFAL[[2]]@data))
    midBFAL_iter@data$ud <- midBFAL_iter@data[midBFAL_randomized_indices,]
    
    ##### Create estUDm ##### 
    imLAAL<-midLAAL_iter
    imBFAL<-midBFAL_iter
    imLAAL_v_imBFAL <- list(imLAAL,imBFAL)
    class(imLAAL_v_imBFAL)<-"estUDm"
    names(imLAAL_v_imBFAL)<-c("imLAAL","imBFAL")
    
    iter_mlmb_UDOI_95 <- kerneloverlaphr(imLAAL_v_imBFAL, method="UDOI", percent=95, conditional=T)
    iter_mlmb_UDOI_50 <- kerneloverlaphr(imLAAL_v_imBFAL, method="UDOI", percent=50, conditional=T)
    iter_mlmb_BA_95 <- kerneloverlaphr(imLAAL_v_imBFAL, method="BA", percent=95, conditional=T) 
    iter_mlmb_BA_50 <- kerneloverlaphr(imLAAL_v_imBFAL, method="BA", percent=50, conditional=T)
    
    
    iter_tally <- iter_tally+1
    print(iter_tally)
    
    resultsUDOI_95_storage[i] <- iter_mlmb_UDOI_95[1,2]
    print(resultsUDOI_95_storage[i])
    if (iter_mlmb_UDOI_95[1,2]<mlmb_UDOI_95_test_stat){
      print(iter_mlmb_UDOI_95)
      print("+1!")
      significance_tallyUDOI_95 <- significance_tallyUDOI_95+1
    }
    
    resultsUDOI_50_storage[i] <- iter_mlmb_UDOI_50[1,2]
    print(resultsUDOI_50_storage[i])
    if (iter_mlmb_UDOI_50[1,2]<mlmb_UDOI_50_test_stat){
      print(iter_mlmb_UDOI_50)
      print("+1!")
      significance_tallyUDOI_50 <- significance_tallyUDOI_50+1
    }
    
    resultsBA_95_storage[i] <- iter_mlmb_BA_95[1,2]
    print(resultsBA_95_storage[i])
    if (iter_mlmb_BA_95[1,2]<mlmb_95_BA_test_stat){
      print(iter_mlmb_BA_95)
      print("+1!")
      significance_tallyBA_95 <- significance_tallyBA_95+1
    }
    
    resultsBA_50_storage[i] <- iter_mlmb_BA_50[1,2]
    print(resultsBA_50_storage[i])
    if (iter_mlmb_BA_50[1,2]<mlmb_50_BA_test_stat){
      print(iter_mlmb_BA_50)
      print("+1!")
      significance_tallyBA_50 <- significance_tallyBA_50+1
    }
  }
  
  
  p_valueUDOI_95 <- significance_tallyUDOI_95/iter
  mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
  sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)
  
  p_valueUDOI_95 # 1
  mean_valueUDOI_95 # 0.032
  sd_valueUDOI_95 # 0.013
  
  p_valueUDOI_50 <- significance_tallyUDOI_50/iter
  mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
  sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)
  
  p_valueUDOI_50 # 0.98
  mean_valueUDOI_50 # 0.001
  sd_valueUDOI_50 # 0.001
  
  p_valueBA_95 <- significance_tallyBA_95/iter
  mean_valueBA_95 <- mean(resultsBA_95_storage)
  sd_valueBA_95 <- sd(resultsBA_95_storage)
  
  p_valueBA_95 # 1
  mean_valueBA_95 # 0.146
  sd_valueBA_95 # 0.028
  
  p_valueBA_50 <- significance_tallyBA_50/iter
  mean_valueBA_50 <- mean(resultsBA_50_storage)
  sd_valueBA_50 <- sd(resultsBA_50_storage)
  
  p_valueBA_50 # 0.987
  mean_valueBA_50 # 0.019
  sd_valueBA_50 # 0.016
  
  
  
  #### ternLAAL v ternBFAL ###
  
  # 95th/50th isopleth and BA
  iter <- 1000
  significance_tallyUDOI_95 <- 0
  significance_tallyUDOI_50 <- 0
  significance_tallyBA_95 <- 0
  significance_tallyBA_50 <- 0
  iter_tally <- 0
  resultsUDOI_95_storage <- vector()
  resultsUDOI_50_storage <- vector()
  resultsBA_95_storage <- vector()
  resultsBA_50_storage <- vector()
  
  for (i in 1:iter){
    ##### Randomize order of normalized/averaged UD grid cells
    # This is doing it over the whole grid - is that right?
    ternLAAL_iter <- ternLAAL_v_ternBFAL[[1]]
    ternLAAL_randomized_indices <- sample(nrow(ternLAAL_v_ternBFAL[[1]]@data))
    ternLAAL_iter@data$ud <- ternLAAL_iter@data[ternLAAL_randomized_indices,]
    
    ternBFAL_iter <- ternLAAL_v_ternBFAL[[2]]
    ternBFAL_randomized_indices <- sample(nrow(ternLAAL_v_ternBFAL[[2]]@data))
    ternBFAL_iter@data$ud <- ternBFAL_iter@data[ternBFAL_randomized_indices,]
    
    ##### Create estUDm ##### 
    itLAAL<-ternLAAL_iter
    itBFAL<-ternBFAL_iter
    itLAAL_v_itBFAL <- list(itLAAL,itBFAL)
    class(itLAAL_v_itBFAL)<-"estUDm"
    names(itLAAL_v_itBFAL)<-c("itLAAL","itBFAL")
    
    iter_tltb_UDOI_95 <- kerneloverlaphr(itLAAL_v_itBFAL, method="UDOI", percent=95, conditional=T)
    iter_tltb_UDOI_50 <- kerneloverlaphr(itLAAL_v_itBFAL, method="UDOI", percent=50, conditional=T)
    iter_tltb_BA_95 <- kerneloverlaphr(itLAAL_v_itBFAL, method="BA", percent=95, conditional=T) 
    iter_tltb_BA_50 <- kerneloverlaphr(itLAAL_v_itBFAL, method="BA", percent=50, conditional=T)
    
    
    iter_tally <- iter_tally+1
    print(iter_tally)
    
    resultsUDOI_95_storage[i] <- iter_tltb_UDOI_95[1,2]
    print(resultsUDOI_95_storage[i])
    if (iter_tltb_UDOI_95[1,2]<tltb_UDOI_95_test_stat){
      print(iter_tltb_UDOI_95)
      print("+1!")
      significance_tallyUDOI_95 <- significance_tallyUDOI_95+1
    }
    
    resultsUDOI_50_storage[i] <- iter_tltb_UDOI_50[1,2]
    print(resultsUDOI_50_storage[i])
    if (iter_tltb_UDOI_50[1,2]<tltb_UDOI_50_test_stat){
      print(iter_tltb_UDOI_50)
      print("+1!")
      significance_tallyUDOI_50 <- significance_tallyUDOI_50+1
    }
    
    resultsBA_95_storage[i] <- iter_tltb_BA_95[1,2]
    print(resultsBA_95_storage[i])
    if (iter_tltb_BA_95[1,2]<tltb_95_BA_test_stat){
      print(iter_tltb_BA_95)
      print("+1!")
      significance_tallyBA_95 <- significance_tallyBA_95+1
    }
    
    resultsBA_50_storage[i] <- iter_tltb_BA_50[1,2]
    print(resultsBA_50_storage[i])
    if (iter_tltb_BA_50[1,2]<tltb_50_BA_test_stat){
      print(iter_tltb_BA_50)
      print("+1!")
      significance_tallyBA_50 <- significance_tallyBA_50+1
    }
  }
  
  
  p_valueUDOI_95 <- significance_tallyUDOI_95/iter
  mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
  sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)
  
  p_valueUDOI_95 # 1
  mean_valueUDOI_95 # 0.036
  sd_valueUDOI_95 # 0.015
  
  p_valueUDOI_50 <- significance_tallyUDOI_50/iter
  mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
  sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)
  
  p_valueUDOI_50 # 0
  mean_valueUDOI_50 # 0.001
  sd_valueUDOI_50 # 0.001
  
  p_valueBA_95 <- significance_tallyBA_95/iter
  mean_valueBA_95 <- mean(resultsBA_95_storage)
  sd_valueBA_95 <- sd(resultsBA_95_storage)
  
  p_valueBA_95 # 1
  mean_valueBA_95 # 0.153
  sd_valueBA_95 # 0.027
  
  p_valueBA_50 <- significance_tallyBA_50/iter
  mean_valueBA_50 <- mean(resultsBA_50_storage)
  sd_valueBA_50 <- sd(resultsBA_50_storage)
  
  p_valueBA_50 # 0
  mean_valueBA_50 # 0.018
  sd_valueBA_50 # 0.017

  
  #### allLAAL v allBFAL ###
  
  # 95th/50th isopleth and BA
  iter <- 1000
  significance_tallyUDOI_95 <- 0
  significance_tallyUDOI_50 <- 0
  significance_tallyBA_95 <- 0
  significance_tallyBA_50 <- 0
  iter_tally <- 0
  resultsUDOI_95_storage <- vector()
  resultsUDOI_50_storage <- vector()
  resultsBA_95_storage <- vector()
  resultsBA_50_storage <- vector()
  
  for (i in 1:iter){
    ##### Randomize order of normalized/averaged UD grid cells
    # This is doing it over the whole grid - is that right?
    allLAAL_iter <- allLAAL_v_allBFAL[[1]]
    allLAAL_randomized_indices <- sample(nrow(allLAAL_v_allBFAL[[1]]@data))
    allLAAL_iter@data$ud <- allLAAL_iter@data[allLAAL_randomized_indices,]
    
    allBFAL_iter <- allLAAL_v_allBFAL[[2]]
    allBFAL_randomized_indices <- sample(nrow(allLAAL_v_allBFAL[[2]]@data))
    allBFAL_iter@data$ud <- allBFAL_iter@data[allBFAL_randomized_indices,]
    
    ##### Create estUDm ##### 
    iaLAAL<-allLAAL_iter
    iaBFAL<-allBFAL_iter
    iaLAAL_v_iaBFAL <- list(iaLAAL,iaBFAL)
    class(iaLAAL_v_iaBFAL)<-"estUDm"
    names(iaLAAL_v_iaBFAL)<-c("iaLAAL","iaBFAL")
    
    iter_alab_UDOI_95 <- kerneloverlaphr(iaLAAL_v_iaBFAL, method="UDOI", percent=95, conditional=T)
    iter_alab_UDOI_50 <- kerneloverlaphr(iaLAAL_v_iaBFAL, method="UDOI", percent=50, conditional=T)
    iter_alab_BA_95 <- kerneloverlaphr(iaLAAL_v_iaBFAL, method="BA", percent=95, conditional=T) 
    iter_alab_BA_50 <- kerneloverlaphr(iaLAAL_v_iaBFAL, method="BA", percent=50, conditional=T)
    
    
    iter_tally <- iter_tally+1
    print(iter_tally)
    
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
  
  
  p_valueUDOI_95 <- significance_tallyUDOI_95/iter
  mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
  sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)
  
  p_valueUDOI_95 # 1
  mean_valueUDOI_95 # 0.046
  sd_valueUDOI_95 # 0.015
  
  p_valueUDOI_50 <- significance_tallyUDOI_50/iter
  mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
  sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)
  
  p_valueUDOI_50 # 0.561
  mean_valueUDOI_50 # 0.001
  sd_valueUDOI_50 # 0.001
  
  p_valueBA_95 <- significance_tallyBA_95/iter
  mean_valueBA_95 <- mean(resultsBA_95_storage)
  sd_valueBA_95 <- sd(resultsBA_95_storage)
  
  p_valueBA_95 # 1
  mean_valueBA_95 # 0.181
  sd_valueBA_95 # 0.026
  
  p_valueBA_50 <- significance_tallyBA_50/iter
  mean_valueBA_50 <- mean(resultsBA_50_storage)
  sd_valueBA_50 <- sd(resultsBA_50_storage)
  
  p_valueBA_50 # 0.568
  mean_valueBA_50 # 0.024
  sd_valueBA_50 # 0.016
  
  
  
  
  
  