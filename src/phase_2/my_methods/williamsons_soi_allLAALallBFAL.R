# Calculation of Williamson's SOI - for allLAAL allBFAL
# This requires a seperate script because the setup is different than individual class rasters
# This requires a separate script because relative density is calculated - for each raster, the value of
# each cell is how many occurrences of a species occur in that cell divided by the total number of occurences. 
# These relative density values are then summed according to Williamson's SOI formula calculation
# Possible that each raster cell will have more than one relocation in it, but not sure - to do 
# this calculation, you will need to use the Raster function instead of your usual KDE formulation

# IMPORTANT - Most recent update 7/7/2021 - There are methods implemented here to load in kernel 
# UD rasters exported from my master_script. These are already averaged and since they are
# probability densities they do not need normalizing. Just load in and run the WSOI calculation
# after converting to a raster layer and ensuring that the 'm' cell count is correct. 


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

# functions and objects
calculate_sp_obj_extent <- function(sp_obj){
  x_range <- sp_obj@bbox[1,2] - sp_obj@bbox[1,1]
  y_range <- sp_obj@bbox[2,2] - sp_obj@bbox[2,1]
  x_range_km <- x_range/1000
  y_range_km <- y_range/1000
  print1<- paste0("the x-range in km is ",x_range_km)
  print2<- paste0("the y-range in km is ",y_range_km)
  grid_calc <- x_range_km/50
  print3<- paste0("for a grid size of 50km^2 use grid parameter ",grid_calc)
  print(print1)
  print(print2)
  print(print3)
  return(grid_calc)
}
pdc_mercator_proj<-"+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
lcea <- "+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs" 
years <- c("2008","2009","2010","2011","2012")

# Load data ---------------------------------------------------------------

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
al <- rbind(lm, lt)
ab <- rbind(bm,bt)



# Import kernelUD rasters -------------------------------------------------

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/final_ud/")
path <- getwd()
file_list <- list.files(path=path)

for (i in 1:length(file_list)){
  load(file_list[i])
}

# converting to raster layer from estUD 

allLAAL <- raster(as(allLAAL,"SpatialPixelsDataFrame"))
allBFAL <- raster(as(allBFAL,"SpatialPixelsDataFrame"))

# Can skip to WSOI calculations - probability density is already normalized and 
# the rasters you loaded in are averaged already from your master_script

rastAL <- allLAAL
rastAB <- allBFAL


# Make density grid rasters  ----------------------------------------------

# Averaging 

# Can easily average through some version of this: 
# ud1 <- stack(lapply(lm_kde, raster))
# plot(raster(lm_kde[[1]]))

lm_list <- split(lm, f=lm$id)
lt_list <- split(lt, f=lt$id)
bm_list <- split(bm, f=bm$id)
bt_list <- split(bt, f=bt$id)

al_list <- append(lm_list,lt_list)
ab_list <- append(bm_list,bt_list)

# now you need to set extent and empty raster, then lapply(lm_avg,rasterize) to create a lot 
# of rasters, then stack those rasters, then avg

all_data_1<-all_data[,1:3]
# set up an 'empty' raster, here via an extent object derived from your data
e <- extent(all_data_1)
# In the next line, setting ext and res lets you not specify other defaults
r <- raster(e, ext=e, res=c(2.7,2.7), crs="+proj=longlat +datum=WGS84" ) # 2.7 degrees is 300km x 300km 
# you need to provide a function 'fun' for when there are multiple points per cell

al_avg1 <- lapply(al_list,"[", c('x','y')) # grabs just x and y columns for each element of list
al_avg2 <- lapply(al_avg1, function(x) rasterize(x,r,fun='count'))
for (i in 1:length(al_avg2)){
  values(al_avg2[[i]]) <- values(al_avg2[[i]])/nrow(al_avg1[[i]])
}
al_stack <- raster::stack(al_avg2)
plot(al_stack)
al_averaged <- raster::calc(al_stack, mean, na.rm=T)
plot(al_averaged) # looks pretty similar to lm_averaged, but confirmed slightly different

ab_avg1 <- lapply(ab_list,"[", c('x','y')) # grabs just x and y columns for each element of list
ab_avg2 <- lapply(ab_avg1, function(x) rasterize(x,r,fun='count'))
for (i in 1:length(ab_avg2)){
  values(ab_avg2[[i]]) <- values(ab_avg2[[i]])/nrow(ab_avg1[[i]])
}
ab_stack <- raster::stack(ab_avg2)
plot(ab_stack)
ab_averaged <- raster::calc(ab_stack, mean, na.rm=T)
plot(ab_averaged)

rastAL <- al_averaged
rastAB <- ab_averaged


# Normalization (not needed for kernelUD rasters) -------------------------

# Normalizing to account for sample size differences in data points. 
# Unless there are two identical max values in a UD, there will be only 1 
# entry that equals 1 after normalizing. 

### normalize lm ###
print(length(values(rastAL)!=0)) # make sure its NA, not 0
al_holder <- values(rastAL)
min_value <- min(al_holder[!is.na(al_holder)])
max_value <- max(al_holder[!is.na(al_holder)])
al_holder <- sapply(al_holder[!is.na(al_holder)], function(x) (x-min_value)/(max_value-min_value))
length(rastAL[!is.na(values(rastAL))])
length(al_holder)
rastAL[!is.na(values(rastAL))]<-al_holder
plot(rastAL)

### normalize bm ### 
print(length(values(rastAB)!=0)) # how many entries don't equal 0? 
ab_holder <- values(rastAB)
min_value <- min(ab_holder[!is.na(ab_holder)])
max_value <- max(ab_holder[!is.na(ab_holder)])
ab_holder <- sapply(ab_holder[!is.na(ab_holder)], function(x) (x-min_value)/(max_value-min_value))
length(rastAB[!is.na(values(rastAB))])
length(ab_holder)
rastAB[!is.na(values(rastAB))]<-ab_holder
plot(rastAB)


# Calculate Williamson's SOI ----------------------------------------------

# AL x AB
oAL <- rastAL
oAB <- rastAB
oAL[is.na(oAL)]<-0
oAB[is.na(oAB)]<-0
m <- 893 # NUMBER OF CELLS - THIS CHANGES BASED ON IF YOU ARE USING DENSITY RASTERS OR KERNEL UD (estUD imports)!
num <- sum((oAL[]*oAB[])*m)
denom <- (sum(oAL[]))*(sum(oAB[]))
AL_AB_SOI <- num/denom # I get 2.801 when using estUD rasters
AL_AB_SOI

plot(oAL)
plot(oAB)

# I did this for the whole study region...from Julia's paper: 
# We first calculated SOI for the entire study region for both 
# years together. Next, we calculated SOI relative to proximity 
# to the 1000 m isobath, in three regions: one value for the region 
# 15 km inshore of the 1000 m isobath, one for the region 15 km 
# offshore of the 1000 m isobath, and one for any region more than 
# 15 km from the 1000 m isobath. The 15 km distance was chosen to 
# match the spatial resolution of the density grids as described in 
# Section 2.5. To assess statistical significance, we calculated SOI
# and compared the observed value to a test distribution of 4999 SOI 
# values obtained by iterating randomized pilot whale and longline 
# density grids (Garrison et al., 2000; Harden and Williard, 2012).






