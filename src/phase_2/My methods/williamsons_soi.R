# Calculation of Williamson's SOI 
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
laea <- "+proj=laea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
lcea <- "+proj=cea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
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

# Import kernelUD rasters -------------------------------------------------
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/final_ud/")
path <- getwd()
file_list <- list.files(path=path)

for (i in 1:length(file_list)){
  load(file_list[i])
}

# converting to raster layer from estUD 

midLAAL <- raster(as(midLAAL,"SpatialPixelsDataFrame"))
midBFAL <- raster(as(midBFAL,"SpatialPixelsDataFrame"))
ternLAAL <- raster(as(ternLAAL,"SpatialPixelsDataFrame"))
ternBFAL <- raster(as(ternBFAL,"SpatialPixelsDataFrame"))

# Can skip to WSOI calculations after this - probability density is already normalized and 
# the rasters you loaded in are averaged already from your master_script

rastLM <- midLAAL
rastBM <- midBFAL
rastLT <- ternLAAL
rastBT <- ternBFAL


# Make density grid rasters -----------------------------------------------

# Averaging 

# Can easily average through some version of this: 
# ud1 <- stack(lapply(lm_kde, raster))
# plot(raster(lm_kde[[1]]))
lm_list <- split(lm, f=lm$id)
lt_list <- split(lt, f=lt$id)
bm_list <- split(bm, f=bm$id)
bt_list <- split(bt, f=bt$id)

# now you need to set extent and empty raster, then lapply(lm_avg,rasterize) to create a lot 
# of rasters, then stack those rasters, then avg

all_data_1<-all_data[,1:3]
# set up an 'empty' raster, here via an extent object derived from your data
e <- extent(all_data_1)
# In the next line, setting ext and res lets you not specify other defaults
r <- raster(e, ext=e, res=c(2.7,2.7), crs="+proj=longlat +datum=WGS84" ) # 2.7 degrees is 300km x 300km 
# you need to provide a function 'fun' for when there are multiple points per cell

# The key here is that you are building individual density grids. Each track gets a density grid (each cell is number of 
# relocations in that cell divided by total number of relocations in that track). Then you are averaging all 
# rasters in the stack to get an average value of each cell. 

lm_avg1 <- lapply(lm_list,"[", c('x','y')) # grabs just x and y columns for each element of list
lm_avg2 <- lapply(lm_avg1, function(x) rasterize(x,r,fun='count'))
for (i in 1:length(lm_avg2)){
  values(lm_avg2[[i]]) <- values(lm_avg2[[i]])/nrow(lm_avg1[[i]])
}
lm_stack <- raster::stack(lm_avg2)
plot(lm_stack)
lm_averaged <- raster::calc(lm_stack, mean, na.rm=T)
plot(lm_averaged)

bm_avg1 <- lapply(bm_list,"[", c('x','y')) # grabs just x and y columns for each element of list
bm_avg2 <- lapply(bm_avg1, function(x) rasterize(x,r,fun='count'))
for (i in 1:length(bm_avg2)){
  values(bm_avg2[[i]]) <- values(bm_avg2[[i]])/nrow(bm_avg1[[i]])
}
bm_stack <- raster::stack(bm_avg2)
plot(bm_stack)
bm_averaged <- raster::calc(bm_stack, mean, na.rm=T)
plot(bm_averaged)

lt_avg1 <- lapply(lt_list,"[", c('x','y')) # grabs just x and y columns for each element of list
lt_avg2 <- lapply(lt_avg1, function(x) rasterize(x,r,fun='count'))
for (i in 1:length(lt_avg2)){
  values(lt_avg2[[i]]) <- values(lt_avg2[[i]])/nrow(lt_avg1[[i]])
}
lt_stack <- raster::stack(lt_avg2)
plot(lt_stack)
lt_averaged <- raster::calc(lt_stack, mean, na.rm=T)
plot(lt_averaged)

bt_avg1 <- lapply(bt_list,"[", c('x','y')) # grabs just x and y columns for each element of list
bt_avg2 <- lapply(bt_avg1, function(x) rasterize(x,r,fun='count'))
for (i in 1:length(bt_avg2)){
  values(bt_avg2[[i]]) <- values(bt_avg2[[i]])/nrow(bt_avg1[[i]])
}
bt_stack <- raster::stack(bt_avg2)
plot(bt_stack)
bt_averaged <- raster::calc(bt_stack, mean, na.rm=T)
plot(bt_averaged)

rastLM <- lm_averaged
rastBM <- bm_averaged
rastLT <- lt_averaged
rastBT <- bt_averaged

# Old method (wrong, doesn't average, just builds density grids)

# rastLM <- rasterize(lm[,2:3], r, fun='count')
# rastBM <- rasterize(bm[,2:3], r, fun='count')
# rastLT <- rasterize(lt[,2:3], r, fun='count')
# rastBT <- rasterize(bt[,2:3], r, fun='count')
# plot(rastLM)
# plot(rastBM)
# plot(rastLT)
# plot(rastBT)
# # Check lengths to calculate density
# # x_vector <-as.numeric(na.omit(values(x)))
# # x_sum <-sum(x_vector)
# # y_vector <-as.numeric(na.omit(values(y)))
# # y_sum <-sum(y_vector)
# 
# values(rastLM) <- rastLM[]/nrow(lm) # 9681
# values(rastBM) <- rastBM[]/nrow(bm) # 4163
# values(rastLT) <- rastLT[]/nrow(lt) # 14799
# values(rastBT) <- rastBT[]/nrow(bt) # 10870
# 
# plot(rastLM)
# plot(rastBM)
# plot(rastLT)
# plot(rastBT)


# Normalization -----------------------------------------------------------

# Normalize each raster, then average those rasters in the next step. 
# Normalizing to account for sample size differences in data points. 
# Unless there are two identical max values in a UD, there will be only 1 
# entry that equals 1 after normalizing. 

# to add in averaging step, you need to make rastLM the averaged raster...
# create rasters for each individual bird, raster stack, raster average

### normalize lm ### 
print(length(values(rastLM)!=0)) # make sure its NA, not 0
lm_holder <- values(rastLM)
min_value <- min(lm_holder[!is.na(lm_holder)])
max_value <- max(lm_holder[!is.na(lm_holder)])
lm_holder <- sapply(lm_holder[!is.na(lm_holder)], function(x) (x-min_value)/(max_value-min_value))
length(rastLM[!is.na(values(rastLM))])
length(lm_holder)
rastLM[!is.na(values(rastLM))]<-lm_holder
plot(rastLM)

### normalize bm ### 
print(length(values(rastBM)!=0)) # how many entries don't equal 0? 
bm_holder <- values(rastBM)
min_value <- min(bm_holder[!is.na(bm_holder)])
max_value <- max(bm_holder[!is.na(bm_holder)])
bm_holder <- sapply(bm_holder[!is.na(bm_holder)], function(x) (x-min_value)/(max_value-min_value))
length(rastBM[!is.na(values(rastBM))])
length(bm_holder)
rastBM[!is.na(values(rastBM))]<-bm_holder
plot(rastBM)

### normalize lt ### 
print(length(values(rastLT)!=0)) # make sure its NA, not 0
lt_holder <- values(rastLT)
min_value <- min(lt_holder[!is.na(lt_holder)])
max_value <- max(lt_holder[!is.na(lt_holder)])
lt_holder <- sapply(lt_holder[!is.na(lt_holder)], function(x) (x-min_value)/(max_value-min_value))
length(rastLT[!is.na(values(rastLT))])
length(lt_holder)
rastLT[!is.na(values(rastLT))]<-lt_holder
plot(rastLT)

### normalize bt ### 
print(length(values(rastBT)!=0)) # how many entries don't equal 0? 
bt_holder <- values(rastBT)
min_value <- min(bt_holder[!is.na(bt_holder)])
max_value <- max(bt_holder[!is.na(bt_holder)])
bt_holder <- sapply(bt_holder[!is.na(bt_holder)], function(x) (x-min_value)/(max_value-min_value))
length(rastBT[!is.na(values(rastBT))])
length(bt_holder)
rastBT[!is.na(values(rastBT))]<-bt_holder
plot(rastBT)



# Calculate Williamson's SOI ----------------------------------------------
## EXAMPLE/TEST ##
# Example first! From Williamson 1993, 
# "Values greater than one represent
# greater overlap than would be expected with uniform
# prey and predator distributions, where the upper limit
# is determined by the number of patches sampled, which 
# is indeed what we see here: 760!
# LM x BM
oLM <- rastLM
oBM <- rastBM
m <- 760
num <- (oLM[358]*oBM[358])*m
denom <- (oLM[358]*oBM[358])
num/denom
###

# Ok, actually: 

# LM x BM
oLM <- rastLM
oBM <- rastBM
oLM[is.na(oLM)]<-0
oBM[is.na(oBM)]<-0
m <- 893 # NUMBER OF CELLS - THIS CHANGES BASED ON IF YOU ARE USING DENSITY RASTERS OR KERNEL UD (estUD imports)!
num <- sum((oLM[]*oBM[])*m)
denom <- (sum(oLM[]))*(sum(oBM[]))
LM_BM_SOI <- num/denom 
LM_BM_SOI # 4.557743
# density grid value 2.25
# non-averaged density grid value 2.64
# Julia got ~4.06

plot(oLM)
plot(oBM)

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

# BM x BT
oBM <- rastBM
oBT <- rastBT
oBM[is.na(oBM)]<-0
oBT[is.na(oBT)]<-0
m <- 893
num <- sum((oBM[]*oBT[])*m)
denom <- (sum(oBM[]))*(sum(oBT[]))
BM_BT_SOI <- num/denom # 2.908
# density grid value I get 2.189821
# non-averaged density grid value 1.714146 
BM_BT_SOI
plot(oBM)
plot(oBT)

# LM x LT
oLM <- rastLM
oLT <- rastLT
oLM[is.na(oLM)]<-0
oLT[is.na(oLT)]<-0
m <- 893
num <- sum((oLM[]*oLT[])*m)
denom <- (sum(oLM[]))*(sum(oLT[]))
LM_LT_SOI <- num/denom # 7.129
# density grid value I get 2.25901
# non-averaged density grid value 4.523776
LM_LT_SOI
plot(oLM)
plot(oLT)


# LT x BT
oLT <- rastLT
oBT <- rastBT
oLT[is.na(oLT)]<-0
oBT[is.na(oBT)]<-0
m <- 893
num <- sum((oLT[]*oBT[])*m)
denom <- (sum(oLT[]))*(sum(oBT[]))
LT_BT_SOI <- num/denom # 1.694
# density grid value I get 1.756759
# non-averaged density grid value 0.9858841
LT_BT_SOI





