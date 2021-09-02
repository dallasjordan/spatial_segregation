# Script to do permutations for "side project" with Julia and Lesley, showing how scale impacts overlap calculations and 
# that previous methods don't necessarily work for our data set
# This script handles WSOI
# Dallas Jordan
# August 12

# This script 1. Loads in data,
#             2. Loads kernelUD (estUD) objects,
#             3. Convert to raster
#             4. 

# Setup -------------------------------------------------------------------

library(adehabitatHR)
library(dplyr)

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
library(raster)

# Load in tracks for species A and species B
# Can skip this if you want, just load in the already averaged and normalized estUD files in the next sections

sppA_read <- file.choose()
sppA <- readRDS(sppA_read)

sppB_read <- file.choose()
sppB <- readRDS(sppB_read)

all_data <- rbind(sppA, sppB)

all_data_1<-all_data[,1:3]
sp::coordinates(all_data_1) <- c("lon", "lat")

# Calculate all estUD for each animalID: href bandwidth, arbitrary grid=100 -> nicer images but longer calculation time
results <- kernelUD(all_data_1, grid=100,same4all=T,extent=0.1)
image(results[[100]])

# seperate by species
names <- names(results)
a_kde <- results[grep("A",names)]
b_kde <- results[grep("B",names)]


# Begin calculations ------------------------------------------------------

# for Mac
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/")
path <- getwd()
file_list <- list.files(path=path)

rastA <- raster(as(sppA_UD,"SpatialPixelsDataFrame"))
rastB <- raster(as(sppB_UD,"SpatialPixelsDataFrame"))

# Can skip to WSOI calculations - probability density is already normalized and 
# the rasters you loaded in are averaged already from your master_script

# Calculate Williamson's SOI ----------------------------------------------

# AL x AB
oA <- rastA
oB <- rastB
oA[is.na(oA)]<-0
oB[is.na(oB)]<-0
m <- 6700 # NUMBER OF CELLS - THIS CHANGES BASED ON IF YOU ARE USING DENSITY RASTERS OR KERNEL UD (estUD imports)!
num <- sum((oA[]*oB[])*m)
denom <- (sum(oA[]))*(sum(oB[]))
A_B_SOI <- num/denom # I get 0.25 when using estUD rasters
A_B_SOI

plot(oA)
plot(oB)

# Randomization tests -----------------------------------------------------

# Randomize, then apply calculations
iA <- oA
iB <- oB

pool <- 1:13400
vector_values_A <- iA@data@values
vector_values_B <- iB@data@values
vector_values_all <- append(vector_values_A, vector_values_B)

a_nums <- sample(pool, 6700)
b_nums <- setdiff(pool,a_nums)

iA <- 




