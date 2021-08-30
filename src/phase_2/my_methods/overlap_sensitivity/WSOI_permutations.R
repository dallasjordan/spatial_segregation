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

# Load in tracks for species A and species B

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
