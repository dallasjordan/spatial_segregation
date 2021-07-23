# Author: Dallas Jordan
# Modified from script written by Eleanor Heywood
# Date last modified: 09-02-2020
# Abstract:
# Using the R-package "ks", this script assumes you've already done some data cleaning and quality control.  It ingests location data
# from tracked animals (e.g. albatross, gulls, etc) and produces a kernel density estimate for each individual. The raster surface outputs
# for each individual are then stacked and averaged, the 95, 90, 75, 50 and 25 contours extracted from this averaged raster surface.


# Load necessary libraries (you will have to install if you have never installed before)
library(ks) #statisticians recommend this as the best kernel density estimator package, you'll be using the function ks::kde()
?ks::kde
library(move)
library(tidyverse)
library(sp)
library(rgdal)
library(raster)
library(ctmm)
library(lubridate)


setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures")
load('BFALdata_tern.Rdata')
data_tern<-BFAL

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures")
load('BFALdata_midway.Rdata')
data_midway<-BFAL

data <- rbind(data_tern,data_midway)

# Define some projections (YOU WILL HAVE TO THINK ABOUT WHAT PROJECTION MAKES SENSE FOR YOUR DATA)
# GREAT PLACE TO LOOKUP Proj4 strings - https://spatialreference.org/ref/
raw_data_proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # helpful to have WGS84 coordinate system as a variable for future calls.
#CHANGE THIS. preserve distance (minimize distance distortion)
target_proj <- "+proj=laea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
# you've used this before: 
# target_proj <- "+proj=utm +zone=1 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # this is a mollweide projection centered on the mean(data$longitude) which was ~ 240.  This is generally a decent projection.

# Kernel Density estimation of tracking data and shapefile export for visualization in Arcmap
# as it is written, it generates KDEs for each individual tag-id, stacks those rasters, and then takes the mean of each raster cell, then generates 
# the contour levels from the averaged density raster

# KS package reverses the contour levels to 100-x, so the 95% is named "5%", etc.
contLevels <- c("5%","50%", "90%") # these are the contour levels that will be generated in the kernel density estimation

# empty df to hold the contour percent associated w/ the polygon(s) and any other attribute table information - You may only need the "cntr_level"
attribute_df <- data.frame(id=character(), cntr_level=character())
# attribute_df <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("id", "cntr_level"))

# Empty list to hold the generated kde estimates for each bird and polygon shapefile outputs
fhats_ind <- list() # holds list of the KDE estimates
polyList <- list() # holds list of polygons

# (If you need to subset your data, do so here)
# E.G. Subset data frame for the post-breeding points 
prod_df <- na.omit(data)

# check by plotting x and y
plot(prod_df$x, prod_df$y) 

# Set target resolution and extent buffer for the subset of data (you should not need to touch this)
resolution <- 300000 # (target resolution (roughly one half degree) in meters and is based off 1 degree Latitude at equater = ~110 km)
#This number must be smaller than the error of my GLS data (~300km)
extent_buffer <- 2.2e6 # in meters 

# Project our prod_df to desired target projection (in this case UTM Zones 54 to 10) 
# Kernel density analysis in package "ks" requires that the coordinates are in measurement units == meteres). 

xy <- data.frame(prod_df$x, prod_df$y)
# Remove any NA values (return to this later!)
xy <- na.omit(xy)
spdf <- SpatialPointsDataFrame(coords=xy, data=prod_df, proj4string = CRS(raw_data_proj))

spdf_UTM1 <- spTransform(spdf, CRS(target_proj)) #transform projection to UTM zone 1

prod_df <- as.data.frame(spdf_UTM1)
prod_df <- prod_df %>% mutate(x = prod_df.x, y = prod_df.y)  %>% dplyr::select(-c(prod_df.x, prod_df.y))


# Set extents and gridsize for all rasters (you should not need to touch this - it assumes the bounds of your data)
x_extent <- c(min(prod_df$x-extent_buffer), max(prod_df$x)+extent_buffer)
x_range <- (max(prod_df$x)+extent_buffer)-(min(prod_df$x)-extent_buffer)

y_extent <- c(min(prod_df$y)-extent_buffer, max(prod_df$y)+extent_buffer)
y_range <- (max(prod_df$y)+extent_buffer) - (min(prod_df$y)-extent_buffer)

gridsz <- c(x_range/resolution, y_range/resolution)
x_min <- c(x_extent[1], y_extent[1])
x_max <- c(x_extent[2], y_extent[2])

# Loop over each individual, generating individual kdes, averaging them, and extracting the contour polys
prod_dfL <- split(prod_df, prod_df$id) #split dataframe by tag-id in order to loop over each individual
# make all IDs the same 'individual'
prod_df$id <- "BFAL"

index <- 1 # index to track each polygon - possible to have several polygons per contour

# empty the lists that store fhat kde objects and fhat$estimate matrices

fhats_ind <- list() 
denestimates_ind <- list()

# normal loop here, but since I'm just doing all LAAL, no loop:

#subset each invidual dataframe from which to draw individuals for kde analysis
curdf_ind <- prod_df
ll= matrix(c(curdf_ind$x,curdf_ind$y), ncol=2) #ll represents x, y (corresponding to long/lat)

# kernel density analysis from package "ks"
fhat <- kde(x=ll,  gridsize = gridsz, xmin = x_min, xmax = x_max, approx.cont = TRUE) # just using the default bandwidth
title(main = paste0(unique(curdf_ind$id), " KDE"))
# Plot output to check the extent buffer in real time and make sure that kde predict is predicting within bounds of the extent set in the kde() call.  
plot(fhat, cont = c(5, 10, 25, 50, 75, 90, 95))
#title(main = paste(unique(curdf_ind$id), unique(curdf_ind$segment)))

# shove this individual's estimate into our list of individual estimates  
fhats_ind <- fhat 


# shove individual density estimates into a list of matrixes for future matrix calculates (mean)
denestimates_ind <- fhat$estimate 



# UD -> pretty much the same as a bivariate kernel density estimate. more digestible/ecological interpretation

# Perform some matrix calculations on the list of individual fhat$estimates.  This will generate an averaged matrix of density estimates on the grid size and extent that we specified above.
mean_ind_est <- apply(simplify2array(denestimates_ind), 1:2, mean)

# CHEAT METHOD::::Replace fhat$estimate with averaged month matrix and extract contours using ks::contourLevels() which calculates the highest density region
i_fhat <- fhats_ind # take the first fhat in the stack of individual as a "shell" fhat kde object
i_fhat$estimate <- mean_ind_est # replace the $estimate matrix with the averaged mean_month_est
i_fhat$x <- matrix(c(prod_df$x,prod_df$y), ncol=2) #replace the coords with the coords for the entire subset
i_fhat$cont <- contourLevels(i_fhat, cont = 1:99, approx = TRUE) #calculate the contour levels using new averaged kde object

plot(i_fhat,cont = c(50, 75, 90, 95)) # plot to double check
title("whatever you want")


# Define some variables for which to name the polygon output of this KDE
SPECIES <- "BFAL"
ID <- paste0("KDE_combined_",SPECIES)  # or something



# Generate the contour levels for export of contour polygons into ESRI shapefile
for (j in 1:length(contLevels)) {
  cl <- contourLines(x=i_fhat$eval.points[[1]],y=i_fhat$eval.points[[2]],
                     z=i_fhat$estimate,levels=i_fhat$cont[contLevels[j]])
  # inner loop over each contour line for the level
  for (i in 1:length(cl)) {
    poly <- data.frame(cl[[i]]$x, cl[[i]]$y)
    poly <- rbind(poly,poly[1,]) #close the polygon
    polyList[index] <- list(Polygons(list(Polygon(poly)),ID=(index)))
    
    # make sure the "to" and "from" columns accurately represent the connectivity of the corridor feature.  If not a corridor, then these columns should be populated with "NA"
    attribute_df <- bind_rows(attribute_df, c(id = ID, cntr_level =  100 - as.numeric(sub("%", "", contLevels[j]))))
    index <- index + 1
  } # end of inner loop for each contour line level
} # end of contour level loop


# make the spatial polygons, and join the dataframe holding the important bits of info
spPoly <- SpatialPolygons(polyList, proj4string= CRS(target_proj))

spPolyDF <- SpatialPolygonsDataFrame(spPoly,attribute_df)
raster::plot(spPolyDF) # can plot polygon to check and make sure it makes sense
shapefilewritepath <- getwd() # points to the folder location to write the shapefile to (DOES NOT INCLUDE the shapefile layer name)
# IFF the species crosses the dateline in the pacific, write out the shapefile in the target projection centered on the mean of the longitude (target_proj).  Once the shapefiles are in ArcMap, they can be reprojected to WGS84 for the MiCO online system.  IFF the species is not pacific centered, the data can be written out in WGS84.  to transform the shapefile from the target_proj to WGS84 use the following code:
# spPolyDF <- spTransform(spPolyDF, CRSobj = CRS(raw_data_proj))
save(fhat, file="BFAL_combined_KDE.Rdata")
rgdal::writeOGR(spPolyDF, shapefilewritepath, layer= ID ,driver="ESRI Shapefile", verbose=TRUE, overwrite=TRUE)
