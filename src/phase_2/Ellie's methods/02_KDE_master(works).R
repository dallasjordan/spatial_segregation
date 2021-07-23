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

YEAR = "2008"
# SET WORKING DIRECTORY FOR RELATIVE PATH (if not working in rproj)
setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",YEAR,"/"))
# RESET WORKING DIRECTORY TO PROJECT DIRECTORY
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/")
load('LAALdata_for_julia.Rdata')

# LOAD YOUR DATA from you project data directory (in this case the read.csv function is skipping the first 6 rows). Run for Midway and Tern folders
# Assumes you have all individuals in one csv.. if that is not the case use code on lines 27&28
# data <- read.csv("./data/GBBG.Individual.Tag.data/GBBG\ Tag\ 403.csv", skip = 6)
files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",YEAR), pattern = "csv")
data <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
#data <- dplyr::slice(data, 1:1230)

# SORT OUT YOUR `datetime` COLUMN and make sure it is in the accurate format
data$dtime <- as.POSIXct(strptime(data$dtime, format = "%m/%d/%y", tz = "UTC"))

# Convert to Pacific centric coords

# Converting +ve/-ve longitude to +ve longitude, ve longitude 
convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}

data <- data %>% dplyr::mutate(Longitude=convlon(Longitude))

# convert longitudes from 0 to 360 back to +/-, don't run this!
# convlon.back <- function(lon){
#   ifelse(lon > 180, lon - 360, lon)
# }

# Define some projections (YOU WILL HAVE TO THINK ABOUT WHAT PROJECTION MAKES SENSE FOR YOUR DATA)
# GREAT PLACE TO LOOKUP Proj4 strings - https://spatialreference.org/ref/
raw_data_proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # helpful to have WGS84 coordinate system as a variable for future calls.
#CHANGE THIS. preserve distance (minimize distance distortion)
target_proj <- "+proj=utm +zone=1 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # this is a mollweide projection centered on the mean(data$longitude) which was ~ 240.  This is generally a decent projection.

# Kernel Density estimation of tracking data and shapefile export for visualization in Arcmap
# as it is written, it generates KDEs for each individual tag-id, stacks those rasters, and then takes the mean of each raster cell, then generates 
# the contour levels from the averaged density raster

# KS package reverses the contour levels to 100-x, so the 95% is named "5%", etc.
contLevels <- c("5%","10%","50%") # these are the contour levels that will be generated in the kernel density estimation

# empty df to hold the contour percent associated w/ the polygon(s) and any other attribute table information - You may only need the "cntr_level"
attribute_df <- data.frame(id=character(), cntr_level=character())
# attribute_df <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("id", "cntr_level"))

# Empty list to hold the generated kde estimates for each bird and polygon shapefile outputs
fhats_ind <- list() # holds list of the KDE estimates
polyList <- list() # holds list of polygons

# (If you need to subset your data, do so here)
# E.G. Subset data frame for the post-breeding points 
prod_df <- na.omit(data)

# Fix an issue with duplicate timestamps - add 1 second to each timestamp if the previous is identical. 
for (i in 4:nrow(prod_df)){
    if (prod_df$dtime[i]==prod_df$dtime[i-1]){
      prod_df$dtime[i] <- prod_df$dtime[i]+seconds(1)
  }
    if (prod_df$dtime[i]==prod_df$dtime[i-2]){
      prod_df$dtime[i] <- prod_df$dtime[i]+seconds(2)
  }
    if (prod_df$dtime[i]==prod_df$dtime[i-3]){
      prod_df$dtime[i] <- prod_df$dtime[i]+seconds(3)
  }
    if (prod_df$dtime[i]==prod_df$dtime[i-3]){
      prod_df$dtime[i] <- prod_df$dtime[i]+seconds(4)
  }
}
  

# check by plotting x and y
plot(prod_df$Longitude, prod_df$Latitude) 

# Set target resolution and extent buffer for the subset of data (you should not need to touch this)
resolution <- 55000 # (target resolution (roughly one half degree) in meters and is based off 1 degree Latitude at equater = ~110 km)
#This number must be smaller than the error of my GLS data (~300km)
extent_buffer <- 2.2e6 # in meters 

# Project our prod_df to desired target projection (in this case UTM Zones 54 to 10) 
# Kernel density analysis in package "ks" requires that the coordinates are in measurement units == meters). 

xy <- data.frame(prod_df$Longitude, prod_df$Latitude)
# Remove any NA values (return to this later!)
xy <- na.omit(xy)
spdf <- SpatialPointsDataFrame(coords=xy, data=prod_df, proj4string = CRS(raw_data_proj))

spdf_UTM1 <- spTransform(spdf, CRS(target_proj)) #transform projection to UTM zone 1

prod_df <- as.data.frame(spdf_UTM1)
prod_df <- prod_df %>% mutate(x = prod_df.Longitude, y = prod_df.Latitude)  %>% dplyr::select(-c(prod_df.Longitude, prod_df.Latitude))


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

index <- 1 # index to track each polygon - possible to have several polygons per contour

# empty the lists that store fhat kde objects and fhat$estimate matrices

fhats_ind <- list() 
denestimates_ind <- list()

for (i in 1:length(prod_dfL)) {
  
  #subset each invidual dataframe from which to draw individuals for kde analysis
  curdf_ind <- prod_dfL[[1]]
  ll= matrix(c(curdf_ind$x,curdf_ind$y), ncol=2) #ll represents x, y (corresponding to long/lat)
  
  # Get MISE Bandwidth - this is a very specific bandwidth that we chose to use with MiCO... you may choose a different bandwidth 
  curdf_move <- move(x = curdf_ind$x, y=curdf_ind$y, as.POSIXct(curdf_ind$dtime, format="%Y-%m-%d %H:%M:%S", tz="UTC"), data = curdf_ind, proj = CRS(target_proj), animal=curdf_ind$id)
  curdf_telemetry <- as.telemetry(curdf_move, projection = CRS(target_proj))
  M.IID <- ctmm.fit(curdf_telemetry)
  UD0 <- akde(curdf_telemetry,M.IID)
  H_MISE = UD0$H
  # Should explore a bandwidth that's appropriate for me! Check ellie's email for Bandwidth reference, read more about it. 
  
  # kernel density analysis from package "ks"
  fhat <- kde(x=ll, H = H_MISE,  gridsize = gridsz, xmin = x_min, xmax = x_max, approx.cont = TRUE)
  title(main = paste0(unique(curdf_ind$id), " KDE"))
  # Plot output to check the extent buffer in real time and make sure that kde predict is predicting within bounds of the extent set in the kde() call.  
  plot(fhat, cont = c(5, 10, 25, 50, 75, 90, 95))
  #title(main = paste(unique(curdf_ind$id), unique(curdf_ind$segment)))
  
  # shove this individual's estimate into our list of individual estimates  
  fhats_ind[[i]] <-  fhat 
  
  
  # shove individual density estimates into a list of matrixes for future matrix calculates (mean)
  denestimates_ind[[i]] <- fhat$estimate 
  
  
} #end individual for loop (l) 

# UD -> pretty much the same as a bivariate kernel density estimate. more digestible/ecological interpretation

# Perform some matrix calculations on the list of individual fhat$estimates.  This will generate an averaged matrix of density estimates on the grid size and extent that we specified above.
mean_ind_est <- apply(simplify2array(denestimates_ind), 1:2, mean)

# CHEAT METHOD::::Replace fhat$estimate with averaged month matrix and extract contours using ks::contourLevels() which calculates the highest density region
i_fhat <- fhats_ind[[1]] # take the first fhat in the stack of individual as a "shell" fhat kde object
i_fhat$estimate <- mean_ind_est # replace the $estimate matrix with the averaged mean_month_est
i_fhat$x <- matrix(c(prod_df$x,prod_df$y), ncol=2) #replace the coords with the coords for the entire subset
i_fhat$cont <- contourLevels(i_fhat, cont = 1:99, approx = TRUE) #calculate the contour levels using new averaged kde object

plot(i_fhat,cont = c(50, 90, 95)) # plot to double check
title("whatever you want")


# Define some variables for which to name the polygon output of this KDE
ID <- paste0("KDE_midway_LAAL","_",YEAR)  # or something



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
rgdal::writeOGR(spPolyDF, shapefilewritepath, layer= ID ,driver="ESRI Shapefile", verbose=TRUE, overwrite=TRUE)
