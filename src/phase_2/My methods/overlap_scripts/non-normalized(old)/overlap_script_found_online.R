#Customised R function to test for spatial segregation
## RandomOverlap ##########################################################
#Function built by Sonia Sánchez on September 2017 – Contact: soniasg9@gmail.com
#Function to investigate differences in area usage between two groups (i.e. spatial segregation). Transferable to tracking data sets with UTM coordinates and already filtered locations to be used in Kernel analysis.
#Copy and paste in your R script to use it

#ARGUMENTS
##### DataGroup = dataframe, datatable or SpatialPointsDataFrame with the tracking data, 
################# with x (UTM coordinates), y (UTM coordinates) and ID (individual trips) as var_seg (grouping variable)
##### Scale = smoothing factor (h) used in the kernel analysis, should be provided in km (default value 1 km)
##### Iteration = number of times DataGroup is iterated to obtained random overlaps (default value 100)
##### UDLev = quantile (percent parameter) to be used in  the Utilisation Distribution (default value 50)
##### res = grid parameter for the Utilisation Distribution
##### method = overlap index method to be used (default value "UDOI" method)

#WHAT DOES THE FUNCTION DO?
# This function iteratively subsamples a dataset (DataGroup) of tracking data and calculates the overlap index between two groups in order to investigate if there are differences in area usage between the two groups (i.e. spatial segregation).
# At each iteration, random trips (ID) are assigned to one of the two levels of the grouping variable (var_seg), and the rest of individuals are assigned to the other level. Random samples keep the ratio of each group level.
# The overlap index is calculated for the observed data and for all the random samples. 

#WHAT DOES THE FUNCTION RETURN?
# (1) PropLower = which is the proportion of random overlaps lower than the observed overlap. If PropLower < 0.05 (or different value decided by user), then we refuse the NULL hypothesis of no differences in the area use between sites.
# (2) Two column data frame, 1st column is the iteration number and 2nd column has the overlap index, with the first value being the observed overlap, followed by all the overlaps obtained by randomisation. 

library(tidyverse)
library(data.table)
library(sf)
library(sp)
library(rgdal)

RandomOverlap <- function(DataGroup, Scale=1, Iteration=100, UDLev = 50, res = 400, method = "UDOI")
{
  require(sp)
  require(geosphere)
  require(rgdal)
  require(adehabitatHR)
  
  if(!"x" %in% names(DataGroup)) stop("x field does not exist")
  if(!"y" %in% names(DataGroup)) stop("y field does not exist")
  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")
  if(!"var_seg" %in% names(DataGroup)) stop("var_seg field does not exist")
  if((class(DataGroup)!= "SpatialPointsDataFrame")[1]) # Convert to SpatialPointsDataFrame and project if DataGroup is not a SpatialPointsDataFrame object
  {
    DataGroup <- SpatialPointsDataFrame(coordinates(cbind(DataGroup$x,    DataGroup$y)), data = DataGroup, proj4string = CRS("+proj=utm +zone=2 +datum=WGS84"))
  } else {DgProj <- DataGroup@proj4string}
  
  variable <- unique(DataGroup$var_seg) # Get the levels of the factor we want to use to group the trips
  UIDs <- unique(DataGroup$ID) # Get the all the unique trip IDs
  nindv1 <- length(unique(DataGroup$ID[ DataGroup$var_seg == variable[1]])) # Get the number of trips for the first level of var-seg
  nindv2 <- length(unique(DataGroup$ID[ DataGroup$var_seg == variable[2]])) # Get the number of trips for the second level of var-seg
  Ntrips <- length(UIDs) # Get the total number of trips
  Output <- data.frame(NIter = 1:(Iteration+1), overlap = rep(0, (Iteration+1))) # Prepare output data frame that will content observed and random overlaps
  KUD <- adehabitatHR::kernelUD(DataGroup[ , "var_seg"], h = Scale*1000, grid = res, same4all = T) # Calculate observed Kernel UD 
  Overobs <- adehabitatHR::kerneloverlaphr( KUD, meth = method, percent = UDLev, conditional = T) # Calculate observed overlap (overlap index indicated by "method") at contour indicated by percent = UDLev
  Overobs <- Overobs[ 1, 2 ] # Get the observed overlap value
  Output[ 1, 2 ] <- Overobs # Put the observed overlap in the first row of the second column of the output data frame
  for (i in 1:Iteration) # Start loop to calculate random overlaps
  {
    RanID1 <- sample(UIDs, nindv1, replace=F) # Get random trips IDs, as many as in level 1 of var_seg
    RanID2 <- UIDs[!( UIDs %in% RanID1 )] # The rest of the trip IDs that aren't in sample RanID1 go to RandID2
    Coord1 <- coordinates(DataGroup[DataGroup$ID %in% RanID1,]) %>% 
      as.data.frame(.) %>% setDT() %>% setnames(c("x", "y")) %>% .[ , var_seg := variable[1]] # Get the coordinates from the trip IDs in RanID1, convert to data frame and incorporate the level 1 of ver_seg
    Coord2 <- coordinates(DataGroup[DataGroup$ID %in% RanID2,]) %>% 
      as.data.frame(.) %>% setDT() %>% setnames(c("x", "y")) %>% .[ , var_seg := variable[2]] # Same than before for the rest of trip IDs
    RanCoord <- rbind(Coord1, Coord2) # Get the random sample together in a dataframe
    RanCoord <- SpatialPointsDataFrame(coordinates(cbind(RanCoord$x, RanCoord$y)), data = RanCoord, proj4string = CRS("+proj=utm +south +zone=55 +datum=WGS84")) # Convert random sample to SpatialPointsDataFrame to calculate kernel UD
    KDE.Surface <- adehabitatHR::kernelUD(RanCoord[ , "var_seg"], h= Scale*1000, grid=res, same4all=T) # Kernel UD of the random sample
    Overran <- adehabitatHR::kerneloverlap( KDE.Surface, meth = method, percent = UDLev, conditional = T) # Overlap of the random sample
    Overran <- Overran[ 1, 2 ] # Get random overlap value
    Output[ i+1, 2] <- Overran # Put the random overlap value in the output data frame
  } # End of loop. It will run as many times as indicated by the value of Iteration argument
  PropLower <- length(Output$overlap[ Output$overlap < Output[1,2]]) /   (length(Output$NIter)-1) # Calculate proportion of random overlaps that are smaller than observed overlap. This will be the empirical P-value.
  print(paste("The proportion of random", method, UDLev, "% UD overlaps lower than observed overlap is PropLower =",  PropLower, sep = " "))
  return(Output) # Function will return the output data frame with the observed overlap and all the random overlaps
}

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_midway.Rdata")
BFALmid <- BFAL
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_tern.Rdata")
BFALtern <- BFAL

BFALmid$ID <- 1
BFALmid$var_seg <- "BFALmid"
BFALtern$ID <- 1
BFALtern$var_seg <- "BFALtern"

both <- rbind(BFALmid,BFALtern)
both <- both[,2:5]


test<-RandomOverlap(both,Scale=1,Iteration=10,UDLev=95,res=300, method="UDOI")



