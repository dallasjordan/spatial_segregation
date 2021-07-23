#can play around with this script (pacific_map is master for this)

library(sp)
library(maps)
library(maptools)   ## map2SpatialLines(), pruneMap()
library(rgdal)      ## CRS(), spTransform()
library(rgeos)      ## readWKT(), gIntersects(), gBuffer(), gDifference() 

## Convert a "maps" map to a "SpatialLines" map
makeSLmap <- function() {
  llCRS <- CRS("+proj=longlat +ellps=WGS84")
  wrld <- map("world", interior = FALSE, plot=FALSE,
              xlim = c(-179, 179), ylim = c(-89, 89))
  wrld_p <- pruneMap(wrld, xlim = c(-179, 179))
  map2SpatialLines(wrld_p, proj4string = llCRS)
}

## Clip SpatialLines neatly along the antipodal meridian
sliceAtAntipodes <- function(SLmap, lon_0) {
  ## Preliminaries
  long_180 <- (lon_0 %% 360) - 180
  llCRS  <- CRS("+proj=longlat +ellps=WGS84")  ## CRS of 'maps' objects
  eqcCRS <- CRS("+proj=eqc")
  ## Reproject the map into Equidistant Cylindrical/Plate Caree projection 
  SLmap <- spTransform(SLmap, eqcCRS)
  ## Make a narrow SpatialPolygon along the meridian opposite lon_0
  L  <- Lines(Line(cbind(long_180, c(-89, 89))), ID="cutter")
  SL <- SpatialLines(list(L), proj4string = llCRS)
  SP <- gBuffer(spTransform(SL, eqcCRS), 10, byid = TRUE)
  ## Use it to clip any SpatialLines segments that it crosses
  ii <- which(gIntersects(SLmap, SP, byid=TRUE))
  # Replace offending lines with split versions
  # (but skip when there are no intersections (as, e.g., when lon_0 = 0))
  if(length(ii)) { 
    SPii <- gDifference(SLmap[ii], SP, byid=TRUE)
    SLmap <- rbind(SLmap[-ii], SPii)  
  }
  return(SLmap)
}

## re-center, and clean up remaining streaks
recenterAndClean <- function(SLmap, lon_0) {
  llCRS <- CRS("+proj=longlat +ellps=WGS84")  ## map package's CRS
  newCRS <- CRS(paste("+proj=eqc +lon_0=", lon_0, sep=""))
  ## Recenter 
  SLmap <- spTransform(SLmap, newCRS)
  ## identify remaining 'scratch-lines' by searching for lines that
  ## cross 2 of 3 lines of longitude, spaced 120 degrees apart
  v1 <-spTransform(readWKT("LINESTRING(-62 -89, -62 89)", p4s=llCRS), newCRS)
  v2 <-spTransform(readWKT("LINESTRING(58 -89, 58 89)",   p4s=llCRS), newCRS)
  v3 <-spTransform(readWKT("LINESTRING(178 -89, 178 89)", p4s=llCRS), newCRS)
  ii <- which((gIntersects(v1, SLmap, byid=TRUE) +
                 gIntersects(v2, SLmap, byid=TRUE) +
                 gIntersects(v3, SLmap, byid=TRUE)) >= 2)
  SLmap[-ii]
}

## Put it all together:
Recenter <- function(lon_0 = -100, grid=FALSE, ...) {                        
  SLmap <- makeSLmap()
  SLmap2 <- sliceAtAntipodes(SLmap, lon_0)
  recenterAndClean(SLmap2, lon_0)
}

test <- Recenter(180)

