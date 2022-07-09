# Updated script to create Figure 5, which is my 4 panel figure that has 95th and 50th contours
# by island/species. This draws from my July 24 2021 scripts/individual KDE mapping, under "Individual_LAALBFAL_mapping", 
# they are titled e.g. midLAAL_contours_mapping
# Dallas Jordan
# Dec 18 2021
# Last updated March 7 2022

# Setup -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sf)

# CRITICAL: 
sf::sf_use_s2(FALSE)
# The above line is necessary as of April 17 2021, when 'sf' switched to using S2 spheroid 
# geometry instead of a GEOS routine that assumed projected coordinates. Until 'mapping' 
# updates its back-end, switching S2 geometry off is necessary for my script to work. 
# See: https://github.com/r-spatial/sf/issues/1649
# Great reference here: https://r-spatial.github.io/sf/articles/sf7.html

lcea <- "+proj=cea +lat_0=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
eqc <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
midway <- data.frame(longitude = -177.3761,
                     latitude = 28.2101)
tern <- data.frame(longitude = -166.284,
                   latitude = 23.870)

midway <- st_as_sf(midway, coords = c("longitude", "latitude"), 
                   crs = 4326, agr = "constant")
tern <- st_as_sf(tern, coords = c("longitude", "latitude"), 
                 crs = 4326, agr = "constant")


# Reference for this is EPSG 6933, WGS 84 / NSIDC EASE-Grid 2.0 Global
#    cea <- "+proj=cea +lat_0=20 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs" # lat_ts (Standard Parallel) describes latitudes where there is no distortion. 
#    _ts means "true scale", the lat or lon where there is no distortion. lon_0 = the natural origin for longitude

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/basemap/")
load("npac_base_res.Rdata")

# Import contours ---------------------------------------------------------

#load contours, as of Jun 22 they are SpatialPolygonsDataFrame
# for mac
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/individual/midLAAL/master_script_contours/")
load("vert95_midLAAL.Rdata")
load("vert50_midLAAL.Rdata")

# Saving one example of old method for future reference:
# ml95c <- st_as_sf(vert95_midLAAL) # 95th UD Contour
# ml95c <- ml95c %>%
#   st_transform(crs = 4326) %>% # transform to WGS84.
#   st_wrap_dateline() %>% # wrap around the dateline
#   st_shift_longitude() %>%
#   st_union(by_feature = TRUE) %>%
#   st_transform(crs = lcea)

ml95c <- st_as_sf(vert95_midLAAL) # 95th UD Contour
ml95c <- ml95c %>%
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

ml50c <- st_as_sf(vert50_midLAAL) # 50th UD Contour
ml50c <- ml50c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)


# to get color scale working: 
ml95c$id <- "foraging_range"
ml95c$cntr_level <- "95"

ml50c$id <- "focal_range"
ml50c$cntr_level <- "50"


# Generate basemap --------------------------------------------------------

# Ptolemy package, formerly "nPacMaps". Very amazing for North Pacific Mapping
# loading in data object from overallLAALBFAL_mapping to get the same proportions as my other figure for this manuscript

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/allLAAL_allBFAL/master_script_rasters/")
load("allBFAL_ud_vol_rast.Rdata")
allBFAL_rast <- allBFAL.ud.vol.raster
allBFAL.sp <- as(allBFAL_rast, "SpatialPixelsDataFrame")
gridded(allBFAL.sp)
spplot(allBFAL.sp, main="raster to sp - SpatialPixelsDataFrame")
allBFAL.rast.stars <- st_as_stars(allBFAL.sp, att=1)
allBFAL.rast.stars
allBFAL.rast.sf <- st_as_sf(allBFAL.rast.stars)
plot(allBFAL.rast.sf)
allBFAL.rast.sf <- allBFAL.rast.sf %>% 
  st_transform(crs = 4326) %>% 
  st_wrap_dateline() %>%
  st_shift_longitude() %>% 
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = lcea)
plot(allBFAL.rast.sf, main="PDCmerc with lat_ts=0") # raster is now sf object in PDC mercator
ab <- allBFAL.rast.sf

npac_base_i <- ptolemy::extract_gshhg(allBFAL.rast.sf, resolution = "i", epsg = NULL, buffer = 5000,
                                      simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS


# Plot  -------------------------------------------------------------------

figure <- ggplot() + 
  geom_sf(data=ml95c, color="black", fill="orangered", size=1.25) +
  geom_sf(data=ml50c, color="black", fill="orangered", size=1.25) +
  geom_sf(data=npac_base_res, fill="grey60") +
  scale_y_continuous(breaks = c(10, 30, 50, 70))+
  scale_x_continuous(breaks = c(-120,-140,-160,180,160,140,120))+
  geom_sf(data=midway, size=3,shape=17,fill="black",color="black")+
  geom_sf(data=tern, size=4,shape=18,fill="black",color="black")+
  coord_sf(expand=F)+
  theme_bw()
  # guides(fill=guide_legend(title="Contour level"))+
  #coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 12000000)) +
  # ggtitle("a")+
  # theme(plot.title = element_text(size=12))+
  # theme(legend.position = "none")+
  # theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
  #       legend.key.height = unit(0.5, 'cm'), #change legend key height
  #       legend.key.width = unit(0.5, 'cm'), #change legend key width
  #       legend.title = element_text(size=10), #change legend title font size
  #       legend.text = element_text(size=8))+ #change legend text font size
  # scalebar(x.min=6547510,x.max=953690,y.min=1803885,y.max=6141532, location = "bottomright", dist = 1500,
  #          dist_unit = "km", transform = FALSE, height = 0.03, st.size=3, st.dist = 0.05)+
  # annotation_north_arrow(location = "br", which_north = "true",
  #                        pad_x = unit(0.7, "in"), pad_y = unit(0.34, "in"),
  #                        style = north_arrow_fancy_orienteering,
  #                        height = unit(.75,"cm"),
#                        width = unit(.75,"cm"))+
# theme(axis.title.x=element_blank(),
#       axis.title.y=element_blank())

figure

