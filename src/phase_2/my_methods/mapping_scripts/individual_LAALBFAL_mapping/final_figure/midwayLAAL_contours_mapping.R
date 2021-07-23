# MIDWAY LAAL 
# Plotting contours from adehabitathr - used for publication figures
# June 30 2021
# Dallas Jordan

# adapting overallLAALBFAL_mapping script to plot my rasters and contours for Midway LAAL
# THERE IS ONE DEPENDENCY IN THIS SCRIPT - MUST LOAD IN THE npac_base Rdata CREATED IN overallLAALBFAL_mapping
# SCRIPT - THIS IS TO KEEP PLOTS LOOKING THE SAME! This object is loaded in right before plotting. 

library(mapproj)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgdal)
library(ggspatial)
library(ggmap)
library(ggthemes)
library(stars)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(dplyr)
library(ggnewscale)
library(ggsn)

# CRITICAL: 
sf::sf_use_s2(FALSE)
# The above line is necessary as of April 17 2021, when 'sf' switched to using S2 spheroid 
# geometry instead of a GEOS routine that assumed projected coordinates. Until 'mapping' 
# updates its back-end, switching S2 geometry off is necessary for my script to work. 
# See: https://github.com/r-spatial/sf/issues/1649
# Great reference here: https://r-spatial.github.io/sf/articles/sf7.html

#############################
######### SETUP #############
#############################

lcea <- "+proj=cea +lat_0=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
# Reference for this is EPSG 6933, WGS 84 / NSIDC EASE-Grid 2.0 Global
#    cea <- "+proj=cea +lat_0=20 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs" # lat_ts (Standard Parallel) describes latitudes where there is no distortion. 
#    _ts means "true scale", the lat or lon where there is no distortion. lon_0 = the natural origin for longitude

#############################
### IMPORT CONTOUR LINES  ###
#############################

#load contours, as of Jun 22 they are SpatialPolygonsDataFrame
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/midLAAL/master_script_contours/")
load("vert95_midLAAL.Rdata")
load("vert50_midLAAL.Rdata")
load("vert10_midLAAL.Rdata")

# Contours in 'sf'
ml95c <- st_as_sf(vert95_midLAAL) # 95th UD Contour
ml95c <- ml95c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)

ml50c <- st_as_sf(vert50_midLAAL) # 50th UD Contour
ml50c <- ml50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)


ml10c <- st_as_sf(vert10_midLAAL) # 10th UD Contour
ml10c <- ml10c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)

# additions July 7 to get viridis scale working
ml95c$id <- "foraging_range"
ml95c$cntr_level <- "95"

ml50c$id <- "focal_range"
ml50c$cntr_level <- "50"

ml10c$id <- "core_range"
ml10c$cntr_level <- "10"


###############################
### IMPORT AND PREP RASTERS ###
###############################

# load rasters
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/midLAAL/master_script_rasters/")
load("midLAAL_ud_vol_rast.Rdata")
midLAAL_rast <- midLAAL.ud.vol.raster
image(midLAAL_rast)
contour(midLAAL_rast,levels=c(95,50,10),add=T)

### LAAL ###

midLAAL.sp <- as(midLAAL_rast, "SpatialPixelsDataFrame")
gridded(midLAAL.sp)
spplot(midLAAL.sp, main="raster to sp - SpatialPixelsDataFrame")
midLAAL.rast.stars <- st_as_stars(midLAAL.sp, att=1)
midLAAL.rast.stars
midLAAL.rast.sf <- st_as_sf(midLAAL.rast.stars)
plot(midLAAL.rast.sf)
midLAAL.rast.sf <- midLAAL.rast.sf %>% 
  st_transform(crs = 4326) %>% 
  st_wrap_dateline() %>%
  st_shift_longitude() %>% 
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)
plot(midLAAL.rast.sf, main="lcea with lat_ts=0") # raster is now sf object in lcea

# now filter into different groups
ml <- midLAAL.rast.sf
ml95 <- ml %>% filter(n<95.00001)
ml50 <- ml %>% filter(between(n,10.0000001,50.00001))
ml10 <- ml %>% filter(n<10.0000001)

#############################
######### PLOTTING ##########
#############################

# Generate basemap 

### TWO METHODS ###
# Ptolemy package, formerly "nPacMaps". Very amazing for North Pacific Mapping

# First method, generate a basemap based on the extent and crs of data
# full resolution (slow)
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/basemap/")
load("npac_base.Rdata")

# for faster test plotting
load("npac_base_res.Rdata")

# the actual one you want to use
load("npac_base_i.Rdata")

# set up necessary color gradients
bluecols <- brewer.pal(9, 'Blues')
pie(rep(1,9), col = bluecols)
newcol <- colorRampPalette(bluecols)
ncols <- 163
bluecols2 <- newcol(ncols)
bluecols2 <- rev(bluecols2)
pie(rep(1, ncols), col = bluecols2, border = NA, labels = NA)

redcols <- brewer.pal(9, 'Reds')
pie(rep(1,9), col = redcols)
newcol <- colorRampPalette(redcols)
ncols <- 163
redcols2 <- newcol(ncols)
redcols2 <- rev(redcols2)
pie(rep(1, ncols), col = redcols2, border = NA, labels = NA)

# ggplot call
# Legends are too difficult to sort out here - just plotting enough to manipulate in Adobe Illustrator

figure <- ggplot() + 
  # LAAL contours
  geom_sf(data=ml95c, aes(fill=cntr_level)) +
  geom_sf(data=ml50c, aes(fill=cntr_level)) +
  geom_sf(data=ml10c, aes(fill=cntr_level)) +
  # base map and other parameters
  geom_sf(data=npac_base_i) +
  viridis::scale_fill_viridis(direction = -1, discrete=T)+
  guides(fill=guide_legend(title="Contour level"))+
  coord_sf(xlim = c(-4537510, 6536980), ylim = c(1463885, 6141532)) +
  theme_bw()+
  ggtitle("a")+
  theme(plot.title = element_text(size=12))+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8))+ #change legend text font size
  scalebar(x.min=6547510,x.max=953690,y.min=1803885,y.max=6141532, location = "bottomright", dist = 1500,
           dist_unit = "km", transform = FALSE, height = 0.03, st.size=3, st.dist = 0.05)+
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.7, "in"), pad_y = unit(0.34, "in"),
                         style = north_arrow_fancy_orienteering,
                         height = unit(.75,"cm"),
                         width = unit(.75,"cm"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

figure

path = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/midLAAL/data_to_load/figure_midLAAL.Rdata"
save(figure, file=path)
setwd(path)
load()















