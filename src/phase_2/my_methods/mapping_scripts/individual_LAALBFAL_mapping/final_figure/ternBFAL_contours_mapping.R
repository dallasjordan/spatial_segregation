# TERN BFAL 
# Plotting contours from adehabitathr - used for publication figures
# June 30 2021
# Dallas Jordan
# Last updated July 24 2021

# adapting overallLAALBFAL_mapping script to plot my rasters and contours for Tern BFAL
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

#############################
### IMPORT CONTOUR LINES  ###
#############################

#load contours, as of Jun 22 they are SpatialPolygonsDataFrame
# for mac
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternBFAL/master_script_contours/")
load("vert95_ternBFAL.Rdata")
load("vert50_ternBFAL.Rdata")
load("vert10_ternBFAL.Rdata")
# for pc
setwd("E:/project_data/spatial_segregation/figures/individual/ternBFAL/master_script_contours/")
load("vert95_ternBFAL.Rdata")
load("vert50_ternBFAL.Rdata")
load("vert10_ternBFAL.Rdata")

# Contours in 'sf'
tb95c <- st_as_sf(vert95_ternBFAL) # 95th UD Contour
tb95c <- tb95c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

tb50c <- st_as_sf(vert50_ternBFAL) # 50th UD Contour
tb50c <- tb50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)


tb10c <- st_as_sf(vert10_ternBFAL) # 10th UD Contour
tb10c <- tb10c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

# additions July 7 to get viridis scale working
tb95c$id <- "foraging_range"
tb95c$cntr_level <- "95"

tb50c$id <- "focal_range"
tb50c$cntr_level <- "50"

tb10c$id <- "core_range"
tb10c$cntr_level <- "10"

#############################
######### PLOTTING ##########
#############################

# load in npac basemap made elsewhere (ptolemy package)
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/basemap/")
load("npac_base.Rdata")

# for faster test plotting
load("npac_base_res.Rdata")

# the actual one you want to use
load("npac_base_i.Rdata")

# second method
# on pc, loading in data object from overallLAALBFAL_mapping to get the same proportions as my other figure for this manuscript
    setwd("E:/project_data/spatial_segregation/figures/allLAAL_allBFAL/master_script_rasters/")
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
      st_transform(crs = 3832)
    plot(allBFAL.rast.sf, main="PDCmerc with lat_ts=0") # raster is now sf object in PDC mercator
    ab <- allBFAL.rast.sf
    
    npac_base_i <- ptolemy::extract_gshhg(allBFAL.rast.sf, resolution = "i", epsg = NULL, buffer = 5000,
                                          simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS

# ggplot call - minimal

figure <- ggplot() + 
  # BFAL contours
  # if you want borders: 
  #   geom_sf(data=tb95c, aes(color=cntr_level, fill=cntr_level)) +
  geom_sf(data=tb95c, aes(fill=cntr_level)) +
  geom_sf(data=tb50c, aes(fill=cntr_level)) +
  geom_sf(data=tb10c, aes(fill=cntr_level)) +
  # base map and other parameters
  geom_sf(data=npac_base_i) +
  viridis::scale_fill_viridis(direction = -1, discrete=T)+
  # if you want borders: 
  #   viridis::scale_color_viridis(direction = -1, discrete=T)+
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 12000000)) +
  theme_bw()+
  ggtitle("d")+
  theme(plot.title = element_text(size=12))+
  theme(legend.position = "none")+
  # theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
  #       legend.key.height = unit(0.5, 'cm'), #change legend key height
  #       legend.key.width = unit(0.5, 'cm'), #change legend key width
  #       legend.title = element_text(size=10), #change legend title font size
  #       legend.text = element_text(size=8))+ #change legend text font size
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

figure


path = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternBFAL/data_to_load/figure_ternBFAL.Rdata"
save(figure, file=path)
setwd(path)
load()