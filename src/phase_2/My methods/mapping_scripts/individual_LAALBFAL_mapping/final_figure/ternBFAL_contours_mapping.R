# TERN BFAL 
# Plotting contours from adehabitathr - used for publication figures
# June 30 2021
# Dallas Jordan

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
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternBFAL/master_script_contours/")
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
  st_transform(lcea)

tb50c <- st_as_sf(vert50_ternBFAL) # 50th UD Contour
tb50c <- tb50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)


tb10c <- st_as_sf(vert10_ternBFAL) # 10th UD Contour
tb10c <- tb10c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)

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
  coord_sf(xlim = c(-4537510, 6536980), ylim = c(1463885, 6141532)) +
  theme_bw()+
  ggtitle("d")+
  theme(plot.title = element_text(size=12))+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8))+ #change legend text font size
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

figure


path = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternBFAL/data_to_load/figure_ternBFAL.Rdata"
save(figure, file=path)
setwd(path)
load()