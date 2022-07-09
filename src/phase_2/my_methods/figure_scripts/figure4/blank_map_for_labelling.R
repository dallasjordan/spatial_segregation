# Script to make a blank map to add location names in Adobe Illustrator

# Dallas Jordan
# May 18 2022

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

npac_base_res <- npac_base_res %>%
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)


# Plot  -------------------------------------------------------------------

figure <- ggplot() + 
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
