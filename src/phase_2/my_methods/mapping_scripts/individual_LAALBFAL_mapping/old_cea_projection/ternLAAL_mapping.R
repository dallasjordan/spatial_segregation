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
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/individual/ternLAAL/master_script_contours/")
load("vert95_ternLAAL.Rdata")
load("vert50_ternLAAL.Rdata")
load("vert10_ternLAAL.Rdata")

# Contours in 'sf'
tl95c <- st_as_sf(vert95_ternLAAL) # 95th UD Contour
tl95c <- tl95c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)

tl50c <- st_as_sf(vert50_ternLAAL) # 50th UD Contour
tl50c <- tl50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)


tl10c <- st_as_sf(vert10_ternLAAL) # 10th UD Contour
tl10c <- tl10c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)

###############################
### IMPORT AND PREP RASTERS ###
###############################

# load rasters
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/individual/ternLAAL/master_script_rasters/")
load("ternLAAL_ud_vol_rast.Rdata")
ternLAAL_rast <- ternLAAL.ud.vol.raster
image(ternLAAL_rast)
contour(ternLAAL_rast,levels=c(95,50,10),add=T)

### LAAL ###

ternLAAL.sp <- as(ternLAAL_rast, "SpatialPixelsDataFrame")
gridded(ternLAAL.sp)
spplot(ternLAAL.sp, main="raster to sp - SpatialPixelsDataFrame")
ternLAAL.rast.stars <- st_as_stars(ternLAAL.sp, att=1)
ternLAAL.rast.stars
ternLAAL.rast.sf <- st_as_sf(ternLAAL.rast.stars)
plot(ternLAAL.rast.sf)
ternLAAL.rast.sf <- ternLAAL.rast.sf %>% 
  st_transform(crs = 4326) %>% 
  st_wrap_dateline() %>%
  st_shift_longitude() %>% 
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)
plot(ternLAAL.rast.sf, main="lcea with lat_ts=0") # raster is now sf object in lcea

# now filter into different groups
tl <- ternLAAL.rast.sf
tl95 <- tl %>% filter(n<95.00001)
tl50 <- tl %>% filter(between(n,10.0000001,50.00001))
tl10 <- tl %>% filter(n<10.0000001)

#############################
######### PLOTTING ##########
#############################

# Generate basemap 

### TWO METHODS ###
# Ptolemy package, formerly "nPacMaps". Very amazing for North Pacific Mapping

# First method, generate a basemap based on the extent and crs of data
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/basemap/")
load("npac_base.Rdata")

# for faster test plotting
load("npac_base_res.Rdata")

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
  # LAAL raster
  geom_sf(data=tl95, aes(fill=n, colors="transparent")) +
  scale_fill_gradientn(colors=bluecols2,na.value = "transparent",
                       breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,100), name="LAAL UD%")+
  new_scale("fill")+
  # LAAL contours
  geom_sf(data=tl95c, color=("#87CEFF"), fill=alpha("#87CEFF",0.2)) +
  geom_sf(data=tl50c, color=("blue"), fill=alpha("blue",0.6)) +
  geom_sf(data=tl10c, color=("darkblue"), alpha=0.8, fill="darkblue") +
  # base map and other parameters
  geom_sf(data=npac_base) +
  coord_sf(xlim = c(-4537510, 6536980), ylim = c(1463885, 6141532)) +
  theme_bw()+
  ggtitle("Tern LAAL Density add contour legend")+
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

path = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/individual/ternLAAL/data_to_load/figure_ternLAAL.Rdata"
save(figure, file=path)
setwd(path)
load()















