# MIDWAY LAAL 
# Plotting contours from adehabitathr - used for publication figures
# June 30 2021
# Dallas Jordan
# Last updated : July 24 2021
  # July 24 2021:
  # I thought this was the script that I needed for my manuscript Figure 2 - it wasn't, this was another version that plotted similar info to Figure 1 script (overalLAALBFAL). 
  # Others in the XXXX_mapping type are not updated with PDCmerc projection. 
  # Dec 18 2021: 
  # I'm returning to the XXXX_mapping type to remake the Figure 2(4 panel plot of individual contours) to make colors
  # consistent with other figures in my manuscript. On July 24, I'm not sure why I thought this wasn't what was needed 
  # for figure 2: it is the only script that makes those images in LCEA projections



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
# for mac
  setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/individual/midLAAL/master_script_contours/")
  load("vert95_midLAAL.Rdata")
  load("vert50_midLAAL.Rdata")
  load("vert10_midLAAL.Rdata")
# for pc
  setwd("E:/project_data/spatial_segregation/figures/individual/midLAAL/master_script_contours/")
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
  st_transform(crs = 3832)

ml50c <- st_as_sf(vert50_midLAAL) # 50th UD Contour
ml50c <- ml50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)


ml10c <- st_as_sf(vert10_midLAAL) # 10th UD Contour
ml10c <- ml10c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

###############################
### IMPORT AND PREP RASTERS ###
###############################

# load rasters
# for mac
  setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/individual/midLAAL/master_script_rasters/")
# for pc 
  setwd("E:/project_data/spatial_segregation/figures/individual/midLAAL/master_script_rasters/")
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
  st_transform(crs = 3832)
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
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/basemap/")
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
  geom_sf(data=ml95, aes(fill=n, colors="transparent")) +
  scale_fill_gradientn(colors=bluecols2,na.value = "transparent",
                       breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,100), name="LAAL UD%")+
  new_scale("fill")+
  # LAAL contours
  geom_sf(data=ml95c, color=("#87CEFF"), fill=alpha("#87CEFF",0.2)) +
  geom_sf(data=ml50c, color=("blue"), fill=alpha("blue",0.6)) +
  #geom_sf(data=ml10c, color=("darkblue"), alpha=0.8, fill="darkblue") +
  # base map and other parameters
  geom_sf(data=npac_base_i) +
  coord_sf(xlim = c(-4537510, 6536980), ylim = c(1463885, 6141532)) +
  theme_bw()+
  ggtitle("Midway LAAL Density add contour legend")+
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

path = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/individual/midLAAL/data_to_load/figure_midLAAL.Rdata"
save(figure, file=path)
setwd(path)
load()















