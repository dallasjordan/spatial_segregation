# Figure 4 script
# based on overalLAALBFAL_mapping.R script, adapted again so that I can standardize figures easier. Trying to make the aesthetics
# match figure 2 and figure 3. 
# Dallas Jordan Feb 28 2022

library(sf)
library(rgdal)
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

# projections 
lcea <- "+proj=cea +lat_0=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
eqc <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# Reference for this is EPSG 6933, WGS 84 / NSIDC EASE-Grid 2.0 Global
#    cea <- "+proj=cea +lat_0=20 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs" # lat_ts (Standard Parallel) describes latitudes where there is no distortion. 
#    _ts means "true scale", the lat or lon where there is no distortion. lon_0 = the natural origin for longitude

# load in basemap 
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/basemap/")
load("npac_base_res.Rdata")

# set up necessary color gradients for plotting
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

#############################
### IMPORT CONTOUR LINES  ###
#############################

#load contours, as of Jun 22 they are SpatialPolygonsDataFrame
# for mac
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/allLAAL_allBFAL/master_script_contours/")
# for pc
setwd("E:/project_data/spatial_segregation/figures/allLAAL_allBFAL/master_script_contours")

load("vert95_allLAAL.Rdata")
load("vert95_allBFAL.Rdata")
load("vert50_allLAAL.Rdata")
load("vert50_allBFAL.Rdata")
load("vert10_allLAAL.Rdata")
load("vert10_allBFAL.Rdata")

# LAAL contours
# One example left for reference of the old way I used to transform these:
# al95c <- st_as_sf(vert95_allLAAL) # All LAAL 95th UD Contour
# al95c <- al95c %>%
#   st_transform(crs = 4326) %>% # transform to WGS84.
#   st_wrap_dateline() %>% # wrap around the dateline
#   st_shift_longitude() %>%
#   st_union(by_feature = TRUE) %>%
#   st_transform(eqc)

al95c <- st_as_sf(vert95_allLAAL) # All LAAL 95th UD Contour
al95c <- al95c %>%
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

al50c <- st_as_sf(vert50_allLAAL) # All LAAL 50th UD Contour
al50c <- al50c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

al10c <- st_as_sf(vert10_allLAAL) # All LAAL 10th UD Contour
al10c <- al10c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

# BFAL contours

ab95c <- st_as_sf(vert95_allBFAL) # All BFAL 95th UD Contour
ab95c <- ab95c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

ab50c <- st_as_sf(vert50_allBFAL) # All BFAL 50th UD Contour
ab50c <- ab50c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)


ab10c <- st_as_sf(vert10_allBFAL) # All BFAL 10th UD Contour
ab10c <- ab10c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

###############################
### IMPORT AND PREP RASTERS ###
###############################

# load rasters
# for mac
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/allLAAL_allBFAL/master_script_rasters/")
# for pc
setwd("E:/project_data/spatial_segregation/figures/allLAAL_allBFAL/master_script_rasters/")
load("allBFAL_ud_vol_rast.Rdata")
load("allLAAL_ud_vol_rast.Rdata")
allLAAL_rast <- allLAAL.ud.vol.raster
allBFAL_rast <- allBFAL.ud.vol.raster
image(allLAAL_rast)
image(allBFAL_rast)
plot(allBFAL_rast)
contour(allBFAL_rast,levels=c(95,50,10),add=T)
contour(allLAAL_rast,levels=c(95,50,10),add=T)

### LAAL ###

allLAAL.sp <- as(allLAAL_rast, "SpatialPixelsDataFrame")
gridded(allLAAL.sp)
spplot(allLAAL.sp, main="raster to sp - SpatialPixelsDataFrame")
allLAAL.rast.stars <- st_as_stars(allLAAL.sp, att=1)
allLAAL.rast.stars
allLAAL.rast.sf <- st_as_sf(allLAAL.rast.stars)
plot(allLAAL.rast.sf)
allLAAL.rast.sf <- allLAAL.rast.sf %>%
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)
plot(allLAAL.rast.sf, main="lcea with lat_ts=0") # raster is now sf object in PDC mercator

# now filter into different groups
al <- allLAAL.rast.sf
al95 <- al %>% filter(n<95.00001)
al50 <- al %>% filter(between(n,10.0000001,50.00001))
al10 <- al %>% filter(n<10.0000001)

### BFAL ###

allBFAL.sp <- as(allBFAL_rast, "SpatialPixelsDataFrame")
gridded(allBFAL.sp)
spplot(allBFAL.sp, main="raster to sp - SpatialPixelsDataFrame")
allBFAL.rast.stars <- st_as_stars(allBFAL.sp, att=1)
allBFAL.rast.stars
allBFAL.rast.sf <- st_as_sf(allBFAL.rast.stars)
plot(allBFAL.rast.sf)
allBFAL.rast.sf <- allBFAL.rast.sf %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)
plot(allBFAL.rast.sf, main="lcea with lat_ts=0") # raster is now sf object in PDC mercator

# now filter into different groups
ab <- allBFAL.rast.sf
ab95 <- ab %>% filter(n<95.00001)
ab50 <- ab %>% filter(between(n,10.0000001,50.00001))
ab10 <- ab %>% filter(n<10.0000001)

#############################
######### PLOTTING ##########
#############################

# ggplot call
# Legends are too difficult to sort out here - just plotting enough to manipulate in Adobe Illustrator

# LAAL RASTER
figure1 <- ggplot() + 
  # LAAL raster
  geom_sf(data=al95, aes(fill=n, colors="transparent")) +
  geom_sf(data=al95c, color=("grey30"), size=1.25, fill=alpha("#E39191",0)) +
  geom_sf(data=al50c, color=("black"), size=1.25, fill=alpha("#E86464",0)) +
  geom_sf(data=npac_base_res, fill='grey60')+
  scale_fill_gradientn(colors=redcols2,na.value = "transparent",
                       breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,100), name="LAAL UD%")+
  theme_bw()+
  coord_sf(expand=F)
figure1

figure2 <- ggplot() + 
  # LAAL raster
  geom_sf(data=ab95, aes(fill=n, colors="transparent")) +
  geom_sf(data=ab95c, color=("grey30"), size=1.25, fill=alpha("#87CEFF",0)) +
  geom_sf(data=ab50c, color=("black"), size=1.25, fill=alpha("blue",0)) +
  geom_sf(data=npac_base_res, fill='grey60')+
  scale_fill_gradientn(colors=bluecols2,na.value = "transparent",
                       breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,100), name="LAAL UD%")+
  theme_bw()+
  coord_sf(expand=F)
figure2


# BFAL RASTER
figure2 <- ggplot() + 
  # BFAL raster
  geom_sf(data=ab95, aes(fill=n, colors="transparent")) +
  scale_fill_gradientn(colors=bluecols2,na.value = "transparent",
                       breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,100), name="BFAL UD%")+
  # base map and other parameters
  geom_sf(data=npac_base_i) +
  coord_sf(expand=F)+
  # for PDC mercator: coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 12000000)) +
  theme_bw()+
  ggtitle("b")+
  theme(plot.title = element_text(size=12))+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8))+ #change legend text font size
  # scalebar(x.min=9947510,x.max=953690,y.min=1803885,y.max=6141532, location = "bottomright", dist = 1500,
  #          dist_unit = "km", transform = FALSE, height = 0.03, st.size=3, st.dist = 0.05)+
  # annotation_north_arrow(location = "br", which_north = "true",
  #                        pad_x = unit(0.5, "in"), pad_y = unit(1, "in"),
  #                        style = north_arrow_fancy_orienteering,
  #                        height = unit(.75,"cm"),
  #                        width = unit(.75,"cm"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
figure2

# LAAL AND BFAL CONTOURS, 50th UD and 95th UD
figure3 <- ggplot() + 
  # LAAL contours
  geom_sf(data=ab95c, color=("#87CEFF"), size=1, fill=alpha("#87CEFF",0.2)) +
  geom_sf(data=ab50c, color=("blue"), size=1, fill=alpha("blue",0.6)) +
  # geom_sf(data=al10c, color=("darkblue"), alpha=0.8, fill="darkblue") +
  # BFAL contours
  geom_sf(data=al95c, color=("#E39191"), size=1, fill=alpha("#E39191",0.2)) +
  geom_sf(data=al50c, color=("#E86464"), size=1, fill=alpha("#E86464",0.6)) +
  # geom_sf(data=ab10c, color=("#E30303"), alpha=0.8, fill="#E30303") +
  # base map and other parameters
  geom_sf(data=npac_base_i) +
  coord_sf(expand=F)+
  # for PDC mercator: coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 12000000)) +
  theme_bw()+
  ggtitle("c")+
  theme(plot.title = element_text(size=12))+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8))+ #change legend text font size
  # scalebar(x.min=9947510,x.max=953690,y.min=1803885,y.max=6141532, location = "bottomright", dist = 1500,
  #          dist_unit = "km", transform = FALSE, height = 0.03, st.size=3, st.dist = 0.05)+
  # annotation_north_arrow(location = "br", which_north = "true",
  #                        pad_x = unit(0.5, "in"), pad_y = unit(1, "in"),
  #                        style = north_arrow_fancy_orienteering,
  #                        height = unit(.75,"cm"),
  #                        width = unit(.75,"cm"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

figure3


# Generate basemap 