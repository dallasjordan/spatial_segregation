# Figure 4 script
# based on overalLAALBFAL_mapping.R script, adapted again so that I can standardize figures easier. Trying to make the aesthetics
# match figure 2 and figure 3. 
# Dallas Jordan Feb 28 2022

### JES EDITS: I made changes to L43: project basemap to eqc. To all your longitude shifts starting on 
# L90, I transformed them to eqc at the end. For your rasters, I just transformed them to eqc and kept 
# them as rasters. I added a bounding box at the start of the plotting section and calculated xmin ymin 
# xmax ymax. In your ggplot code, I used layer_spatial for the rasters. I edited the range of the color 
# gradient to go from 0,95. I added a ton to the coord_sf to plot everything in eqc with the same 
# minmax and used ndiscr to get two x axes. Couldn't solve getting you longitude labels every 20 degrees 
# or so as that is a major problem with sf.

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

midway <- data.frame(longitude = -177.3761,
                     latitude = 28.2101)
tern <- data.frame(longitude = -166.284,
                   latitude = 23.870)

midway <- st_as_sf(midway, coords = c("longitude", "latitude"), 
                   crs = 4326, agr = "constant")
tern <- st_as_sf(tern, coords = c("longitude", "latitude"), 
                 crs = 4326, agr = "constant")

# load in basemap 
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/basemap/")
load("npac_base_res.Rdata")

npac_base_res_wgs <- st_transform(npac_base_res, crs=eqc) #Set to eqc to match the rest of the files

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
  st_transform(crs=4326) %>%
  st_shift_longitude() %>%
  st_transform(crs=eqc)

al50c <- st_as_sf(vert50_allLAAL) # All LAAL 50th UD Contour
al50c <- al50c %>% 
  st_transform(crs=4326) %>%
  st_shift_longitude() %>%
  st_transform(crs=eqc)

al10c <- st_as_sf(vert10_allLAAL) # All LAAL 10th UD Contour
al10c <- al10c %>% 
  st_transform(crs=4326) %>%
  st_shift_longitude() %>%
  st_transform(crs=eqc)

# BFAL contours

ab95c <- st_as_sf(vert95_allBFAL) # All BFAL 95th UD Contour
ab95c <- ab95c %>% 
  st_transform(crs=4326) %>%
  st_shift_longitude() %>%
  st_transform(crs=eqc)

ab50c <- st_as_sf(vert50_allBFAL) # All BFAL 50th UD Contour
ab50c <- ab50c %>% 
  st_transform(crs=4326) %>%
  st_shift_longitude() %>%
  st_transform(crs=eqc)


ab10c <- st_as_sf(vert10_allBFAL) # All BFAL 10th UD Contour
ab10c <- ab10c %>% 
  st_transform(crs=4326) %>%
  st_shift_longitude() %>%
  st_transform(crs=eqc)

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
### JES NOTES: I've played with the stars package before and it's always messed up my rasters. So instead, just project the raster to eqc and use the layer_spatial() call in ggplot instead. It's so much easier. It took me a LONG time to find this function!!!
allLAAL.rast.jes <- projectRaster(allLAAL_rast, crs=eqc)

# allLAAL.sp <- as(allLAAL.rast.jes, "SpatialPixelsDataFrame")
# gridded(allLAAL.sp)
# spplot(allLAAL.sp, main="raster to sp - SpatialPixelsDataFrame")
# allLAAL.rast.stars <- st_as_stars(allLAAL.sp, att=1)
# allLAAL.rast.stars
# allLAAL.rast.sf <- st_as_sf(allLAAL.rast.stars)
# plot(allLAAL.rast.sf)
# allLAAL.rast.sf <- allLAAL.rast.sf %>%
#   st_transform(crs=4326) %>%
#   st_shift_longitude() %>%
# st_transform(eqc)
# plot(allLAAL.rast.sf, main="lcea with lat_ts=0") # raster is now sf object in PDC mercator
# 
# # now filter into different groups
# al <- allLAAL.rast.sf
# al95 <- al %>% filter(n<95.00001)
# al50 <- al %>% filter(between(n,10.0000001,50.00001))
# al10 <- al %>% filter(n<10.0000001)

### BFAL ###
allBFAL.rast.jes <- projectRaster(allBFAL_rast, crs=eqc)

# allBFAL.sp <- as(allBFAL.rast.jes, "SpatialPixelsDataFrame")
# gridded(allBFAL.sp)
# spplot(allBFAL.sp, main="raster to sp - SpatialPixelsDataFrame")
# allBFAL.rast.stars <- st_as_stars(allBFAL.sp, att=1)
# allBFAL.rast.stars
# allBFAL.rast.sf <- st_as_sf(allBFAL.rast.stars)
# plot(allBFAL.rast.sf)
# allBFAL.rast.sf <- allBFAL.rast.sf %>% 
#   st_transform(crs=4326) %>%
#   st_shift_longitude() %>%
# st_transform(crs=eqc)
# plot(allBFAL.rast.sf, main="lcea with lat_ts=0") # raster is now sf object in PDC mercator
# 
# # now filter into different groups
# ab <- allBFAL.rast.sf
# ab95 <- ab %>% filter(n<95.00001)
# ab50 <- ab %>% filter(between(n,10.0000001,50.00001))
# ab10 <- ab %>% filter(n<10.0000001)

#############################
######### PLOTTING ##########
#############################

# ggplot call
###JES NOTES: I changed this code around a bit. We want to first make a bounding box so all our 
# extents will be the same. I outline that below here. Then in the gpplot calls, we use 
# layer_spatial() to bring in our raster grid. I also set the "limits" in the colors to c(0,95) 
# to remove the pink shading. Feel free to change as you see fit. Last, I futzed with the coord_sf() 
# call by adding in a crs which will transform everything for the plot which makes it so you don't 
# have to pre-transform all your stuff. But remember all your layers have to be in the same projection 
# for input! Because that no longer uses the nice limits you had laid out before, I cropped it to 
# 10,70 x -120,120 like you had before. I also tried to add x axlis limits but it just did not work at 
# all... seems like it's an underlying sf issue. If you add 

## First, make bounding box in WGS and convert to eqc
## I made this for 10-70 lat x -120 to 120 longitudinal...
## Order of points goes...
## xmax, ymin
## xmax, ymax
## xmin, ymax
## xmin, ymin
## xmax, ymin <- you have to repeat this from the first coords set to "close" the box
your_extent <- st_sfc(st_polygon(list(matrix(c(123,5.278694,
                                               123,78.761556,
                                               -110,78.761556,
                                               -110,5.278694,
                                               123,5.278694), ncol=2, byrow = T))))
st_crs(your_extent) <- 4326 #set the crs to wgs
your_extent <- st_transform(your_extent, crs=eqc) #transform the crs
#ggplot(data=your_extent) + geom_sf() # test to make sure coords are correct

#you will use these next two calls in your coords_sf() call in your ggplots to set the same extent in each plot
xlim_prj <- c(min(st_coordinates(your_extent)[,1]), max(st_coordinates(your_extent)[,1])) #pull xmin and xmax
ylim_prj <- c(min(st_coordinates(your_extent)[,2]), max(st_coordinates(your_extent)[,2])) #pull ymin and ymax
x_breaks <- (abs(xlim_prj[1]) + abs(xlim_prj[2])) / 9 #calculate break width for plotting axis labels
y_breaks <- (abs(ylim_prj[1]) + abs(ylim_prj[2])) / 7

# Legends are too difficult to sort out here - just plotting enough to manipulate in Adobe Illustrator

# LAAL RASTER
figure1 <- ggplot() + 
  # LAAL raster
  #geom_sf(data=al95, aes(fill=n, colors="transparent")) +
  layer_spatial(allLAAL.rast.jes) +
  geom_sf(data=al95c, color=("grey30"), size=1.25, fill=alpha("#E39191",0), inherit.aes = FALSE) +
  geom_sf(data=al50c, color=("black"), size=1.25, fill=alpha("#E86464",0), inherit.aes = FALSE) +
  geom_sf(data=npac_base_res, fill='grey60')+
  scale_y_continuous(breaks = c(10, 30, 50, 70))+
  scale_x_continuous(breaks = c(-120,-140,-160,180,160,140,120))+
  geom_sf(data=midway, size=3,shape=17,fill="black",color="black")+
  geom_sf(data=tern, size=4,shape=18,fill="black",color="black")+
  scale_fill_gradientn(colors=redcols2, na.value = NA,
                       breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,95), name="LAAL UD%")+
  theme_bw()+
  coord_sf(expand=F, crs=st_crs(eqc), xlim=xlim_prj, ylim=ylim_prj, ndiscr=5) 
figure1

figure2 <- ggplot() + 
  # LAAL raster
  #geom_sf(data=ab95, aes(fill=n, colors="transparent")) +
  layer_spatial(allBFAL.rast.jes) + 
  geom_sf(data=ab95c, color=("grey30"), size=1.25, fill=alpha("#87CEFF",0)) +
  geom_sf(data=ab50c, color=("black"), size=1.25, fill=alpha("blue",0)) +
  geom_sf(data=npac_base_res, fill='grey60')+
  scale_y_continuous(breaks = c(10, 30, 50, 70))+
  scale_x_continuous(breaks = c(-120,-140,-160,180,160,140,120))+
  geom_sf(data=midway, size=3,shape=17,fill="black",color="black")+
  geom_sf(data=tern, size=4,shape=18,fill="black",color="black")+
  scale_fill_gradientn(colors=bluecols2, na.value = NA,
                       breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,95), name="LAAL UD%")+
  theme_bw()+
  coord_sf(expand=F, crs=st_crs(eqc), xlim=xlim_prj, ylim=ylim_prj, ndiscr=5) 
figure2


# BFAL RASTER
figure2 <- ggplot() + 
  # BFAL raster
  #geom_sf(data=ab95, aes(fill=n, colors="transparent")) +
  layer_spatial(allLAAL.rast.jes) +
  scale_fill_gradientn(colors=bluecols2,na.value = "transparent",
                       breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,100), name="BFAL UD%")+
  # base map and other parameters
  geom_sf(data=npac_base_res_wgs) +
  coord_sf(expand=F, crs=st_crs(eqc), xlim = xlim_prj, ylim=ylim_prj, ndiscr=5) +
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
  geom_sf(data=npac_base_res_wgs) +
  coord_sf(expand=F, crs= st_crs(eqc), xlim=xlim_prj, ylim=ylim_prj, ndiscr=5)+
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