# Plotting contours from adehabitathr - used for publication figures
# Most up to date mapping script June 21 2021
# Dallas Jordan

# adapting master_mapping script to plot my rasters and contours for overallLAAL/overallBFAL
# notes and organization is better here than master_mapping script. Will need to adapt that script or just use this one
# primarily. Currently afraid to touch that script because it is able to produce spp/island comparison maps using ks
# generated shape files

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
  setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_contours/")
# for pc
  setwd("E:/project_data/spatial_segregation/figures/allLAAL_allBFAL/master_script_contours")

  load("vert95_allLAAL.Rdata")
  load("vert95_allBFAL.Rdata")
  load("vert50_allLAAL.Rdata")
  load("vert50_allBFAL.Rdata")
  load("vert10_allLAAL.Rdata")
  load("vert10_allBFAL.Rdata")

# LAAL contours
al95c <- st_as_sf(vert95_allLAAL) # All LAAL 95th UD Contour
al95c <- al95c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

al50c <- st_as_sf(vert50_allLAAL) # All LAAL 50th UD Contour
al50c <- al50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)


al10c <- st_as_sf(vert10_allLAAL) # All LAAL 10th UD Contour
al10c <- al10c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)

# BFAL contours

ab95c <- st_as_sf(vert95_allBFAL) # All BFAL 95th UD Contour
ab95c <- ab95c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)

ab50c <- st_as_sf(vert50_allBFAL) # All BFAL 50th UD Contour
ab50c <- ab50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)


ab10c <- st_as_sf(vert10_allBFAL) # All BFAL 10th UD Contour
ab10c <- ab10c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)

###############################
### IMPORT AND PREP RASTERS ###
###############################

# load rasters
# for mac
  setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/master_script_rasters/")
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
  st_transform(crs = 4326) %>% 
  st_wrap_dateline() %>%
  st_shift_longitude() %>% 
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)
plot(allLAAL.rast.sf, main="lcea with lat_ts=0") # raster is now sf object in lcea

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
  st_transform(crs = 4326) %>% 
  st_wrap_dateline() %>%
  st_shift_longitude() %>% 
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)
plot(allBFAL.rast.sf, main="lcea with lat_ts=0") # raster is now sf object in lcea

# now filter into different groups
ab <- allBFAL.rast.sf
ab95 <- ab %>% filter(n<95.00001)
ab50 <- ab %>% filter(between(n,10.0000001,50.00001))
ab10 <- ab %>% filter(n<10.0000001)

#############################
######### PLOTTING ##########
#############################

# Generate basemap 

### TWO METHODS ###
# Ptolemy package, formerly "nPacMaps". Very amazing for North Pacific Mapping

# First method, generate a basemap based on the extent and crs of data
npac_base<- ptolemy::extract_gshhg(allBFAL.rast.sf, resolution = "h", epsg = NULL, buffer = 5000,
                                    simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS

# for faster test plotting
npac_base_res <- ptolemy::extract_gshhg(allBFAL.rast.sf, resolution = "l", epsg = NULL, buffer = 5000,
                                        simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS

npac_base_i <- ptolemy::extract_gshhg(allBFAL.rast.sf, resolution = "i", epsg = NULL, buffer = 5000,
                                        simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS

# Second method, baked-in function for PDC mercator Npac basemap
npac_base <- ptolemy::npac()

npac_plot <- ggplot() + 
  geom_sf(data = npac_base,
          fill = "grey60", size = 0.2) +
  ggtitle('North Pacific Basemap (epsg:LCEA)')
npac_plot
            
            # Second method again, used built-in npac basemap and transform to the projection I need (default for this one is
            # WGS84 / PDC Mercator, but I need LCEA)
            npac_base2 <- ptolemy::npac(epsg = lcea) # EPSG code 6933 is Lambert Cylindrical Equal Area. This function will throw some errors, just ignore. Default is some other CRS
            
            npac_plot <- ggplot() + 
              geom_sf(data = npac_base,
                      fill = "grey60", size = 0.2) +
              ggtitle('North Pacific Basemap (epsg:LCEA)')+
              coord_sf(xlim = c(-4537510, 6536980), ylim = c(1063885, 5654736)) 
            npac_plot

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

# blank for adding placenames
figure <- ggplot() + 
  # base map and other parameters
  geom_sf(data=npac_base_i) +
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 12000000)) +
  theme_bw()+
  ggtitle("a")+
  theme(plot.title = element_text(size=12))+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8))+ #change legend text font size
  scalebar(x.min=9447510,x.max=953690,y.min=1673885,y.max=6141532, location = "bottomright", dist = 1500,
           dist_unit = "km", transform = FALSE, height = 0.03, st.size=3, st.dist = 0.05)+
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.5, "in"), pad_y = unit(1, "in"),
                         style = north_arrow_fancy_orienteering,
                         height = unit(.75,"cm"),
                         width = unit(.75,"cm"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
figure


# LAAL RASTER
figure1 <- ggplot() + 
  # LAAL raster
  geom_sf(data=al95, aes(fill=n, colors="transparent")) +
  scale_fill_gradientn(colors=bluecols2,na.value = "transparent",
                      breaks=c(0,50,100),labels=c(0,50,100),
                      limits=c(0,100), name="LAAL UD%")+
  new_scale("fill")+
  # BFAL raster
  # geom_sf(data=ab95, aes(fill=n, colors="transparent")) +
  # scale_fill_gradientn(colors=redcols2,na.value = "transparent",
  #                       breaks=c(0,50,100),labels=c(0,50,100),
  #                      limits=c(0,100), name="BFAL UD%")+
  # LAAL contours
  # geom_sf(data=al95c, color=("#87CEFF"), fill=alpha("#87CEFF",0.2)) +
  # geom_sf(data=al50c, color=("blue"), fill=alpha("blue",0.6)) +
  # geom_sf(data=al10c, color=("darkblue"), alpha=0.8, fill="darkblue") +
  # BFAL contours
  # geom_sf(data=ab95c, color=("#E39191"), fill=alpha("#E39191",0.2)) +
  # geom_sf(data=ab50c, color=("#E86464"), fill=alpha("#E86464",0.6)) +
  # geom_sf(data=ab10c, color=("#E30303"), alpha=0.8, fill="#E30303") +
  # base map and other parameters
  geom_sf(data=npac_base_i) +
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 12000000)) +
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

figure1


# BFAL RASTER
figure2 <- ggplot() + 
  # BFAL raster
  geom_sf(data=ab95, aes(fill=n, colors="transparent")) +
  scale_fill_gradientn(colors=redcols2,na.value = "transparent",
                        breaks=c(0,50,100),labels=c(0,50,100),
                       limits=c(0,100), name="BFAL UD%")+
# base map and other parameters
  geom_sf(data=npac_base_i) +
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 12000000)) +
  theme_bw()+
  ggtitle("d")+
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
  geom_sf(data=al95c, color=("#87CEFF"), size=1.1, fill=alpha("#87CEFF",0.2)) +
  geom_sf(data=al50c, color=("blue"), size=1.1, fill=alpha("blue",0.6)) +
  # geom_sf(data=al10c, color=("darkblue"), alpha=0.8, fill="darkblue") +
  # BFAL contours
  geom_sf(data=ab95c, color=("#E39191"), size=1.1, fill=alpha("#E39191",0.2)) +
  geom_sf(data=ab50c, color=("#E86464"), size=1.1, fill=alpha("#E86464",0.6)) +
  # geom_sf(data=ab10c, color=("#E30303"), alpha=0.8, fill="#E30303") +
  # base map and other parameters
  geom_sf(data=npac_base_i) +
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 12000000)) +
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

figure3





# cowplot::plot_grid(figure, figure2, figure3, labelsize=12, ncol=1, nrow=3)



path = "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/figures/allLAAL_allBFAL/data_to_load/figure_allLAALBFAL_A.Rdata"
save(figure, file=path)
setwd(path)
load()






vvc("#E39191", "#E86464", "#E30303")
c("#FF0000", "#EE0000", "#CD0000")
colourpicker:::colourPickerAddin()
c("#87CEFF", "#FFFFFF", "#FFFFFF")c("#FFFFFF")

#############################
###### STUDY SITE MAP #######
#############################

# Making study site map 

midway <- c(182.63, 28.2)
midway <- fortify(midway)
tern <- c(193.716, 23.87)
islands <- data.frame(ID = c("midway","tern"),
                      x = c(182.63,193.716),
                      y = c(28.2, 23.87))

islands = fortify(islands)


figure <- ggplot() + 
  geom_point(data=islands, aes(x=x,y=y),col="red", size=5) +
  geom_sf(data=sf_worldmp)+ #I just called a discrete color scale here. Custom colors can get really tricky in sf when you have multiple layers.
  coord_sf(xlim = c(120, 240), ylim = c(10, 70)) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()+

# change colors 
figure












###############################################################################################################################################################

###############################################
### PROBABLY DON'T NEED ANYTHING BELOW THIS ###
###############################################

##################################################################################
### OLD METHOD OF LOADING IN RASTERS (that were exported from another script) ###
##################################################################################

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures")
load('LAALdata_tern.Rdata')
data_tern<-LAAL

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures")
load('LAALdata_midway.Rdata')
data_midway<-LAAL

data <- rbind(data_tern,data_midway)# load LAAL point data
LAAL <- data

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/kde/combined/LAAL/")
LAALshape <- readOGR(".","KDE_combined_LAAL") # read shapefile in, change to wherever it is stored, don't have .shp at the end of pathway

################################
### TRANSFORMING DATA FOR SF ###
################################

# USEFUL REFERENCE FOR THE FUTURE: 

# These are relics - don't need these if using the Ptolemy North Pacific basemaps
sf_worldmp <- st_as_sf(maps::map("world2",plot=FALSE,fill=TRUE)) #basemap
sf_worldmp <- sf_worldmp %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(lcea)
sf_LAALpts <- st_as_sf(LAAL, coords=c("x","y"), crs=4326) #points
sf_kde <- st_as_sf(LAALshape) #KDE, original projection

sf_kde_fin <- sf_kde %>% #here I call sf_kde, my starting point. 
  st_transform("+proj=longlat +datum=WGS84") %>% #transform to WGS84.
  st_wrap_dateline() %>% #wrap around the dateline
  st_shift_longitude() %>% #shift. somehow this shifts the longitude to where we need it to be. I don't understand the details. 
  st_union(by_feature = TRUE) #union. We want this by feature so we can include that part of the call in here. And then union by feature combines those features back to their original shape. If you just run this above line without the st_union call there is a vertical line at the date line.

################################
### MANIPULATING BBOX IN SF  ###
################################
# test_bbox <- npac_bbox
# old_bbox <- ptolemy::npac_bbox
# new_bb = c(-2231390, 2253424, 11136950, 21121680)
# names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
# attr(new_bb, "class") = "bbox"
# attr(test_bbox, "bbox") = new_bb


#############################
### NOT SURE WHAT THIS IS ###
#############################

###NOTE: projecting this into WGS effs it up and makes it wrap around weird.
#LAALshape <- spTransform(LAALshape,CRS("+proj=longlat +datum=WGS84 +no_defs"))
#LAALshape <- st_make_valid(LAALshape) ## This is for SF objects!
#LAALshape <- maptools::nowrapRecenter(LAALshape) #this is for spatial objects
#LAALshape <- spTransform(LAALshape,CRS("+proj=longlat +lon_0=-150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +no_defs")) ## Not the projection we want
#LAALshape2 <- fortify(LAALshape) ##Not relevant in the sf space
###This is the old plot with the spatialpoints/lines/polygons
ggplot() + 
  geom_polygon(data = fortify(maps::map("world2",plot=FALSE,fill=TRUE)), aes(x=long, y = lat, group=group)) +
  coord_fixed(1.3) + 
  coord_map(projection = "laea") +
  coord_equal(xlim = c(120, 240), ylim = c(10, 70)) +
  geom_point(data=LAAL, aes(x=x,y=y), size=0.1) 
# geom_path(data=LAALshape2, aes(x=long, y=lat), size=3, color="red")



##### WHAT YOU ORIGINALLY HAD BEFORE YOU HAD TO FIX THE WEIRD COLORINGS/LAYERING ISSUE
figure <- ggplot() + 
  # geom_sf(data=sf_LAALpts, size=0.8) +
  geom_sf(data=sf_kde_fin, aes(fill=rev(cntr_level)), alpha=0.8) +
  geom_sf(data=sf_worldmp)+
  viridis::scale_fill_viridis(discrete=T)+ #I just called a discrete color scale here. Custom colors can get really tricky in sf when you have multiple layers.
  coord_sf(xlim = c(120, 240), ylim = c(10, 70)) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()

# change colors 
figure


#############################
### NOT SURE WHAT THIS IS ###
#############################


# export the figure!
import <- file.path("LAAL_combined_pointdens_raster.Rdata")
stars_r <- read_stars("import", RasterIO = list(nBufXSize = 600, nBufYSize = 600), proxy=T)
test_stars <- st_as_stars()
figure + geom_stars(data=test_spdf)+coord_equal() +
  facet_wrap(~band) +
  theme_void() +
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))

figure <- ggplot() + 
  geom_tile(data=test_df, aes(x=x,y=y,fill = value, alpha=0.8)) +
  scale_fill_gradient(low = 'white', high = 'blue') +
  #coord_sf(xlim = c(120, 240), ylim = c(10, 70)) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()

figure

gplot(r) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'white', high = 'blue') +
  coord_equal()

p <- ggplot() +
  geom_tile(data=test_df, aes(x=x,y=y,fill = value)) +
  scale_fill_viridis_c(option = "B", direction = -1) +
  geom_sf(data=sf_kde1, aes(color=(cntr_level)), alpha=0.8, fill=NA) +
  labs(title = "Likelihood of swinging and missing on a fastball",
       y = "spin rate (rpm)") +
  theme_light()

p


#############################
### NOT SURE WHAT THIS IS ###
#############################

# CONTOURS NOT WORKING FEB 22, REGULAR PLOTTING DOES

library(data.table)
library(rnaturalearthdata)
library(sp)
library(sf)
library(ggplot2)

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata.Rdata") # load LAAL point data

LAALshape <- readOGR(".","KDE_midway_LAAL") # read shapefile in, change to wherever it is stored, dont have .shp at the end of pathway
LAALshape <- spTransform(LAALshape,CRS("+proj=longlat +datum=WGS84 +no_defs"))


sf_worldmp <- st_as_sf(maps::map("world2",plot=FALSE,fill=TRUE))
sf_LAALpts <- st_as_sf(LAAL, coords=c("x","y"), crs=4326)
sf_kde <- st_as_sf(LAALshape)
sf_kde <- st_make_valid(sf_kde)

ggplot() + 
  geom_sf(data=sf_worldmp)+
  geom_sf(data=sf_LAALpts) +
  geom_sf(data=sf_kde) +
  coord_sf(xlim = c(120, 240), ylim = c(10, 70)) 


# old, testing another way: 

# LAALshape <- st_make_valid(LAALshape)
# LAALshape <- maptools::nowrapRecenter(LAALshape)
# 
# LAALshape <- spTransform(LAALshape,CRS("+proj=longlat +lon_0=-150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +no_defs"))
# LAALshape2 <- fortify(LAALshape)
# 
# ggplot() + 
#   geom_polygon(data = fortify(maps::map("world2",plot=FALSE,fill=TRUE)), aes(x=long, y = lat, group=group)) +
#   coord_fixed(1.3) + 
#   coord_map(projection = "laea") +
#   coord_equal(xlim = c(120, 240), ylim = c(10, 70)) +
#   geom_point(data=LAAL, aes(x=x,y=y), size=0.1) 
#   # geom_path(data=LAALshape2, aes(x=long, y=lat), size=3, color="red")

#############################
# ANOTHER MAPPING REFERENCE #
#############################
ggplot(atmos, aes(lon, lat)) +
  geom_world() +
  geom_contour(aes(z = slp, color = ..level..), binwidth = 4) +
  scale_color_viridis_c("Sea level pressure") +
  
  new_scale_color() +   # same as `new_scale("color")`
  
  geom_contour(aes(z = air, color = ..level..), binwidth = 4) +
  scale_color_distiller("Air Temperature", palette = "Spectral")  +
  
  scale_x_longitude(limits = c(-150, 0)) +
  scale_y_latitude(ticks = 15) +
  ggalt::coord_proj("+proj=moll +lon_0=-75", 
                    ylim = c(-60, 0), xlim = c(-150, 0))




