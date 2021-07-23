# Plotting contours from adehabitathr - used for publication figures
# Most up to date mapping script as of July 2021 is the overallLAALBFAL_mapping and individual mapping scripts
# This script is preserved because it can handle 'ks' objects generated using the script Ellie helped me adapt
# Potentially useful but not used for my Master's thesis
# Dallas Jordan



# for reading in the map, when you load the rnaturalearthdata library it already has all these 
# maps in there. If you need something specific and customized, there is code to download those 
# (like oceans data specifically, or islands, etc) from the naturalearthdata library
# So for your whole read zip remove zip etc, it's just here...

sf::sf_extSoftVersion()
# Put this in terminal if you need to add GEOS support: 
# R CMD INSTALL sf_0.9-7.tar.gz --configure-args='--with-gdal-config=/Library/Frameworks/GDAL.framework/Versions/3.2/unix/bin/gdal-config -with-geos-config=/Library/Frameworks/GEOS.framework/Versions/3D/unix/bin/geos-config'
#######
# WORKING

library(mapproj)
library(data.table)
library(rnaturalearthdata)
library(sf)
sf::sf_use_s2(FALSE)
# The above line is necessary as of April 17 2021, when 'sf' switched to using S2 spheroid 
# geometry instead of a GEOS routine that assumed projected coordinates. Until 'mapping' 
# updates its back-end, switching S2 geometry off is necessary for my script to work. 
# See: https://github.com/r-spatial/sf/issues/1649
# Great reference here: https://r-spatial.github.io/sf/articles/sf7.html
library(rgdal)
library(ggspatial)
library(ggmap)
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

###Now let's try it in the sf space
#first, convert everything to sf objects
sf_worldmp <- st_as_sf(maps::map("world2",plot=FALSE,fill=TRUE)) #basemap
sf_LAALpts <- st_as_sf(LAAL, coords=c("x","y"), crs=4326) #points
sf_kde <- st_as_sf(LAALshape) #KDE, original projection

# these might be helpful one day
#   st_wrap_dateline(c("WRAPDATELINE=YES", "DATELINEOFFSET=40")) %>%
#   st_zm(drop = TRUE) drops a "z" or "m" attribute if its erroneous

sf_kde_fin <- sf_kde %>% #here I call sf_kde, my starting point. I use the %>% sign to indicate that I'm gonna keep doing processing in this pipe
  st_transform("+proj=longlat +datum=WGS84") %>% #transform to WGS84.
  st_wrap_dateline() %>% #wrap around the dateline
  st_shift_longitude() %>% #shift. somehow this shifts the longitude to where we need it to be. I don't understand the details. 
  st_union(by_feature = TRUE) #union. We want this by feature so we can include that part of the call in here. And then union by feature combines those features back to their original shape. If you just run this above line without the st_union call there is a vertical line at the date line.

##### fucked up shit to reorder geom_sf by layer

# CHANGE THE ROWS CORRESPONDING TO EACH UD LEVEL!!

sf_kde1 <- sf_kde_fin[1:3,]
sf_kde2 <- sf_kde_fin[4,]
sf_kde3 <- sf_kde_fin[5,]

# add raster
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures")
load("LAAL_combined_pointdens_raster.Rdata")
r_df <- as.data.frame(r, xy = TRUE) 
test_spdf <- as(r, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")
test_sf<-st_as_sf(test_df,as_points=FALSE)

figure <- ggplot() + 
  # geom_sf(data=sf_LAALpts, size=0.8) +
  geom_sf(data=sf_kde1, aes(color=(cntr_level)), alpha=0.8, fill=NA) +
  geom_sf(data=sf_kde2, aes(color=(cntr_level)), alpha=0.8, fill=NA) +
  geom_sf(data=sf_kde3, aes(color=(cntr_level)), alpha=0.8, fill=NA) +
  geom_sf(data=sf_worldmp)+
  viridis::scale_fill_viridis(direction = -1, discrete=T)+ #I just called a discrete color scale here. Custom colors can get really tricky in sf when you have multiple layers.
  coord_sf(xlim = c(120, 240), ylim = c(10, 70)) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()

figure

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
  # geom_sf(data=sf_LAALpts, size=0.8) +
  # geom_sf(data=sf_kde1, aes(color=(cntr_level)), fill=NA) +
  # geom_sf(data=sf_kde2, aes(color=(cntr_level)), fill=NA) +
  # geom_sf(data=sf_kde3, aes(color=(cntr_level)), fill=NA) +
  #geom_sf(data=sf_worldmp)+
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
  theme_bw()

# change colors 
figure


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




  
