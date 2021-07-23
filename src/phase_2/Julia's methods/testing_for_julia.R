library(ggmap)
library(mapdata)
library(maptools)
library(sp)
library(adehabitatHR)
library(ggplot2)
theme_set(theme_bw())
library(sf)
library('rnaturalearth')
library('rnaturalearthdata')
library(rgeos)
library(ggspatial)
library(lwgeom)


#load LAALdata_for_julia.Rdata, this is a giant data frame of all my modeled LAAL positions across 5 years

###########################
# BEST SOLUTION SO FAR: 
library(ggplot2)
library(mapdata)

mp1 <- fortify(map(fill=TRUE, plot=FALSE))
mp2 <- mp1
mp2$long <- mp2$long + 360
mp2$group <- mp2$group + max(mp2$group) + 1
mp <- rbind(mp1, mp2)
best<- ggplot() + 
  geom_path(aes(x = long, y = lat, group = group), data = mp) + 
  scale_x_continuous(limits = c(120, 240)) + 
  scale_y_continuous(limits = c(-10, 75)) +
  geom_density_2d(data = LAAL, aes(x,y))


# issues with this: 
# contour levels are generated internally by geom_density_2d, which relies 
# on the MASS:kde2d() function. I can't easily change these levels to correspond 
# to what would be the 50%. 75%, 90%, and 95% contour lines of the density. 

###########################
#using adehabitatHR to make home ranges. Here, I can specify specific contour lines

### Use adehabitatHR to make contour lines
pdc_mercator_proj<-sf::st_crs(3832)
testLAAL <- LAAL[,c("x","y")]
testLAAL <- sp::SpatialPoints(testLAAL, proj4strinßßg = CRS("+proj=longlat +datum=WGS84"))
testLAAL <- sp::spTransform(testLAAL,pdc_mercator_proj$proj4string)

kernel.ref <- kernelUD(testLAAL)
image(kernel.ref)
ud <- getverticeshr(kernel.ref, percent = 95, unin = "m", unout = "km2") # can do more levels here
plot(ud)
df <- fortify(ud)

# ok now there is a dataframe of the contour lines...but issues:
  # it looks weird. don't really know if it is accuract
  # data fed into this is in proj4string = CRS("+proj=longlat +datum=WGS84") - no idea if this is an issue or not?
  # what do I do with it at this point? Convert to spatialLines object and plot with ggplot?


##########################
# another plotting attempt: 

sf_extSoftVersion()
# generate world map


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

world <- as(world, "Spatial")
world <- maptools::nowrapRecenter(world)
world <- st_as_sf(world)

world = st_intersection(world, st_make_valid(world))
# world<- st_shift_longitude(world)

ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(120.00, 235.00), ylim = c(0, 70), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()

mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))

ggplot() +
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group)) 

# This works but has weird polygon fill
###########
# https://www.r-bloggers.com/2019/04/zooming-in-on-maps-with-sf-and-ggplot2/

worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
zoom_to <- c(177.37, 28.2)  # ~ center of Midway
zoom_level <- 1.5

# Lambert azimuthal equal-area projection around center of interest
target_crs <- sprintf('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                      zoom_to[1], zoom_to[2])

C <- 40075016.686   # ~ circumference of Earth in meters
x_span <- 360 / 2^zoom_level
y_span <- 180 / 2^(zoom_level+1)


zoom_to_xy <- st_transform(st_sfc(st_point(zoom_to), crs = 4326),
                           crs = target_crs)
zoom_to_xy

disp_window <- st_sfc(
  st_point(st_coordinates(zoom_to_xy - c(x_span / 2, y_span / 2))),
  st_point(st_coordinates(zoom_to_xy + c(x_span / 2, y_span / 2))),
  crs = target_crs
)

ggplot() + geom_sf(data = worldmap) +
  geom_sf(data = zoom_to_xy, color = 'red') +
  coord_sf(xlim = st_coordinates(disp_window)[,'X'],
           ylim = c(10,70),
           crs = target_crs, datum = target_crs) +
  theme_bw()

# issues: 
# cut off...won't map the pacific. If I could get this working to show the whole pacific, then I could take the contours from the 
# previous chunk, try to convert to Spatial Lines, and then plot them here?

#####################################################
# another plotting attempt: 
  sf_extSoftVersion()
# generate world map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
# world <- as(world, "Spatial")
# world <- maptools::nowrapRecenter(world)
world <- st_as_sf(world)
world = st_intersection(world, st_make_valid(world))
world<- st_shift_longitude(world)
ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(120.00, 235.00), ylim = c(0, 70), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()

# seems broken

##################################################

# https://www.r-bloggers.com/2019/04/zooming-in-on-maps-with-sf-and-ggplot2/
worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
zoom_to <- c(177.37, 28.2)  # ~ center of Midway
zoom_level <- 1.5
# Lambert azimuthal equal-area projection around center of interest
target_crs <- sprintf('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                      zoom_to[1], zoom_to[2])
C <- 40075016.686   # ~ circumference of Earth in meters
x_span <- 360 / 2^zoom_level
y_span <- 180 / 2^(zoom_level+1)
zoom_to_xy <- st_transform(st_sfc(st_point(zoom_to), crs = 4326),
                           crs = target_crs)
zoom_to_xy
disp_window <- st_sfc(
  st_point(st_coordinates(zoom_to_xy - c(x_span / 2, y_span / 2))),
  st_point(st_coordinates(zoom_to_xy + c(x_span / 2, y_span / 2))),
  crs = target_crs
)
ggplot() + geom_sf(data = worldmap) +
  geom_sf(data = zoom_to_xy, color = 'red') +
  coord_sf(xlim = st_coordinates(disp_window)[,'X'],
           ylim = c(10,70),
           crs = target_crs, datum = target_crs) +
  theme_bw()
