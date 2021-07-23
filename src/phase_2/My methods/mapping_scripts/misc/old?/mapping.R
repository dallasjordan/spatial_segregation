# plotting testing

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

# SOME IMPORTANT NOTES: 
# sf is most up to date spatial system
# tmap is one mapping package; you can use ggmap
# To make an sf object from a dataframe:
st_as_sf(my_data, coords=c("x","y"))

# generate dataset

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
calculate_sp_obj_extent <- function(sp_obj){
  x_range <- sp_obj@bbox[1,2] - sp_obj@bbox[1,1]
  y_range <- sp_obj@bbox[2,2] - sp_obj@bbox[2,1]
  x_range_km <- x_range/1000
  y_range_km <- y_range/1000
  print1<- paste0("the x-range in km is ",x_range_km)
  print2<- paste0("the y-range in km is ",y_range_km)
  grid_calc <- x_range_km/300
  print3<- paste0("for a grid size of 300km^2 use grid parameter ",grid_calc)
  print(print1)
  print(print2)
  print(print3)
  return(grid_calc)
}
years <- c("2008","2009","2010","2011","2012")

###################################################
#COMBINE ALL INDIVIDUALS IN A YEAR, DON'T PRESERVE INDIVIDUAL ID:

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  colnames(data_loop)= c("dtime","x","y")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","x","y")]


# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  colnames(data_loop)= c("dtime","x","y")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","x","y")]


# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)

data.spL <- split(data.sp, data.sp$id)

LAAL <-rbind(data.spL[["LAAL_2008"]], data.spL[["LAAL_2009"]], data.spL[["LAAL_2010"]], data.spL[["LAAL_2011"]], data.spL[["LAAL_2012"]])
BFAL <-rbind(data.spL[["BFAL_2008"]], data.spL[["BFAL_2009"]], data.spL[["BFAL_2010"]], data.spL[["BFAL_2011"]], data.spL[["BFAL_2012"]])


### Use adehabitatHR to make contour lines
pdc_mercator_proj<-sf::st_crs(3832)
testLAAL <- LAAL[,c("x","y")]
testLAAL <- sp::SpatialPoints(testLAAL, proj4string = CRS("+proj=longlat +datum=WGS84"))
testLAAL <- sp::spTransform(testLAAL,pdc_mercator_proj$proj4string)

kernel.ref <- kernelUD(testLAAL)
image(kernel.ref)
ud <- getverticeshr(kernel.ref, percent = 95, unin = "m", unout = "km2") 
plot(ud)
df <- fortify(ud)

### Doing everything with ggplot2
# Internally, uses MASS::kde2d()
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

g <- ggplot(df, aes(x = long, y = lat, fill = id, group = group)) +
      geom_polygon(alpha = .4) +
      coord_equal() +
      theme_void()
g

g2 <- ggplot() +
        geom_sf(data=world) +
        geom_point(data=LAAL, aes(x,y)) + # you have geom_points with a geom_sf - does this work? Are the projections the same?
        geom_density_2d(data = LAAL, aes(x,y))
        coord_sf(crs = "+proj=laea +lat_0=10 +lon_0=-165 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs ") 
      
g2

g3 <- ggplot() +
        geom_sf(data=world) +
        geom_point(data=LAAL, aes(x,y)) + # you have geom_points with a geom_sf - does this work? Are the projections the same?
        geom_density_2d(data = LAAL, aes(x,y)) +
        coord_sf(crs = "+proj=longlat +zone=2 +datum=WGS84 +no_defs") 

g3

########
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

#################
# THIS IS IT. THIS IS YOUR SOLUTION!

mp1 <- fortify(map(fill=TRUE, plot=FALSE))
mp2 <- mp1
mp2$long <- mp2$long + 360
mp2$group <- mp2$group + max(mp2$group) + 1
mp <- rbind(mp1, mp2)

LAALcoords <- LAAL[,c("x", "y")]
kd <- ks::kde(LAALcoords, compute.cont=TRUE)
contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["5%"])[[1]])
contour_95 <- data.frame(contour_95)

BFALcoords <- BFAL[,c("x","y")]
kd2 <- ks::kde(BFALcoords, compute.cont=TRUE)
BFALcontour_95 <- with(kd2, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                      z=estimate, levels=cont["5%"])[[1]])
BFALcontour_95 <- data.frame(BFALcontour95)

best2 <- ggplot() + 
            geom_polygon(aes(x = long, y = lat, group = group), data = mp, fill="grey80") + 
            scale_x_continuous(limits = c(130, 240)) + 
            scale_y_continuous(limits = c(10, 75)) +
            # geom_point(data=LAALcoords, aes(x,y)) +
            geom_path(aes(x, y), data=contour_95) +
            geom_path(aes(x,y,), data=BFALcontour_95, color = "blue") + 
            # geom_point(data=BFALcoords, aes(x,y)) + 
            theme_bw()

best2


require(sp)
require(maptools)
IDs <- sapply(strsplit(mp$region, ":"), function(x) x[1])


poly <- map2SpatialPolygons(mp, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))





library(rgeos)
library(ggspatial)
# gene world map
ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(120.00, 235.00), ylim = c(0, 70), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()



# example of pacific centered map in base R
        for (i in 1:length(LAAL$y)){
          if (LAAL$y[i]<0){
            LAAL$y[i] = LAAL$y[i]+360
          }
        }
        data(world2HiresMapEnv)
        graphics::plot(LAAL$x,LAAL$y,type = "p", ylab="Latitude", xlab="Longitude",
             xlim= c(120,240), ylim = c(-10, 75), bty = "n")
        map('world2Hires', xlim= c(120,240), ylim = c(-10, 75),  col = "grey", add = T)


###########################################################################################################
# Mapping

# hypothetical steps: 
# convert data to sf
# use ggplot( data = sf object)
# use the st_shift_longitude() function

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)



ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("test map") +
  coord_sf(xlim = c(-200, -140.12), ylim = c(7.65, 33.97), expand = FALSE, crs = "+proj=laea +lat_0=10 +lon_0=-165 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs ") 

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-200, -140.12), ylim = c(7.65, 33.97), expand = FALSE)


Lon180to360 <- function(lon){
  lon %% 360
}

library(maptools)
plot(Lon180to360(m$lon), Lon180to360(m$lat), ylim=c(minlat, maxlat), c(minlon, maxlon)) # you will need to define max and min lons and lats for your extent
plot(wrld_simpl, add = TRUE, col="grey")
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE, col="grey")







# Not sure if this will ever work: 

pdc_mercator_proj<-sf::st_crs(3832)
LEA <- sf::st_crs(3035)
coordinates(LAAL) <- c("x","y")
proj4string(LAAL) <- CRS("+proj=laea +lat_0=10 +lon_0=-165 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")
#LAAL <- spTransform(current_LAAL,pdc_mercator_proj$proj4string)

coordinates(BFAL) <- c("x","y")
proj4string(BFAL) <- CRS("+proj=laea +lat_0=10 +lon_0=-165 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")
#BFAL<- spTransform(current_BFAL,pdc_mercator_proj$proj4string)

grid_input <- calculate_sp_obj_extent(LAAL)
kernel.ref <- kernelUD(LAAL, same4all = T, grid=38)
data.kernel.poly <- getverticeshr(kernel.ref, percent = 90, unin = "m", unout = "km2") 
data.kernel.poly$id
data.kernel.poly$area

kernel.ref <- kernelUD(current_BFAL, same4all = T)
data.kernel.poly <- getverticeshr(kernel.ref, percent = 90, unin = "m", unout = "km2") 
data.kernel.poly$id
data.kernel.poly$area
