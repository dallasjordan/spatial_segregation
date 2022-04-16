# Land mask creation and testing script for the CRW simulations of real data. I created this script
# to help make an appropriate land mask solution for filtering out CRW points that fall over land, 
# in my CRW_real_data script.
# Jan 29 2022
# Dallas Jordan

library(sf)
sf::sf_use_s2(FALSE)
library(dplyr)

# Solution chosen, using Ptolemy to get a land mask. I set a custom bbox to encompass a large region
# of the North Pacific. Anything outside of this will get filtered by unreasonable lat/lons. 
land_mask_bbox<-st_bbox(c(xmin = -3008182, xmax = 10844089, ymax = 18027535, ymin = 185132), crs = st_crs(3832))
land_mask_sfc <- st_as_sfc(land_mask_bbox)
land_mask <- ptolemy::extract_gshhg(land_mask_sfc, resolution = "i", epsg = NULL, buffer = 5000,
                                    simplify = FALSE) 
land_mask$type = "land"

# testing
plot(land_mask)
point <- st_as_sfc(c("POINT(-160.476561 57.766494)")) %>% # this point is over ocean
  st_sf(type = 'point')
st_crs(point)<- 4326
point <- point %>% st_transform(crs=3832)
st_join(point, land_mask, join = st_intersects)

point <- st_as_sfc(c("POINT(-157.428771 59.749996)")) %>% # this point is over land
  st_sf(type = 'point')
st_crs(point)<- 4326
point <- point %>% st_transform(crs=3832)
st_join(point, land_mask, join = st_intersects)

proj4string(CRS("+init=epsg:3832")) 

# Other unused solutions to create land mask ------------------------------
# using rnaturalearthdata polygon
land <- st_read("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/environment_data/ne_110m_land/ne_110m_land.shp")
land <- land %>% 
  rename(id='featurecla') %>% 
  dplyr::select(id) 
plot(st_geometry(land))

point <- st_as_sfc(c("POINT(-160.476561 57.766494)"))%>% # this point is over ocean
  st_sf(type = 'point')
st_crs(point)<- crs(land)

in_land <- st_join(point, land, join = st_intersects)
in_land

# using world from spData
library(spData)
world <- world %>% dplyr::select('name_long')
st_join(point, world, join = st_intersects)
plot(st_geometry(world))
plot(st_geometry(point), pch=16,add=T)

# original code used for testing st_join
poly <- st_as_sfc(c("POLYGON((0 0 , 0 1 , 1 1 , 1 0, 0 0))")) %>% 
  st_sf(ID = "poly1")    

pts <- st_as_sfc(c("POINT(0.5 0.5)",
                   "POINT(0.6 0.6)",
                   "POINT(3 3)")) %>%
  st_sf(ID = paste0("point", 1:3))

st_join(pts, poly, join = st_intersects)




# figuring out lat/lon bounds of reasonable points in PDC mercator using bounds of actual data
point <- st_as_sfc(c("POINT(238.43096 10.65519)")) %>% # this point is over ocean
  st_sf(type = 'point')
st_crs(point)<- 4326
point <- point %>% st_transform(crs=3832)
point

# xmin : -1708182
# xmax: 9844089
# ymin: 1185132
# ymax: 9310913

point <- st_as_sfc(c("POINT(-3000000 185132)")) %>% # this point is over ocean
  st_sf(type = 'point')
st_crs(point)<- 3832
point <- point %>% st_transform(crs=4326)
point


# what are the bounds of my actual data?
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
all_data_sf <- st_as_sf(all_data,coords=c("x","y"),crs=4326)
st_bbox(all_data_sf)
# bbox for all_data:
# xmin : 134.65514 
# xmax: 238.43096
# ymin: 10.65519
# ymax: 63.99840 


