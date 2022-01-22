# replace section 3a with this to get working aeqd plot

# Midway LAAL x Tern LAAL
# Midway LAAL
c3a_midLAAL_points <- CRW_sim_allpoints %>% filter(grepl('LAAL_MID', Animal_ID))
c3a_midLAAL_points <- c3a_midLAAL_points %>% filter(iter==1)
# over_pole <- function(x){
#   if_else(x>6363885, x-(x-6363885),x)
# }
# over_pole <- function(x,y){
#   if(x>6363885){
#     x<-x-(x-6363885);
#     y<-
#   } 
# }
# c3a_midLAAL_points$Latitude <- over_pole(c3a_midLAAL_points$Latitude, c3a_midLAAL_points$Longitude)

# fix latitude to cap out at 90 N
max_lat <- function(x){
  if_else(x>6339452,6339452,x) 
}
c3a_midLAAL_points$Latitude <- max_lat(c3a_midLAAL_points$Latitude)
ind <- with(c3a_midLAAL_points, (Latitude == 6339452))
c3a_midLAAL_points <- c3a_midLAAL_points[!ind,]

#convert to sf for easy ggplot2 plotting
c3a_midLAAL_points_sf <- st_as_sf(c3a_midLAAL_points, coords = c("Longitude","Latitude"), remove = TRUE, na.fail = TRUE, crs=lcea)
c3a_midLAAL_points_sf <- c3a_midLAAL_points_sf %>%
  st_transform(3349)
# 
# g=st_sfc(st_point(c(-425223.3,6303822)),crs=lcea)
# g <- g %>% st_transform(crs=3349)
# 
# g=st_sfc(st_point(c(-5843408.58,6363885)),crs=lcea)
# g <- g %>% st_transform(crs=4326)
# g
# 
# g=st_sfc(st_point(c(180,49)),crs=4326)
# g <- g%>% st_transform(crs=lcea)
# g
# 


c3a_midLAAL_points <- c3a_midLAAL_points[,c(2,3,4)]
sp::coordinates(c3a_midLAAL_points) <- c("Longitude", "Latitude")
proj4string(c3a_midLAAL_points) <- CRS(lcea) # placeholder
c3a_midLAAL_points <- spTransform(c3a_midLAAL_points,aeqd)
grid_input <- calculate_sp_obj_extent(c3a_midLAAL_points,0.1)

c3a_midLAAL_ud <-  kernelUD(c3a_midLAAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)

# average into 1 UD
mid_iter_holder <- numeric(length = 1764)
for (d in 1:length(c3a_midLAAL_ud)) {
  c3a_midLAAL_ud[[d]]@data$ud[is.na(c3a_midLAAL_ud[[d]]@data$ud)] <- 0
  add <- c3a_midLAAL_ud[[d]]@data$ud
  mid_iter_holder <- mid_iter_holder+add
}
mid_iter_holder <- mid_iter_holder/length(c3a_midLAAL_ud)
mid_iter_holder

##### modify existing estUD object with averaged values, then rename
mid_iter_avg <- c3a_midLAAL_ud[[1]]
mid_iter_avg@data$ud <- mid_iter_holder
mid_iter_avg@data$ud[is.na(mid_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
image(mid_iter_avg)
CRW_ml95c <- getverticeshr.estUD(mid_iter_avg,percent = 95)
# take SpatialPolygonsDataframe, reproject...this is where it fails
CRW_ml95c  <- st_as_sf(CRW_ml95c)
st_is_valid(CRW_ml95c)
# sf2 <- st_transform(CRW_ml95c, crs = "+proj=longlat +datum=WGS84" )
# st_is_valid(sf2, reason=T)
# sf2 <- st_make_valid(sf2$geometry)
# sf2 <- sf2 %>%
#   st_wrap_dateline() %>% # wrap around the dateline
#   st_shift_longitude() %>%
#   st_union(by_feature = TRUE) %>%
#   st_transform(crs = 3349)
# 
# 
# CRW_ml95c  <- CRW_ml95c %>%
#   st_transform(crs = 4326) 
# CRW_ml95c <- spTransform(CRW_ml95c,CRS("+proj=longlat +datum=WGS84 +no_defs"))
# CRW_ml95c  <- st_as_sf(CRW_ml95c) # 95th UD Contour
# CRW_ml95c <- st_make_valid(CRW_ml95c$geometry)
# 
# CRW_ml95c  <- CRW_ml95c %>%
#   st_transform(crs = 4326) %>% # transform to WGS84.
#   st_wrap_dateline() %>% # wrap around the dateline
#   st_shift_longitude() %>%
#   st_union(by_feature = TRUE) %>%
#   st_transform(crs = 3349)

# Tern LAAL
c3a_ternLAAL_points <- CRW_sim_allpoints %>% filter(grepl('LAAL_TERN', Animal_ID))
c3a_ternLAAL_points <- c3a_ternLAAL_points %>% filter(iter==1) 

c3a_ternLAAL_points$Latitude <- max_lat(c3a_ternLAAL_points$Latitude)
ind <- with(c3a_ternLAAL_points, (Latitude == 6339452))
c3a_ternLAAL_points <- c3a_ternLAAL_points[!ind,]

#convert to sf for easy ggplot2 plotting
c3a_ternLAAL_points_sf <- st_as_sf(c3a_ternLAAL_points, coords = c("Longitude","Latitude"), remove = TRUE, na.fail = TRUE, crs=lcea)
c3a_ternLAAL_points_sf <- c3a_ternLAAL_points_sf %>%
  st_transform(aeqd)

c3a_ternLAAL_points <- c3a_ternLAAL_points[,c(2,3,4)]
sp::coordinates(c3a_ternLAAL_points) <- c("Longitude", "Latitude")
proj4string(c3a_ternLAAL_points) <- CRS(lcea) # placeholder
c3a_ternLAAL_points <- spTransform(c3a_ternLAAL_points,aeqd)
grid_input <- calculate_sp_obj_extent(c3a_ternLAAL_points,0.1)

c3a_ternLAAL_ud <-  kernelUD(c3a_ternLAAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)

# average into 1 UD
tern_iter_holder <- numeric(length = 4148)
for (d in 1:length(c3a_ternLAAL_ud)) {
  c3a_ternLAAL_ud[[d]]@data$ud[is.na(c3a_ternLAAL_ud[[d]]@data$ud)] <- 0
  add <- c3a_ternLAAL_ud[[d]]@data$ud
  tern_iter_holder <- tern_iter_holder+add
}
tern_iter_holder <- tern_iter_holder/length(c3a_ternLAAL_ud)
tern_iter_holder

##### modify existing estUD object with averaged values, then rename
tern_iter_avg <- c3a_ternLAAL_ud[[1]]
tern_iter_avg@data$ud <- tern_iter_holder
tern_iter_avg@data$ud[is.na(tern_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
image(tern_iter_avg)
CRW_tl95c <- getverticeshr(tern_iter_avg,percent = 95)

CRW_tl95c  <- st_as_sf(CRW_tl95c)
st_is_valid(CRW_tl95c)


# CRW_tl95c  <- CRW_tl95c %>%
#   st_transform(crs = 4326) %>% # transform to WGS84.
#   st_wrap_dateline() %>% # wrap around the dateline
#   st_shift_longitude() %>%
#   st_union(by_feature = TRUE) %>%
#   st_transform(crs = 3832)
# 
# plot c3a
npac_base_i <- ptolemy::extract_gshhg(c3a_ternLAAL_points_sf, resolution = "c", epsg = NULL, buffer = 5000,
                                      simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS

midLAALternLAALCRWrandom <- ggplot() + 
  # base map and other parameters
  #geom_sf(data=c3a_midLAAL_points_sf, color="orange",size=0.005)+
  #geom_sf(data=c3a_ternLAAL_points_sf, color="firebrick",size=0.005)+
  geom_sf(data=CRW_ml95c, color="orange",fill="orange",alpha=0.5)+
  geom_sf(data=CRW_tl95c, color="firebrick",fill="firebrick",alpha=0.5)+
  geom_sf(data=npac_base_i, fill="grey60")
#coord_sf(xlim = c(-8000000, 8000000), ylim = c(0, 13000000),expand=F)
midLAALternLAALCRWrandom

setwd("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3")
save(midLAALternLAALCRWrandom,file="3a.Rdata")
load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/3a.Rdata")
midLAALternLAALCRWrandom