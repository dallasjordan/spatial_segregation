# Mapping 9-panel plot showing randomizations and original data. 
# Three cases: Midway LAAL x Tern LAAL - CRW shows strong overlap
#              Midway BFAL x Tern BFAL - CRW shows medium overlap 
#              Tern LAAL x Tern BFAL - CRW agrees with track perm, sig. segregation 
# Plot original data (points), plot one random of track perm, plot one random of CRW, each of these 
# have contours plotted over randomization points
# Ostensibly will be figure 3
# Dallas Jordan
# Dec 23 2021

# Setup -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(sf)
library(adehabitatHR)

calculate_sp_obj_extent <- function(sp_obj, extent_parameter){
  x_range <- sp_obj@bbox[1,2] - sp_obj@bbox[1,1]
  x_min <- sp_obj@bbox[1,1]-(extent_parameter*x_range)
  x_max <- sp_obj@bbox[1,2]+(extent_parameter*x_range)
  ade_x_range <- x_max-x_min
  y_range <- sp_obj@bbox[2,2] - sp_obj@bbox[2,1]
  y_min <- sp_obj@bbox[2,1]-(extent_parameter*y_range)
  y_max <- sp_obj@bbox[2,2]+(extent_parameter*y_range)
  ade_y_range <- y_max-y_min
  x_range_km <- ade_x_range/1000
  y_range_km <- ade_y_range/1000
  print1<- paste0("the x-range in km is ",x_range_km)
  print2<- paste0("the y-range in km is ",y_range_km)
  grid_calc <- x_range_km/300
  print3<- paste0("for a grid size of 300km x 300km use grid parameter ",grid_calc)
  print(print1)
  print(print2)
  print(print3)
  return(grid_calc)
}

# CRITICAL: 
sf::sf_use_s2(FALSE)
# The above line is necessary as of April 17 2021, when 'sf' switched to using S2 spheroid 
# geometry instead of a GEOS routine that assumed projected coordinates. Until 'mapping' 
# updates its back-end, switching S2 geometry off is necessary for my script to work. 
# See: https://github.com/r-spatial/sf/issues/1649
# Great reference here: https://r-spatial.github.io/sf/articles/sf7.html

lcea <- "+proj=cea +lat_0=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

# import points
    # for mac
    load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
    lm <- all_data[grep("lm", all_data$id), ]
    lt <- all_data[grep("lt", all_data$id), ]
    bm <- all_data[grep("bm", all_data$id), ]
    bt <- all_data[grep("bt", all_data$id), ]
    
    # for pc 
    load("G:/Other computers/My MacBook Pro/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
    lm <- all_data[grep("lm", all_data$id), ]
    lt <- all_data[grep("lt", all_data$id), ]
    bm <- all_data[grep("bm", all_data$id), ]
    bt <- all_data[grep("bt", all_data$id), ]

# import contours 
    # for mac
    setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/individual/midLAAL/master_script_contours/")
    files <- dir("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/individual/", recursive=TRUE, full.names=TRUE, pattern="\\.Rdata$")
    result <- lapply(files, load,envir = .GlobalEnv)
    
    #for pc
    setwd("G:/Other computers/My MacBook Pro/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/individual/midLAAL/master_script_contours/")
    files <- dir("G:/Other computers/My MacBook Pro/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/individual/", recursive=TRUE, full.names=TRUE, pattern="\\.Rdata$")
    result <- lapply(files, load,envir = .GlobalEnv)

# Loading in low res base map since final figure gonna be huge
    # for mac
    setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/basemap/")
    load("npac_base_res.Rdata")
    
    #for pc
    setwd("G:/Other computers/My MacBook Pro/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/contours_rasters_figureData/basemap/")
    load("npac_base_res.Rdata")

npac_base_res <- npac_base_res %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3349)

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

tl95c <- st_as_sf(vert95_ternLAAL) # 95th UD Contour
tl95c <- tl95c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

tl50c <- st_as_sf(vert50_ternLAAL) # 50th UD Contour
tl50c <- tl50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

mb95c <- st_as_sf(vert95_midBFAL) # 95th UD Contour
mb95c <- mb95c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

mb50c <- st_as_sf(vert50_midBFAL) # 50th UD Contour
mb50c <- mb50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

tb95c <- st_as_sf(vert95_ternBFAL) # 95th UD Contour
tb95c <- tb95c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

tb50c <- st_as_sf(vert50_ternBFAL) # 50th UD Contour
tb50c <- tb50c %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)


# Column 1 - original points with contours overlapping  -----------------------------------------------
# Midway LAAL x Tern LAAL 
lm_points = st_as_sf(lm, coords = c("x","y"), remove = FALSE, crs=4326)
lm_points <- lm_points %>%
  st_transform(crs=3349)

lt_points = st_as_sf(lt, coords = c("x","y"), remove = FALSE, crs=4326)
lt_points <- lt_points %>%
  st_transform(crs=3349)

midLAALternLAAL <- ggplot() + 
  # base map and other parameters
  #geom_sf(data=lm_points, color="orange",size=0.005)+
  #geom_sf(data=lt_points, color="firebrick",size=0.005)+
  geom_sf(data=ml95c, color="orange",fill="orange",alpha=0.5)+
  geom_sf(data=tl95c, color="firebrick",fill="firebrick",alpha=0.5)+
  geom_sf(data=npac_base_res, fill="grey60")+
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 13000000),expand=F)
midLAALternLAAL

# Midway BFAL x Tern BFAL

bm_points = st_as_sf(bm, coords = c("x","y"), remove = FALSE, crs=4326)
bm_points <- bm_points %>%
  st_transform(crs=3349)

bt_points = st_as_sf(bt, coords = c("x","y"), remove = FALSE, crs=4326)
bt_points <- bt_points %>%
  st_transform(crs=3349)

midBFALternBFAL <- ggplot() + 
  # base map and other parameters
  #geom_sf(data=lm_points, color="orange",size=0.005)+
  #geom_sf(data=lt_points, color="firebrick",size=0.005)+
  geom_sf(data=mb95c, color="turquoise2",fill="turquoise2",alpha=0.5)+
  geom_sf(data=tb95c, color="royalblue4",fill="royalblue4",alpha=0.5)+
  geom_sf(data=npac_base_res, fill="grey60")+
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 13000000),expand=F)
midBFALternBFAL

# Tern LAAL x Tern BFAL
lt_points = st_as_sf(lt, coords = c("x","y"), remove = FALSE, crs=4326)
lt_points <- lt_points %>%
  st_transform(crs=3349)

bt_points = st_as_sf(bt, coords = c("x","y"), remove = FALSE, crs=4326)
bt_points <- bt_points %>%
  st_transform(crs=3349)

ternLAALternBFAL <- ggplot() + 
  # base map and other parameters
  #geom_sf(data=lm_points, color="orange",size=0.005)+
  #geom_sf(data=lt_points, color="firebrick",size=0.005)+
  geom_sf(data=tl95c, color="firebrick",fill="firebrick",alpha=0.5)+
  geom_sf(data=tb95c, color="royalblue4",fill="royalblue4",alpha=0.5)+
  geom_sf(data=npac_base_res, fill="grey60")+
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 13000000),expand=F)
ternLAALternBFAL


# Column 2 - track permutation random examples ----------------------------
# After running once, you can skip to the bottom and just load these in. The code to recreate 
# these random iterations is left here in case you want to grab different randomized
# iterations. HAVE to run permutation_xx_xx.R scripts when called for to process mid_iter_avg
# and mid/tern_nums

# you are getting randomized 95th contours from the permutation_XX_XX.R scripts - you are running the
# randomization procedure like 20 times and just pulling one of the xxx_iter_avg. You run it 
# in the other script to get it in the global environment, then getverticeshr here. For the points, you still 
# need the xxx_iter_avg, because you need to pull the numbers used for that contour to get the matching underlying
# points used to make that contour

# Midway LAAL x Tern LAAL
# Randomized points: 
############ now run other script to get xxx_iter_avg ############

mid_iter_points <- resample_all_tracks[mid_nums]
tern_iter_points <- resample_all_tracks[tern_nums]

mid_iter_points <- as.data.frame(do.call(rbind, mid_iter_points)) 
tern_iter_points <- as.data.frame(do.call(rbind, tern_iter_points)) 

mid_iter_points = st_as_sf(mid_iter_points, coords = c("x","y"), remove = FALSE, crs=4326)
mid_iter_points <- mid_iter_points %>%
  st_transform(crs=3349)

tern_iter_points = st_as_sf(tern_iter_points, coords = c("x","y"), remove = FALSE, crs=4326)
tern_iter_points <- tern_iter_points %>%
  st_transform(crs=3349)

# Randomized contours: 
ml95cr <- getverticeshr(mid_iter_avg)
ml95cr <- st_as_sf(ml95cr ) # 95th UD Contour
ml95cr <- ml95cr %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

tl95cr <- getverticeshr(tern_iter_avg)
tl95cr <- st_as_sf(tl95cr ) # 95th UD Contour
tl95cr <- tl95cr %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

midLAALternLAALtrackpermrandom <- ggplot() + 
  # base map and other parameters
  geom_sf(data=mid_iter_points, color="orange",size=0.005)+
  geom_sf(data=tern_iter_points, color="firebrick",size=0.005)+
  geom_sf(data=ml95cr, color="orange",fill="orange",alpha=0.5)+
  geom_sf(data=tl95cr, color="firebrick",fill="firebrick",alpha=0.5)+
  geom_sf(data=npac_base_res, fill="grey60")+
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 13000000),expand=F)
midLAALternLAALtrackpermrandom

# save this so you don't save to re-generate from scripts
save(midLAALternLAALtrackpermrandom,file="2a.Rdata")
load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/2a.Rdata")
midLAALternLAALtrackpermrandom

# Midway BFAL x Tern BFAL
# Randomized points: 
############ now run other script to get xxx_iter_avg ############

mid_iter_points <- resample_all_tracks[mid_nums]
tern_iter_points <- resample_all_tracks[tern_nums]

mid_iter_points <- as.data.frame(do.call(rbind, mid_iter_points)) 
tern_iter_points <- as.data.frame(do.call(rbind, tern_iter_points)) 

mid_iter_points = st_as_sf(mid_iter_points, coords = c("x","y"), remove = FALSE, crs=4326)
mid_iter_points <- mid_iter_points %>%
  st_transform(crs=3349)

tern_iter_points = st_as_sf(tern_iter_points, coords = c("x","y"), remove = FALSE, crs=4326)
tern_iter_points <- tern_iter_points %>%
  st_transform(crs=3349)

# Randomized contours: 
mb95cr <- getverticeshr(mid_iter_avg)
mb95cr <- st_as_sf(mb95cr ) # 95th UD Contour
mb95cr <- mb95cr %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

tb95cr <- getverticeshr(tern_iter_avg)
tb95cr <- st_as_sf(tb95cr ) # 95th UD Contour
tb95cr <- tb95cr %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

midBFALternBFALtrackpermrandom <- ggplot() + 
  # base map and other parameters
  geom_sf(data=mid_iter_points, color="turquoise2",size=0.005)+
  geom_sf(data=tern_iter_points, color="royalblue4",size=0.005)+
  geom_sf(data=mb95cr, color="turquoise2",fill="turquoise2",alpha=0.5)+
  geom_sf(data=tb95cr, color="royalblue4",fill="royalblue4",alpha=0.5)+
  geom_sf(data=npac_base_res, fill="grey60")+
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 13000000),expand=F)
midBFALternBFALtrackpermrandom

save(midBFALternBFALtrackpermrandom,file="2b.Rdata")
load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/2b.Rdata")
midBFALternBFALtrackpermrandom

# Tern LAAL x Tern BFAL
# Randomized points: 
############ now run other script to get xxx_iter_avg ############

ternLAAL_iter_points <- resample_all_tracks[mid_nums]
ternBFAL_iter_points <- resample_all_tracks[tern_nums]

ternLAAL_iter_points <- as.data.frame(do.call(rbind, ternLAAL_iter_points)) 
ternBFAL_iter_points <- as.data.frame(do.call(rbind, ternBFAL_iter_points)) 

ternLAAL_iter_points = st_as_sf(ternLAAL_iter_points, coords = c("x","y"), remove = FALSE, crs=4326)
ternLAAL_iter_points <- ternLAAL_iter_points %>%
  st_transform(crs=3349)

ternBFAL_iter_points = st_as_sf(ternBFAL_iter_points, coords = c("x","y"), remove = FALSE, crs=4326)
ternBFAL_iter_points <- ternBFAL_iter_points %>%
  st_transform(crs=3349)

# Randomized contours: 
tl95cr <- getverticeshr(mid_iter_avg)
tl95cr <- st_as_sf(tl95cr) # 95th UD Contour
tl95cr <- tl95cr %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

tb95cr <- getverticeshr(tern_iter_avg)
tb95cr <- st_as_sf(tb95cr) # 95th UD Contour
tb95cr <- tb95cr %>% 
  st_transform(crs = 4326) %>% # transform to WGS84.
  st_wrap_dateline() %>% # wrap around the dateline
  st_shift_longitude() %>%  
  st_union(by_feature = TRUE) %>% 
  st_transform(crs = 3832)

ternLAALternBFALtrackpermrandom <- ggplot() + 
  # base map and other parameters
  geom_sf(data=ternLAAL_iter_points, color="firebrick",size=0.005)+
  geom_sf(data=ternBFAL_iter_points, color="royalblue4",size=0.005)+
  geom_sf(data=tl95cr, color="firebrick",fill="firebrick",alpha=0.5)+
  geom_sf(data=tb95cr, color="royalblue4",fill="royalblue4",alpha=0.5)+
  geom_sf(data=npac_base_res, fill="grey60")+
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 13000000),expand=F)
ternLAALternBFALtrackpermrandom

save(ternLAALternBFALtrackpermrandom,file="2c.Rdata")
load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/2c.Rdata")
ternLAALternBFALtrackpermrandom


# Column 3 - CRW randomization examples -----------------------------------
# import points from CRW_sim_allpoints.csv. Pick an iteration for each class you need. Convert those to an sp dataframe of points. Change 
# the projection. Create contours using getverticeshr, which will require running kernelUD on the spdf. 

CRW_sim_allpoints <- read.csv("G:/Other computers/My MacBook Pro/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/julia_simulations/CRW_permutation_results/real_data/CRW_sim_allpoints.csv")

# Midway LAAL x Tern LAAL
  # Midway LAAL
      c3a_midLAAL_points <- CRW_sim_allpoints %>% filter(grepl('LAAL_MID', Animal_ID))
      c3a_midLAAL_points <- c3a_midLAAL_points %>% filter(iter==1) 
      c3a_midLAAL_points <- c3a_midLAAL_points[,c(2,3,4)]
      
      c3a_midLAAL_sf = st_as_sf(c3a_midLAAL_points)
      
      sp::coordinates(c3a_midLAAL_points) <- c("Longitude", "Latitude")
      grid_input <- calculate_sp_obj_extent(c3a_midLAAL_points,0.1)
      
      c3a_midLAAL_ud <-  kernelUD(c3a_midLAAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      mid_iter_holder <- numeric(length = 2622)
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
      CRW_ml95c <- getverticeshr(mid_iter_avg,percent = 95)
      CRW_ml95c  <- st_as_sf(CRW_ml95c) # 95th UD Contour
      CRW_ml95c <- CRW_ml95c %>% st_set_crs(lcea)
      CRW_ml95c  <- CRW_ml95c %>%
        st_transform(crs = 4326) %>% # transform to WGS84.
        st_wrap_dateline() %>% # wrap around the dateline
        st_shift_longitude() %>%
        st_union(by_feature = TRUE) %>%
        st_transform(crs = 3832)
      
  # Tern LAAL
      c3a_ternLAAL_points <- CRW_sim_allpoints %>% filter(grepl('LAAL_TERN', Animal_ID))
      c3a_ternLAAL_points <- c3a_ternLAAL_points %>% filter(iter==1) 
      c3a_ternLAAL_points <- c3a_ternLAAL_points[,c(2,3,4)]
      
      c3a_ternLAAL_sf = st_as_sf(c3a_ternLAAL_points)
      
      sp::coordinates(c3a_ternLAAL_points) <- c("Longitude", "Latitude")
      grid_input <- calculate_sp_obj_extent(c3a_ternLAAL_points,0.1)
      
      c3a_ternLAAL_ud <-  kernelUD(c3a_ternLAAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      tern_iter_holder <- numeric(length = 3304)
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
      CRW_tl95c  <- st_as_sf(CRW_tl95c) # 95th UD Contour
      CRW_tl95c <- CRW_tl95c %>% st_set_crs(lcea)
      # CRW_tl95c  <- CRW_tl95c %>%
      #   st_transform(crs = 4326) %>% # transform to WGS84.
      #   st_wrap_dateline() %>% # wrap around the dateline
      #   st_shift_longitude() %>%
      #   st_union(by_feature = TRUE) %>%
      #   st_transform(crs = 3832)
      # 
      # plot c3a
      midLAALternLAALCRWrandom<- ggplot() + 
        # base map and other parameters
        geom_sf(data=c3a_midLAAL_sf, color="orange",size=0.005)+
        geom_sf(data=c3a_ternLAAL_sf, color="firebrick",size=0.005)
        geom_sf(data=CRW_ml95c, color="orange",fill="orange",alpha=0.5)+
        geom_sf(data=CRW_tl95c, color="firebrick",fill="firebrick",alpha=0.5)+
        #geom_sf(data=npac_base_res, fill="grey60")
        coord_sf(xlim = c(-9e06, 9e06), ylim = c(-5.5e06, 1.25e07),expand=F)
      midLAALternLAALCRWrandom
      
      save(midLAALternLAALCRWrandom,file="3a.Rdata")
      load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/3a.Rdata")
      midLAALternLAALCRWrandom
      
# Midway BFAL x Tern BFAL
    # Midway BFAL
      c3b_midBFAL_points <- CRW_sim_allpoints %>% filter(grepl('BFAL_MID', Animal_ID))
      c3b_midBFAL_points <- c3b_midBFAL_points %>% filter(iter==1) 
      c3b_midBFAL_points <- c3b_midBFAL_points[,c(2,3,4)]
      
      c3b_midBFAL_sf = st_as_sf(c3b_midBFAL_points, coords = c("Longitude","Latitude"), remove = FALSE)
      
      sp::coordinates(c3b_midBFAL_points) <- c("Longitude", "Latitude")
      grid_input <- calculate_sp_obj_extent(c3b_midBFAL_points,0.1)
      
      c3b_midBFAL_ud <-  kernelUD(c3b_midBFAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      mid_iter_holder <- numeric(length = 1598)
      for (d in 1:length(c3b_midBFAL_ud)) {
        c3b_midBFAL_ud[[d]]@data$ud[is.na(c3b_midBFAL_ud[[d]]@data$ud)] <- 0
        add <- c3b_midBFAL_ud[[d]]@data$ud
        mid_iter_holder <- mid_iter_holder+add
      }
      mid_iter_holder <- mid_iter_holder/length(c3b_midBFAL_ud)
      mid_iter_holder
      
      ##### modify existing estUD object with averaged values, then rename
      mid_iter_avg <- c3b_midBFAL_ud[[1]]
      mid_iter_avg@data$ud <- mid_iter_holder
      mid_iter_avg@data$ud[is.na(mid_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
      image(mid_iter_avg)
      CRW_mb95c <- getverticeshr(mid_iter_avg,percent = 95)
      CRW_mb95c  <- st_as_sf(CRW_mb95c ) # 95th UD Contour
      # CRW_ml95c  <- CRW_ml95c %>% 
      #   st_set_crs(3832) %>%
      #   st_transform(crs = 4326) %>% # transform to WGS84.
      #   st_wrap_dateline() %>% # wrap around the dateline
      #   st_shift_longitude() %>%  
      #   st_union(by_feature = TRUE) %>% 
      #   st_transform(crs = 3832)
      
    # Tern BFAL
      c3b_ternBFAL_points <- CRW_sim_allpoints %>% filter(grepl('BFAL_TERN', Animal_ID))
      c3b_ternBFAL_points <- c3b_ternBFAL_points %>% filter(iter==1) 
      c3b_ternBFAL_points <- c3b_ternBFAL_points[,c(2,3,4)]
      
      c3b_ternBFAL_sf = st_as_sf(c3b_ternBFAL_points, coords = c("Longitude","Latitude"), remove = FALSE)
      
      
      sp::coordinates(c3b_ternBFAL_points) <- c("Longitude", "Latitude")
      grid_input <- calculate_sp_obj_extent(c3b_ternBFAL_points,0.1)
      
      c3b_ternBFAL_ud <-  kernelUD(c3b_ternBFAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      tern_iter_holder <- numeric(length = 3240)
      for (d in 1:length(c3b_ternBFAL_ud)) {
        c3b_ternBFAL_ud[[d]]@data$ud[is.na(c3b_ternBFAL_ud[[d]]@data$ud)] <- 0
        add <- c3b_ternBFAL_ud[[d]]@data$ud
        tern_iter_holder <- tern_iter_holder+add
      }
      tern_iter_holder <- tern_iter_holder/length(c3b_ternBFAL_ud)
      tern_iter_holder
      
      ##### modify existing estUD object with averaged values, then rename
      tern_iter_avg <- c3b_ternBFAL_ud[[1]]
      tern_iter_avg@data$ud <- tern_iter_holder
      tern_iter_avg@data$ud[is.na(tern_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
      image(tern_iter_avg)
      CRW_tb95c <- getverticeshr(tern_iter_avg,percent = 95)
      CRW_tb95c  <- st_as_sf(CRW_tb95c ) # 95th UD Contour
      # CRW_tl95c  <- CRW_tl95c %>% 
      #   st_set_crs(3832) %>%
      #   st_transform(crs = 4326) %>% # transform to WGS84.
      #   st_wrap_dateline() %>% # wrap around the dateline
      #   st_shift_longitude() %>%  
      #   st_union(by_feature = TRUE) %>% 
      #   st_transform(crs = 3832)
      
      # plot c3b
      midBFALternBFALCRWrandom<- ggplot() + 
        # base map and other parameters
        geom_sf(data=c3b_midBFAL_sf, color="orange",size=0.005)+
        geom_sf(data=c3b_ternBFAL_sf, color="firebrick",size=0.005)+
        geom_sf(data=CRW_mb95c, color="orange",fill="orange",alpha=0.5)+
        geom_sf(data=CRW_tb95c, color="firebrick",fill="firebrick",alpha=0.5)+
      #geom_sf(data=npac_base_res, fill="grey60")
        coord_sf(xlim = c(-9e06, 9e06), ylim = c(-5.5e06, 1.25e07),expand=F)
      midBFALternBFALCRWrandom
      
      save(midBFALternBFALCRWrandom,file="3b.Rdata")
      load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/3b.Rdata")
      midBFALternBFALCRWrandom
      
# Tern LAAL x Tern BFAL
  # Tern LAAL
      c3c_ternLAAL_points <- CRW_sim_allpoints %>% filter(grepl('LAAL_TERN', Animal_ID))
      c3c_ternLAAL_points <- c3c_ternLAAL_points %>% filter(iter==1) 
      c3c_ternLAAL_points <- c3c_ternLAAL_points[,c(2,3,4)]
      
      c3c_ternLAAL_sf = st_as_sf(c3c_ternLAAL_points, coords = c("Longitude","Latitude"), remove = FALSE)
      
      sp::coordinates(c3c_ternLAAL_points) <- c("Longitude", "Latitude")
      grid_input <- calculate_sp_obj_extent(c3c_ternLAAL_points,0.1)
      
      c3c_ternLAAL_ud <-  kernelUD(c3c_ternLAAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      mid_iter_holder <- numeric(length = 3304)
      for (d in 1:length(c3c_ternLAAL_ud)) {
        c3c_ternLAAL_ud[[d]]@data$ud[is.na(c3c_ternLAAL_ud[[d]]@data$ud)] <- 0
        add <- c3c_ternLAAL_ud[[d]]@data$ud
        mid_iter_holder <- mid_iter_holder+add
      }
      mid_iter_holder <- mid_iter_holder/length(c3c_ternLAAL_ud)
      mid_iter_holder
      
      ##### modify existing estUD object with averaged values, then rename
      mid_iter_avg <- c3c_ternLAAL_ud[[1]]
      mid_iter_avg@data$ud <- mid_iter_holder
      mid_iter_avg@data$ud[is.na(mid_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
      image(mid_iter_avg)
      CRW_tl95c <- getverticeshr(mid_iter_avg,percent = 95)
      CRW_tl95c  <- st_as_sf(CRW_tl95c ) # 95th UD Contour
      # CRW_ml95c  <- CRW_ml95c %>% 
      #   st_set_crs(3832) %>%
      #   st_transform(crs = 4326) %>% # transform to WGS84.
      #   st_wrap_dateline() %>% # wrap around the dateline
      #   st_shift_longitude() %>%  
      #   st_union(by_feature = TRUE) %>% 
      #   st_transform(crs = 3832)
      
  # Tern BFAL
      c3c_ternBFAL_points <- CRW_sim_allpoints %>% filter(grepl('BFAL_TERN', Animal_ID))
      c3c_ternBFAL_points <- c3c_ternBFAL_points %>% filter(iter==1) 
      c3c_ternBFAL_points <- c3c_ternBFAL_points[,c(2,3,4)]
      
      c3c_ternBFAL_sf = st_as_sf(c3c_ternBFAL_points, coords = c("Longitude","Latitude"), remove = FALSE)
      
      
      sp::coordinates(c3c_ternBFAL_points) <- c("Longitude", "Latitude")
      grid_input <- calculate_sp_obj_extent(c3c_ternBFAL_points,0.1)
      
      c3c_ternBFAL_ud <-  kernelUD(c3c_ternBFAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      tern_iter_holder <- numeric(length = 3240)
      for (d in 1:length(c3c_ternBFAL_ud)) {
        c3c_ternBFAL_ud[[d]]@data$ud[is.na(c3c_ternBFAL_ud[[d]]@data$ud)] <- 0
        add <- c3c_ternBFAL_ud[[d]]@data$ud
        tern_iter_holder <- tern_iter_holder+add
      }
      tern_iter_holder <- tern_iter_holder/length(c3c_ternBFAL_ud)
      tern_iter_holder
      
      ##### modify existing estUD object with averaged values, then rename
      tern_iter_avg <- c3c_ternBFAL_ud[[1]]
      tern_iter_avg@data$ud <- tern_iter_holder
      tern_iter_avg@data$ud[is.na(tern_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
      image(tern_iter_avg)
      CRW_tb95c <- getverticeshr(tern_iter_avg,percent = 95)
      CRW_tb95c  <- st_as_sf(CRW_tb95c ) # 95th UD Contour
      # CRW_tl95c  <- CRW_tl95c %>% 
      #   st_set_crs(3832) %>%
      #   st_transform(crs = 4326) %>% # transform to WGS84.
      #   st_wrap_dateline() %>% # wrap around the dateline
      #   st_shift_longitude() %>%  
      #   st_union(by_feature = TRUE) %>% 
      #   st_transform(crs = 3832)
      
      # plot c3c
      ternLAALternBFALCRWrandom<- ggplot() + 
        # base map and other parameters
        geom_sf(data=c3c_ternLAAL_sf, color="orange",size=0.005)+
        geom_sf(data=c3c_ternBFAL_sf, color="firebrick",size=0.005)+
        geom_sf(data=CRW_tl95c, color="orange",fill="orange",alpha=0.5)+
        geom_sf(data=CRW_tb95c, color="firebrick",fill="firebrick",alpha=0.5)+
      #geom_sf(data=npac_base_res, fill="grey60")
        coord_sf(xlim = c(-9e06, 9e06), ylim = c(-5.5e06, 1.25e07),expand=F)
      ternLAALternBFALCRWrandom
      
      save(ternLAALternBFALCRWrandom,file="3c.Rdata")
      load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/3c.Rdata")
      ternLAALternBFALCRWrandom
      

# Plot everything ---------------------------------------------------------

      library(ggpubr)
      combined_figure <- ggarrange(midLAALternLAAL, midLAALternLAALtrackpermrandom,midLAALternLAALCRWrandom,
                                   midBFALternBFAL, midBFALternBFALtrackpermrandom,midBFALternBFALCRWrandom,
                                   ternLAALternBFAL,ternLAALternBFALtrackpermrandom,ternLAALternBFALCRWrandom,
                                   labels = c("a", "","" ,"b","","", "c","", ""), hjust=-5,
                                   ncol = 3, nrow = 3)
      combined_figure
      
