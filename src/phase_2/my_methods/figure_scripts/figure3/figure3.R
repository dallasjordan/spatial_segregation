# Mapping 9-panel plot showing randomizations and original data. 
# Three cases: Midway LAAL x Tern LAAL - CRW shows strong overlap
#              Midway BFAL x Tern BFAL - CRW shows medium overlap 
#              Tern LAAL x Tern BFAL - CRW agrees with track perm, sig. segregation 
# Plot original data (points), plot one random of track perm, plot one random of CRW, each of these 
# have contours plotted over randomization points
# Ostensibly will be figure 3 (ended up being figure 4)
# NOTE - You did not end up making a 9 panel figure!
# You instead just used this script to make MidwayLAAL x TernLAAL comparisons for observed data,
# track ID switch contours, and CRW generated contours to show that the CRW contours made more sense.
# Dallas Jordan
# May 2 2022
  # Changing panel C to include the new 1000 iteration random walks - I have to redo this panel 
  # because the CRWs now are drawn from a distribution of distances and turning angles, so they 
  # look different

# May 17 2022 updates
  # Just added icons for Tern and Midway colony sites and changed x/y axis spacing/labels to 20 degrees

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

midway <- data.frame(longitude = -177.3761,
                    latitude = 28.2101)
tern <- data.frame(longitude = -166.284,
                   latitude = 23.870)

midway <- st_as_sf(midway, coords = c("longitude", "latitude"), 
                   crs = 4326, agr = "constant")
tern <- st_as_sf(tern, coords = c("longitude", "latitude"), 
                   crs = 4326, agr = "constant")

# CRITICAL: 
sf::sf_use_s2(FALSE)
# The above line is necessary as of April 17 2021, when 'sf' switched to using S2 spheroid 
# geometry instead of a GEOS routine that assumed projected coordinates. Until 'mapping' 
# updates its back-end, switching S2 geometry off is necessary for my script to work. 
# See: https://github.com/r-spatial/sf/issues/1649
# Great reference here: https://r-spatial.github.io/sf/articles/sf7.html

lcea <- "+proj=cea +lat_0=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
aeqd = "+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
eqc <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
pdc_mercator_proj<-"+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
load("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3/figure3_old/basemap_figure3_v2")

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
# One example left for reference of the old way I used to transform these:
# ml95c <- ml95c %>%
#   st_transform(crs = 4326) %>% # transform to WGS84.
#   st_wrap_dateline() %>% # wrap around the dateline
#   st_shift_longitude() %>%
#   st_union(by_feature = TRUE) %>%
#   st_transform(crs = 3832)

ml95c <- ml95c %>%
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

ml50c <- st_as_sf(vert50_midLAAL) # 50th UD Contour
ml50c <- ml50c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

tl95c <- st_as_sf(vert95_ternLAAL) # 95th UD Contour
tl95c <- tl95c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

tl50c <- st_as_sf(vert50_ternLAAL) # 50th UD Contour
tl50c <- tl50c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

mb95c <- st_as_sf(vert95_midBFAL) # 95th UD Contour
mb95c <- mb95c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

mb50c <- st_as_sf(vert50_midBFAL) # 50th UD Contour
mb50c <- mb50c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

tb95c <- st_as_sf(vert95_ternBFAL) # 95th UD Contour
tb95c <- tb95c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

tb50c <- st_as_sf(vert50_ternBFAL) # 50th UD Contour
tb50c <- tb50c %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

# Column 1 - original points with contours overlapping  -----------------------------------------------
# Midway LAAL x Tern LAAL 
lm_points = st_as_sf(lm, coords = c("x","y"), remove = FALSE, crs=4326)
lm_points <- lm_points %>%
  st_transform(eqc)

lt_points = st_as_sf(lt, coords = c("x","y"), remove = FALSE, crs=4326)
lt_points <- lt_points %>%
  st_transform(eqc)

midLAALternLAAL <- ggplot() + 
  # base map and other parameters
  #geom_sf(data=lm_points, color="orange",size=0.005)+
  #geom_sf(data=lt_points, color="firebrick",size=0.005)+
  geom_sf(data=ml95c, color="chocolate2",fill="chocolate2",alpha=0.5,size=1.25)+
  geom_sf(data=tl95c, color="darkred",fill="darkred",alpha=0.5,size=1.25)+
  #geom_sf(data=CRW_ml95c, color=NA,alpha=0)+
  geom_sf(data=npac_base_res, fill="grey60")+
  theme_bw()+
  scale_y_continuous(breaks = c(10, 30, 50, 70))+
  scale_x_continuous(breaks = c(-120,-140,-160,180,160,140,120))+
  geom_sf(data=midway, size=3,shape=17,fill="black",color="black")+
  geom_sf(data=tern, size=4,shape=18,fill="black",color="black")+
  coord_sf(expand=F)
  #coord_sf(xlim = c(-10246822, 4164347), ylim = c(-6563332, 15190864),expand=F)
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
############ now run other script to get xxx_iter_avg, can go to Randomized contours line after ############

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
ml95cr <- st_as_sf(ml95cr) # 95th UD Contour
ml95cr <- ml95cr %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

tl95cr <- getverticeshr(tern_iter_avg)
tl95cr <- st_as_sf(tl95cr) # 95th UD Contour
tl95cr <- tl95cr %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)

midLAALternLAALtrackpermrandom <- ggplot() + 
  # base map and other parameters
  #geom_sf(data=mid_iter_points, color="orange",size=0.005)+
  #geom_sf(data=tern_iter_points, color="firebrick",size=0.005)+
  geom_sf(data=ml95cr, color="chocolate2",fill="chocolate2",alpha=0.5,size=1.25)+
  geom_sf(data=tl95cr, color="darkred",fill="darkred",alpha=0.5,size=1.25)+
  #geom_sf(data=CRW_ml95c, color=NA,alpha=0)+
  geom_sf(data=npac_base_res, fill="grey60")+
  theme_bw()+
  geom_sf(data=midway, size=3,shape=17,fill="black",color="black")+
  geom_sf(data=tern, size=4,shape=18,fill="black",color="black")+
  scale_y_continuous(breaks = c(10, 30, 50, 70))+
  scale_x_continuous(breaks = c(-120,-140,-160,180,160,140,120))+
  coord_sf(expand=F)
  #coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 13000000),expand=F)
midLAALternLAALtrackpermrandom

# save this so you don't save to re-generate from scripts
save(midLAALternLAALtrackpermrandom,file="2a_v4.Rdata")
load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/2a_v4.Rdata")
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3/figure3_old/2a.Rdata")
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

# May 2 2022 update: new CRW file of CRW drawn from distributions
load(file="/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/CRW_v4_all.Rdata")
spp_col_list_all <- save_spp_col_list
list_to_df_add_iteration_number <- function(spp_col_list){
  iteration_vector <- 1:1000
  num <- length(spp_col_list)/1000
  if (num==3){
    spp_col_one <- spp_col_list[1:1000]
    spp_col_two <- spp_col_list[1001:2000]
    spp_col_three <- spp_col_list[2001:3000]
    for (a in 1:1000){
      spp_col_one[[a]]$iteration_number <- a
    }
    for (b in 1:1000){
      spp_col_two[[b]]$iteration_number <- b
    }
    for (c in 1:1000){
      spp_col_three[[c]]$iteration_number <- c
    }
    merge_list <- c(spp_col_one, spp_col_two, spp_col_three)
  }
  if (num==4){
    spp_col_one <- spp_col_list[1:1000]
    spp_col_two <- spp_col_list[1001:2000]
    spp_col_three <- spp_col_list[2001:3000]
    spp_col_four<- spp_col_list[3001:4000]
    for (a in 1:1000){
      spp_col_one[[a]]$iteration_number <- a
    }
    for (b in 1:1000){
      spp_col_two[[b]]$iteration_number <- b
    }
    for (c in 1:1000){
      spp_col_three[[c]]$iteration_number <- c
    }
    for (d in 1:1000){
      spp_col_four[[d]]$iteration_number <- d
    }
    merge_list <- c(spp_col_one, spp_col_two, spp_col_three, spp_col_four)
  }
  output <- bind_rows(merge_list)
  return(output)
}

# convert to large dataframes
CRW_sim_allpoints <- list_to_df_add_iteration_number(spp_col_list_all)
#CRW_sim_allpoints <- read.csv("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/one_iteration_of_CRW_pts.csv")
#CRW_sim_allpoints <- read.csv("G:/Other computers/My MacBook Pro/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/julia_simulations/CRW_permutation_results/real_data/CRW_pts.csv")
#CRW_sim_allpoints <- read.csv("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/julia_simulations/CRW_permutation_results/real_data/CRW_pts.csv")

# Midway LAAL x Tern LAAL
  # Midway LAAL
      c3a_midLAAL_points <- CRW_sim_allpoints %>% filter(grepl('LAAL_MID', Animal_ID))
      c3a_midLAAL_points <- c3a_midLAAL_points %>% filter(iteration_number==1)
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
      # max_lat <- function(x){
      #   if_else(x>6339452,6339452,x) 
      # }
      # c3a_midLAAL_points$Latitude <- max_lat(c3a_midLAAL_points$Latitude)
      # ind <- with(c3a_midLAAL_points, (Latitude == 6339452))
      # c3a_midLAAL_points <- c3a_midLAAL_points[!ind,]
      
      #convert to sf for easy ggplot2 plotting
      c3a_midLAAL_points_sf <- st_as_sf(c3a_midLAAL_points, coords = c("Longitude","Latitude"), remove = TRUE, na.fail = TRUE, crs=pdc_mercator_proj)
      # c3a_midLAAL_points_sf <- c3a_midLAAL_points_sf %>%
      #   st_transform(3349)
      # 
      # g=st_sfc(st_point(c(-425223.3,6303822)),crs=lcea)
      # g <- g %>% st_transform(crs=3349)
      # 
      # g=st_sfc(st_point(c(-5843408.58,6363885)),crs=lcea)
      # g <- g %>% st_transform(crs=4326)
      # g
      # 
      # g=st_sfc(st_point(c(180,85)),crs=4326)
      # g <- g%>% st_transform(crs=lcea)
      # g

      
      
      c3a_midLAAL_points <- c3a_midLAAL_points[,c(1,2,3)]
      sp::coordinates(c3a_midLAAL_points) <- c("Longitude", "Latitude")
      #proj4string(c3a_midLAAL_points) <- CRS(lcea) # placeholder
      #c3a_midLAAL_points <- spTransform(c3a_midLAAL_points,CRS("+init=epsg:3349"))
      grid_input <- calculate_sp_obj_extent(c3a_midLAAL_points,0.1)
      
      c3a_midLAAL_ud <-  kernelUD(c3a_midLAAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      mid_iter_holder <- numeric(length = 1472)
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
      st_crs(CRW_ml95c) <- 3349
      st_is_valid(CRW_ml95c)
      CRW_ml95c <- CRW_ml95c %>%
        st_transform(4326) %>%
        st_shift_longitude() %>%
        st_transform(eqc)
      
      
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
      c3a_ternLAAL_points <- c3a_ternLAAL_points %>% filter(iteration_number==1)
    
      # c3a_ternLAAL_points$Latitude <- max_lat(c3a_ternLAAL_points$Latitude)
      # ind <- with(c3a_ternLAAL_points, (Latitude == 6339452))
      # c3a_ternLAAL_points <- c3a_ternLAAL_points[!ind,]
      
      #convert to sf for easy ggplot2 plotting
      c3a_ternLAAL_points_sf <- st_as_sf(c3a_ternLAAL_points, coords = c("Longitude","Latitude"), remove = TRUE, na.fail = TRUE, crs=pdc_mercator_proj)
      # c3a_ternLAAL_points_sf <- c3a_ternLAAL_points_sf %>%
      #   st_transform(3349)
      
      c3a_ternLAAL_points <- c3a_ternLAAL_points[,c(1,2,3)]
      sp::coordinates(c3a_ternLAAL_points) <- c("Longitude", "Latitude")
      # proj4string(c3a_ternLAAL_points) <- CRS(lcea) # placeholder
      # c3a_ternLAAL_points <- spTransform(c3a_ternLAAL_points,CRS("+init=epsg:3349"))
      grid_input <- calculate_sp_obj_extent(c3a_ternLAAL_points,0.1)
      
      c3a_ternLAAL_ud <-  kernelUD(c3a_ternLAAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      tern_iter_holder <- numeric(length = 1376)
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
      st_crs(CRW_tl95c) <- 3349
      st_is_valid(CRW_tl95c)
      
      CRW_tl95c  <- CRW_tl95c %>%
        st_transform(4326) %>%
        st_shift_longitude() %>%
        st_transform(eqc)
      
      # 
      # plot c3a
      # npac_base_i <- ptolemy::extract_gshhg(CRW_ml95c, resolution = "c", epsg = NULL, buffer = 5000,
      #                                       simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS
      # save(npac_base_i,file="basemap_figure3_v2")
      # load("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3/figure3_old/basemap_figure3_v2")
      midLAALternLAALCRWrandom <- ggplot() + 
        # base map and other parameters
        #geom_sf(data=c3a_midLAAL_points_sf, color="orange",size=0.005)+
        #geom_sf(data=c3a_ternLAAL_points_sf, color="firebrick",size=0.005)+
        geom_sf(data=CRW_ml95c, color="chocolate2",fill="chocolate2",alpha=0.5, size=1.25)+
        geom_sf(data=CRW_tl95c, color="darkred",fill="darkred",alpha=0.5, size=1.25)+
        #geom_sf(data=land_mask, fill="grey60")+
        geom_sf(data=npac_base_res,fill="grey60")+
        theme_bw()+
        geom_sf(data=midway, size=3,shape=17,fill="black",color="black")+
        geom_sf(data=tern, size=4,shape=18,fill="black",color="black")+
        scale_y_continuous(breaks = c(10, 30, 50, 70))+
        scale_x_continuous(breaks = c(-120,-140,-160,180,160,140,120))+
        coord_sf(expand=F)
        # coord_sf(expand=F)
      midLAALternLAALCRWrandom
      
  setwd("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3")
  save(midLAALternLAALCRWrandom,file="3a_v5.Rdata")
      load("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3/3a_v3.Rdata")
      load("E:/spatial_segregation/src/phase_2/my_methods/figure_scripts/figure3/3a.Rdata")
      midLAALternLAALCRWrandom
      
      
# Midway BFAL x Tern BFAL
    # Midway BFAL
      c3b_midBFAL_points <- CRW_sim_allpoints %>% filter(grepl('BFAL_MID', Animal_ID))
      c3b_midBFAL_points <- c3b_midBFAL_points %>% filter(iter==1) 
      c3b_midBFAL_points <- c3b_midBFAL_points[,c(2,3,4)]
      
      c3b_midBFAL_sf = st_as_sf(c3b_midBFAL_points, coords = c("Longitude","Latitude"), remove = TRUE, na.fail = TRUE, crs=pdc_mercator_proj)
      
      sp::coordinates(c3b_midBFAL_points) <- c("Longitude", "Latitude")
      grid_input <- calculate_sp_obj_extent(c3b_midBFAL_points,0.1)
      
      c3b_midBFAL_ud <-  kernelUD(c3b_midBFAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      mid_iter_holder <- numeric(length = 1334)
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
      st_crs(CRW_mb95c) <- 3349
      st_is_valid(CRW_mb95c)
      
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
      
      c3b_ternBFAL_sf = st_as_sf(c3b_ternBFAL_points, coords = c("Longitude","Latitude"), remove = TRUE, na.fail = TRUE, crs=pdc_mercator_proj)
      
      
      sp::coordinates(c3b_ternBFAL_points) <- c("Longitude", "Latitude")
      grid_input <- calculate_sp_obj_extent(c3b_ternBFAL_points,0.1)
      
      c3b_ternBFAL_ud <-  kernelUD(c3b_ternBFAL_points, grid=grid_input,same4all=T,extent=0.1,h=150000)
      
      # average into 1 UD
      tern_iter_holder <- numeric(length = 1551)
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
      st_crs(CRW_tb95c) <- 3349
      st_is_valid(CRW_tb95c)
      
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
        #geom_sf(data=c3b_midBFAL_sf, color="orange",size=0.005)+
        #geom_sf(data=c3b_ternBFAL_sf, color="firebrick",size=0.005)+
        geom_sf(data=CRW_mb95c, color="orange",fill="orange",alpha=0.5)+
        geom_sf(data=CRW_tb95c, color="firebrick",fill="firebrick",alpha=0.5)+
        geom_sf(data=land_mask, fill="grey60")+
        coord_sf(expand=F)
      midBFALternBFALCRWrandom
      
      setwd("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3")
      save(midBFALternBFALCRWrandom,file="3b_v3.Rdata")
        load("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3/3b_v3.Rdata")
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
      mid_iter_holder <- numeric(length = 1551)
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
      st_crs(CRW_tl95c) <- 3349
      st_is_valid(CRW_tl95c)
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
      tern_iter_holder <- numeric(length = 1551)
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
      st_crs(CRW_tb95c) <- 3349
      st_is_valid(CRW_tb95c)
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
        #geom_sf(data=c3c_ternLAAL_sf, color="orange",size=0.005)+
        #geom_sf(data=c3c_ternBFAL_sf, color="firebrick",size=0.005)+
        geom_sf(data=CRW_tl95c, color="orange",fill="orange",alpha=0.5)+
        geom_sf(data=CRW_tb95c, color="firebrick",fill="firebrick",alpha=0.5)+
        geom_sf(data=land_mask, fill="grey60")+
        coord_sf(expand=F)
      ternLAALternBFALCRWrandom
      
      setwd("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3")
      save(ternLAALternBFALCRWrandom,file="3c_v3.Rdata")
        load("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/src/phase_2/my_methods/figure_scripts/figure3/3c_v2")
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
      
