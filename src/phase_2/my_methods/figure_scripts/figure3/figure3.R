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

# CRITICAL: 
sf::sf_use_s2(FALSE)
# The above line is necessary as of April 17 2021, when 'sf' switched to using S2 spheroid 
# geometry instead of a GEOS routine that assumed projected coordinates. Until 'mapping' 
# updates its back-end, switching S2 geometry off is necessary for my script to work. 
# See: https://github.com/r-spatial/sf/issues/1649
# Great reference here: https://r-spatial.github.io/sf/articles/sf7.html

lcea <- "+proj=cea +lat_0=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

# import points
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

# import contours 
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/individual/midLAAL/master_script_contours/")
files <- dir("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/individual/", recursive=TRUE, full.names=TRUE, pattern="\\.Rdata$")
result <- lapply(files, load,envir = .GlobalEnv)

# Loading in low res base map since final figure gonna be huge
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/contours_rasters_figureData/basemap/")
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

library(ggpubr)
combined_figure <- ggarrange(midLAALternLAAL, midBFALternBFAL, ternLAALternBFAL, 
                             labels = c("a", "b","c"),
                             ncol = 1, nrow = 3)
combined_figure


# Column 2 - track permutation random examples ----------------------------
# you are getting randomized 95th contours from the permutation_XX_XX.R scripts - you are running the
# randomization procedure like 20 times and just pulling one of the xxx_iter_avg. You run it 
# in the other script to get it in the global environment, then getverticeshr here. 

# First, code for randomized tracks to get points and contours: 
# Midway LAAL x Tern LAAL
# Randomized points: 
add_lm <-lm
add_lt <-lt
add_lm$island <- "LAALmid"
add_lt$island <- "LAALtern"
resample_this <- rbind(add_lm,add_lt)
resample_this_island <- split(resample_this, resample_this$island)
resample_this_mid_track <- split(resample_this_island[[1]],resample_this_island[[1]]$track)
resample_this_tern_track <- split(resample_this_island[[2]],resample_this_island[[2]]$track)
resample_all_tracks <- c(resample_this_mid_track,resample_this_tern_track)
available_mid_tracks <- unique(resample_this_island[[1]]$track)
available_tern_tracks <- unique(resample_this_island[[2]]$track)
n_amt <- length(available_mid_tracks) # you are calculating these to preserve sample size in randomizations
n_att <- length(available_tern_tracks)
n_all <- n_amt+n_att
counts <- resample_this %>% count(id,track)


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
  #geom_sf(data=lm_points, color="orange",size=0.005)+
  #geom_sf(data=lt_points, color="firebrick",size=0.005)+
  geom_sf(data=ml95cr, color="orange",fill="orange",alpha=0.5)+
  geom_sf(data=tl95cr, color="firebrick",fill="firebrick",alpha=0.5)+
  geom_sf(data=npac_base_res, fill="grey60")+
  coord_sf(xlim = c(-2000000, 10000000), ylim = c(1464000, 13000000),expand=F)
midLAALternLAALtrackpermrandom
