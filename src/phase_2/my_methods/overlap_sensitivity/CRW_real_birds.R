# This script is CRW code for simulations for REAL BIRD TRACKS. This script was created at a 'stopping point'
# of the sim_loops_CRW_BIRDS script, where the data was left in LCEA projection, with no changes to the simulation
# portion of the code. Here, I change the data to PDC mercator and implement the land mask and lat/lon filters to 
# get more reasonable simulations of birds
# Jan 29 2022
# Dallas Jordan

library(geosphere)
library(beepr)
library(lubridate)
library(maptools)
library(reproj)
library(dplyr)
library(tidyr)
library(ggplot2)
library(adehabitatLT)
library(adehabitatHR)
library(sf)
library(reshape2)

## Read in Spp_A and Spp_B from "simulation_exploration.R"
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/import_for_CRW/BFAL/")
d_spps_bfal_tern <- list.files(pattern = ".csv$", recursive=T)
d_spps_bfal_tern <- lapply(d_spps_bfal_tern, function(x) read.csv(x, header=T))
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/import_for_CRW/LAAL/")
d_spps_laal_tern <- list.files(pattern = ".csv$", recursive=T)
d_spps_laal_tern <- lapply(d_spps_laal_tern, function(x) read.csv(x, header=T))

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/BFAL/")
d_spps_bfal_mid <- list.files(pattern = ".csv$", recursive=T)
d_spps_bfal_mid <- lapply(d_spps_bfal_mid, function(x) read.csv(x, header=T))
setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/LAAL/")
d_spps_laal_mid <- list.files(pattern = ".csv$", recursive=T)
d_spps_laal_mid <- lapply(d_spps_laal_mid, function(x) read.csv(x, header=T))

#BFAL TERN
d_spps_bfal_ternID <- c()
for(i in 1:38){
  w <- d_spps_bfal_tern[[i]]
  w$ID <- paste("BFAL_TERN_", i, sep="")
  if(is.na(ymd(w$GMT_YYYY.MM.DD)) == TRUE){
    w$GMT <- mdy_hm(paste(w$GMT_YYYY.MM.DD, w$GMT_HH.mm.ss, sep=" "))
  }else{
    w$GMT <- ymd_hm(paste(w$GMT_YYYY.MM.DD, w$GMT_HH.mm.ss, sep=" "))
  }
  d_spps_bfal_ternID <- rbind(d_spps_bfal_ternID, w)
}
d_spps_bfal_ternID <- d_spps_bfal_ternID %>% dplyr::select(ID, GMT, Longitude, Latitude )

#LAAL TERN
d_spps_laal_ternID <- c()
for(i in 1:45){
  w <- d_spps_laal_tern[[i]]
  w$ID <- paste("LAAL_TERN_", i, sep="")
  if(is.na(ymd(w$GMT_YYYY.MM.DD)) == TRUE){
    w$GMT <- mdy_hm(paste(w$GMT_YYYY.MM.DD, w$GMT_HH.mm.ss, sep=" "))
  }else{
    w$GMT <- ymd_hm(paste(w$GMT_YYYY.MM.DD, w$GMT_HH.mm.ss, sep=" "))
  }
  d_spps_laal_ternID <- rbind(d_spps_laal_ternID, w)
}
d_spps_laal_ternID <- d_spps_laal_ternID %>% dplyr::select(ID, GMT, zm.lon, zm.lat )
names(d_spps_laal_ternID) <- names(d_spps_bfal_ternID)

#BFAL MID
d_spps_bfal_midID <- c()
for(i in 1:23){
  w <- d_spps_bfal_mid[[i]]
  w$ID <- paste("BFAL_MID_", i, sep="")
  w$dtime <- ymd_hms(w$dtime)
  ee <- w %>% dplyr::select(ID, dtime, Longitude, Latitude)
  colnames(ee) <- colnames(d_spps_bfal_ternID)
  d_spps_bfal_midID <- rbind(d_spps_bfal_midID, ee)
}

#LAAL MID
d_spps_laal_midID <- c()
for(i in 1:37){
  w <- d_spps_laal_mid[[i]]
  w$ID <- paste("LAAL_MID_", i, sep="")
  w$dtime <- ymd_hms(w$dtime)
  ee <- w %>% dplyr::select(ID, dtime, Longitude, Latitude)
  colnames(ee) <- colnames(d_spps_bfal_ternID)
  d_spps_laal_midID <- rbind(d_spps_laal_midID, ee)
}

d_birds <- data.frame(rbind(d_spps_bfal_midID, d_spps_bfal_ternID, d_spps_laal_midID, d_spps_laal_ternID))
newLL <- reproj(cbind(d_birds$Longitude, d_birds$Latitude), target = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", source = "+proj=longlat +datum=WGS84 +no_defs")
d_birds$Longitude2 <- newLL[,1]
d_birds$Latitude2 <- newLL[,2]
d_birds$match <- paste(d_birds$ID, d_birds$GMT, sep=".")

d_birds_traj <- as.ltraj(cbind(d_birds$Longitude2, d_birds$Latitude2), date=d_birds$GMT, id = d_birds$ID, typeII = T,
                         proj4string = CRS("+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
d_birds_sp <- ltraj2spdf(d_birds_traj)
d_birds_df <- data.frame(d_birds_sp)
d_birds_df <- merge(d_birds_df, d_birds[,5:7], by.x="pkey", by.y="match")
colnames(d_birds_df)[13] <- "Longitude"
#d_birds_df <- data.frame(cbind(d_birds_df, d_birds$Latitude2))
colnames(d_birds_df)[14] <- "Latitude"
nms <- data.frame(do.call(rbind, strsplit(d_birds_df$pkey, split="_", fixed=T)))
nms2 <- data.frame(do.call(rbind, strsplit(nms$X3, split=".", fixed=T)))
nms <- data.frame(nms, nms2)
colnames(nms) <- c("SPP", "tag_loc", "ptID", "track", "trackpt")

d_birds_df <- data.frame(d_birds_df, nms)
d_birds_df$Animal_ID <- paste(d_birds_df$SPP, d_birds_df$tag_loc, d_birds_df$track, sep="_")

angles_LAAL <- d_birds_df$abs.angle[d_birds_df$SPP == "LAAL"]
angles_LAAL <- angles_LAAL[!is.na(angles_LAAL)]
angles_LAAL <- (angles_LAAL * 180) / pi # convert to degrees from radians

angles_BFAL <- d_birds_df$abs.angle[d_birds_df$SPP == "BFAL"]
angles_BFAL <- angles_BFAL[!is.na(angles_BFAL)]
angles_BFAL <- (angles_BFAL * 180) / pi

dist_LAAL <- d_birds_df$dist[d_birds_df$SPP == "LAAL"]
dist_LAAL <- dist_LAAL[!is.na(dist_LAAL)]
dist_BFAL <- d_birds_df$dist[d_birds_df$SPP == "BFAL"]
dist_BFAL <- dist_BFAL[!is.na(dist_BFAL)]

comb <- d_birds_df

### pull first row of each, simulate from there
first_pts <- NULL
for(i in unique(comb$Animal_ID)){
  p <- comb[comb$Animal_ID == i,]
  p <- p[which(p$trackpt == min(p$trackpt)),]
  first_pts <- rbind(first_pts, p)
}

### Do each species separately ----
rand_out_pts_birds <- NULL # save both species points
rand_out_plot_birds <- NULL # save ggplot
rand_out_blobs_birds <- NULL # save udois
rand_out_table_birds <- NULL # save udoi, ba, other value
rand_out_plot_all <- NULL
rand_out_plot_mid <- NULL
rand_out_plot_tern <- NULL

rand_out_poly_all <- NULL
rand_out_poly_mid <- NULL
rand_out_poly_tern <- NULL

setwd("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity")

##### simulations #####

library(sf)
sf::sf_use_s2(FALSE)
library(dplyr)
land_mask_bbox<-st_bbox(c(xmin = -3008182, xmax = 10844089, ymax = 18027535, ymin = 185132), crs = st_crs(3832))
land_mask_sfc <- st_as_sfc(land_mask_bbox)
land_mask <- ptolemy::extract_gshhg(land_mask_sfc, resolution = "i", epsg = NULL, buffer = 5000,
                                    simplify = FALSE) 
land_mask$type = "land"

simtype <- "birds"
big_storage <- NULL
# initialize empty list of length 100, then shove sim_CRW_pts into the list to eliminate the highest level 
# of the for loop
# make functions for intermediate steps

# just doing 1 iteration - you can bump this up if you want, but when you plot figure 3
# you just plot 1 iteration anyway
for(j in 34:34){
  sim_CRW_pts <- NULL
  store_here <- NULL
  #for(i in unique(first_pts$Animal_ID[first_pts$SPP == "LAAL"])){
  for(i in unique(first_pts$Animal_ID[first_pts$SPP == "LAAL"])){  
    a <- first_pts[first_pts$Animal_ID == i,]
    
    temp_pts <- a
    temp_pts$new_dis <- 0
    temp_pts$new_ang <- 0
    temp_pts$date <- ymd_hms(temp_pts$date)
    temp_pts$num <- 0
    temp_pts$iter <- 0
    
    store_here <- NULL
    
    lenn <- nrow(d_birds_df[d_birds_df$Animal_ID == i,])
    
    for(z in 1:1){
      bloop <- temp_pts %>% dplyr::select(Animal_ID, Longitude, Latitude,date, new_dis, new_ang, num)
      laal_sim_counter <- 0
      while (laal_sim_counter<(lenn)){
        new_dis <- sample(dist_LAAL,1)
        new_ang <- sample(angles_LAAL,1)
          new_x <- bloop$Longitude[laal_sim_counter+1] + new_dis*cos(new_ang)
          new_y <- bloop$Latitude[laal_sim_counter+1] + new_dis*sin(new_ang)
          timee <- ymd_hms(bloop$date[laal_sim_counter+1]) + (laal_sim_counter+1)
          
          new_point <- st_point(x = c(new_x, new_y), dim = "XY") %>% 
            st_sfc(crs=st_crs(land_mask)) %>% 
            st_sf()
          in_land <- st_join(new_point, land_mask, join = st_intersects)
          if (is.na(in_land$type) & (-1708182<new_x) & (new_x<9844089) & (1185132<new_y) & (new_y<9310913)){
            e <- data.frame("Animal_ID" = i, "Longitude" = new_x, "Latitude" = new_y, 
                            "date" = timee, "new_dis" = new_dis, "new_ang" = new_ang, 
                            "num" = laal_sim_counter+1)
            bloop <- rbind(bloop, e)
            laal_sim_counter<- laal_sim_counter+1
            }
        }
        bloop$iter <- j
        store_here <- rbind(store_here, bloop)
        
    }
    #print(i)
    sim_CRW_pts <- rbind(sim_CRW_pts, store_here)
  }
  
  # xmin : -1708182
  # xmax: 9844089
  # ymin: 1185132
  # ymax: 9310913
  
  
  # test plotting
  # bloop_sf = st_as_sf(bloop, coords = c("Longitude", "Latitude"),
  #                  crs = 3349)
  # 
  # ggplot(data=land_mask)+
  #   geom_sf()+
  #   geom_sf(data=point)
    
  
  sim_CRW_pts2 <- NULL
  store_here <- NULL
  for(i in unique(first_pts$Animal_ID[first_pts$SPP == "BFAL"])){
    a <- first_pts[first_pts$Animal_ID == i,]
    
    temp_pts <- a
    temp_pts$new_dis <- 0
    temp_pts$new_ang <- 0
    temp_pts$date <- ymd_hms(temp_pts$date)
    temp_pts$num <- 0
    temp_pts$iter <- 0
    
    store_here <- NULL
    
    lenn <- nrow(d_birds_df[d_birds_df$Animal_ID == i,])
    
    for(z in 1:1){
      bloop <- temp_pts %>% dplyr::select(Animal_ID, Longitude, Latitude,date, new_dis, new_ang, num)
      bfal_sim_counter <-0
      while (bfal_sim_counter<(lenn)){
        new_dis <- sample(dist_BFAL,1)
        new_ang <- sample(angles_BFAL,1)
        new_x <- bloop$Longitude[bfal_sim_counter+1] + new_dis*cos(new_ang)
        new_y <- bloop$Latitude[bfal_sim_counter+1] + new_dis*sin(new_ang)
        timee <- ymd_hms(bloop$date[bfal_sim_counter+1]) + (bfal_sim_counter+1)
        
        new_point <- st_point(x = c(new_x, new_y), dim = "XY") %>% 
          st_sfc(crs=st_crs(land_mask)) %>% 
          st_sf()
        in_land <- st_join(new_point, land_mask, join = st_intersects)
        if (is.na(in_land$type) & (-1708182<new_x) & (new_x<9844089) & (1185132<new_y) & (new_y<9310913)){
          e <- data.frame("Animal_ID" = i, "Longitude" = new_x, "Latitude" = new_y, 
                          "date" = timee, "new_dis" = new_dis, "new_ang" = new_ang, 
                          "num" = bfal_sim_counter+1)
          bloop <- rbind(bloop, e)
          bfal_sim_counter<- bfal_sim_counter+1
        }
      }
      bloop$iter <- j
      store_here <- rbind(store_here, bloop)
    }
    #print(i)
    sim_CRW_pts2 <- rbind(sim_CRW_pts2, store_here)
  }
  sim_CRW_pts <- rbind(sim_CRW_pts, sim_CRW_pts2)
  big_storage <- rbind(big_storage, sim_CRW_pts)
  
  ### Assign new names to simulated points, convert to adehabitat points, combine with original points
  ### In a for loop create a new column with resampled track IDs (but still maintaining track IDs)
  ### Select all the real animals, calculate UDOI and SOI. Paste out to table
  
  sim_CRW_pts$Animal_ID <- paste( "sim_", sim_CRW_pts$Animal_ID, sep="")
  
  sim_CRW_sp <- sim_CRW_pts[!is.na(sim_CRW_pts$Longitude) & !is.na(sim_CRW_pts$Latitude),]
  sim_CRW_sp <- SpatialPointsDataFrame(coords = cbind(sim_CRW_pts$Longitude, sim_CRW_pts$Latitude), data = sim_CRW_pts)
  #sim_pts <- sim_CRW_pts[,c(1:4)]
  #sppA_pts <- spp1_shp[,c(1:4)]
  #sppB_pts <- spp2_shp[,c(1:4)]
  #comb_spp <- spRbind(sim_pts,sppA_pts)
  
  beep()
  
  # Dallas code starts here -------
  #all_data_1 <- sim_CRW_sp[sim_CRW_sp$Latitude < 90,1]
  all_data_1 <- sim_CRW_sp[,1]
  save(all_data_1, file="checkpoint.Rdata")
  # functions and objects
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
  pdc_mercator_proj<-"+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  # lcea <- "+proj=cea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" # what I used to use - is the 20, 180 necessary?
  lcea <- "+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs" 
  years <- c("2008","2009","2010","2011","2012")
  
  #sp::coordinates(all_data_1) <- c("x", "y")
  #proj4string(all_data_1) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # placeholder
  proj4string(all_data_1) <- CRS(pdc_mercator_proj)
  #all_data_1 <- spTransform(all_data_1,lcea)
  # Next line is critical for 300km resolution of KDEs
  grid_input <- calculate_sp_obj_extent(all_data_1,0.1) 
  
  # Generate individual KDE -------------------------------------------------
  
  # kernelUD output currently is probability density: absolute density (events per unit area) divided by total number of events
  results <- kernelUD(all_data_1, grid=grid_input, same4all=T, extent=0.1, h=150000)# grid-input guarantees AT LEAST 300km x 300km - in actuality, it is very slightly more than 300km x 300km.
  
  
  # Alter lines 87-90 with the results of "names" to see individual tracks. You might not want to do this because you have thousands of tracks
  names <- names(results)
  lm_kde <- results[grep("LAAL_MID",names)]
  lt_kde <- results[grep("LAAL_TERN",names)]
  bm_kde <- results[grep("BFAL_MID",names)]
  bt_kde <- results[grep("BFAL_TERN",names)]
  all_LAAL <- append(lm_kde, lt_kde)
  all_BFAL <- append(bm_kde, bt_kde)
  
  # Averaging -> new KDEs, account for track length ------------------------
  
  # "0" values are preserved here, not converted to NA
  # lt_holder length will change here if your cell size changes (calc_sp_obj_extent)
  len <- nrow(results$sim_LAAL_MID_1@data)
  
  # average lm #
  lm_holder <- numeric(length = len)
  for (i in 1:length(lm_kde)) {
    lm_kde[[i]]@data$ud[is.na(lm_kde[[i]]@data$ud)] <- 0
    add <- lm_kde[[i]]@data$ud
    lm_holder <- lm_holder+add
  }
  lm_holder <- lm_holder/length(lm_kde)
  
  
  ##### modify existing estUD object with averaged values, then rename
  lm_averaged_estUD <- lm_kde[[1]]
  lm_averaged_estUD@data$ud <- lm_holder
  lm_averaged_estUD@data$ud[is.na(lm_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
  
  # average lt #
  lt_holder <- numeric(length = len)
  for (i in 1:length(lt_kde)) {
    lt_kde[[i]]@data$ud[is.na(lt_kde[[i]]@data$ud)] <- 0
    add <- lt_kde[[i]]@data$ud
    lt_holder <- lt_holder+add
  }
  lt_holder <- lt_holder/length(lt_kde)
  
  
  ##### modify existing estUD object with averaged values, then rename
  lt_averaged_estUD <- lt_kde[[1]]
  lt_averaged_estUD@data$ud <- lt_holder
  lt_averaged_estUD@data$ud[is.na(lt_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
  
  # average bm #
  bm_holder <- numeric(length = len)
  for (i in 1:length(bm_kde)) {
    bm_kde[[i]]@data$ud[is.na(bm_kde[[i]]@data$ud)] <- 0
    add <- bm_kde[[i]]@data$ud
    bm_holder <- bm_holder+add
  }
  bm_holder <- bm_holder/length(bm_kde)
  
  ##### modify existing estUD object with averaged values, then rename
  bm_averaged_estUD <- bm_kde[[1]]
  bm_averaged_estUD@data$ud <- bm_holder
  bm_averaged_estUD@data$ud[is.na(bm_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
  
  # average bt #
  bt_holder <- numeric(length = len)
  for (i in 1:length(bt_kde)) {
    bt_kde[[i]]@data$ud[is.na(bt_kde[[i]]@data$ud)] <- 0
    add <- bt_kde[[i]]@data$ud
    bt_holder <- bt_holder+add
  }
  bt_holder <- bt_holder/length(bt_kde)
  
  ##### modify existing estUD object with averaged values, then rename
  bt_averaged_estUD <- bt_kde[[1]]
  bt_averaged_estUD@data$ud <- bt_holder
  bt_averaged_estUD@data$ud[is.na(bt_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
  
  # average all_LAAL #
  all_LAAL_holder <- numeric(length = len)
  for (i in 1:length(all_LAAL)) {
    all_LAAL[[i]]@data$ud[is.na(all_LAAL[[i]]@data$ud)] <- 0
    add <- all_LAAL[[i]]@data$ud
    all_LAAL_holder <- all_LAAL_holder+add
  }
  all_LAAL_holder <- all_LAAL_holder/length(all_LAAL)
  
  ##### modify existing estUD object with averaged values, then rename
  all_LAAL_averaged_estUD <- all_LAAL[[1]]
  all_LAAL_averaged_estUD@data$ud <- all_LAAL_holder
  all_LAAL_averaged_estUD@data$ud[is.na(all_LAAL_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
  
  # average all_BFAL #
  all_BFAL_holder <- numeric(length = 893)
  for (i in 1:length(all_BFAL)) {
    all_BFAL[[i]]@data$ud[is.na(all_BFAL[[i]]@data$ud)] <- 0
    add <- all_BFAL[[i]]@data$ud
    all_BFAL_holder <- all_BFAL_holder+add
  }
  all_BFAL_holder <- all_BFAL_holder/length(all_BFAL)
  
  ##### modify existing estUD object with averaged values, then rename
  all_BFAL_averaged_estUD <- all_BFAL[[1]]
  all_BFAL_averaged_estUD@data$ud <- all_BFAL_holder
  all_BFAL_averaged_estUD@data$ud[is.na(all_BFAL_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
  
  # Generate class KDE ------------------------------------------------------
  
  
  allLAAL <- all_LAAL_averaged_estUD
  allBFAL <- all_BFAL_averaged_estUD
  midLAAL <- lm_averaged_estUD
  midBFAL <- bm_averaged_estUD
  ternLAAL <- lt_averaged_estUD
  ternBFAL <- bt_averaged_estUD
  
  allLAAL_v_allBFAL <- list(allLAAL,allBFAL)
  class(allLAAL_v_allBFAL)<-"estUDm"
  names(allLAAL_v_allBFAL)<-c("allLAAL","allBFAL")
  
  # LAAL between islands
  midLAAL_v_ternLAAL <- list(midLAAL,ternLAAL)
  class(midLAAL_v_ternLAAL)<-"estUDm"
  names(midLAAL_v_ternLAAL)<-c("midLAAL","ternLAAL")
  
  # BFAL between islands
  midBFAL_v_ternBFAL <- list(midBFAL,ternBFAL)
  class(midBFAL_v_ternBFAL)<-"estUDm"
  names(midBFAL_v_ternBFAL)<-c("midBFAL","ternBFAL")
  
  # within Midway
  midLAAL_v_midBFAL <- list(midLAAL,midBFAL)
  class(midLAAL_v_midBFAL)<-"estUDm"
  names(midLAAL_v_midBFAL)<-c("midLAAL","midBFAL")
  
  # within Tern
  ternLAAL_v_ternBFAL <- list(ternLAAL,ternBFAL)
  class(ternLAAL_v_ternBFAL)<-"estUDm"
  names(ternLAAL_v_ternBFAL)<-c("ternLAAL","ternBFAL")
  
  # Overlap calculations ----------------------------------------------------
  # all LAAL vs all BFAL
  alab_UDOI_95 <- kerneloverlaphr(allLAAL_v_allBFAL, method="UDOI", percent=95, conditional=T)
  alab_UDOI_50 <- kerneloverlaphr(allLAAL_v_allBFAL, method="UDOI", percent=50, conditional=T)
  alab_PHR_95 <- kerneloverlaphr(allLAAL_v_allBFAL, method="PHR", percent=95, conditional=T)
  alab_PHR_50 <- kerneloverlaphr(allLAAL_v_allBFAL, method="PHR", percent=50, conditional=T)
  alab_95_BA <- kerneloverlaphr(allLAAL_v_allBFAL, method="BA", percent=95, conditional = T) 
  alab_50_BA <- kerneloverlaphr(allLAAL_v_allBFAL, method="BA", percent=50, conditional = T)
  alab_UDOI_95_test_stat <- alab_UDOI_95[1,2]
  alab_UDOI_50_test_stat <- alab_UDOI_50[1,2]
  alab_PHR_95_test_stat <- alab_PHR_95[1,2]
  alab_PHR_50_test_stat <- alab_PHR_50[1,2]
  alab_95_BA_test_stat <- alab_95_BA[1,2]
  alab_50_BA_test_stat <- alab_50_BA[1,2]
  
  # LAAL between islands
  mltl_UDOI_95 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="UDOI", percent=95, conditional=T)
  mltl_UDOI_50 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="UDOI", percent=50, conditional=T)
  mltl_PHR_95 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="PHR", percent=95, conditional=T)
  mltl_PHR_50 <- kerneloverlaphr(midLAAL_v_ternLAAL, method="PHR", percent=50, conditional=T)
  mltl_95_BA <- kerneloverlaphr(midLAAL_v_ternLAAL, method="BA", percent=95, conditional = T) 
  mltl_50_BA <- kerneloverlaphr(midLAAL_v_ternLAAL, method="BA", percent=50, conditional = T)
  mltl_UDOI_95_test_stat <- mltl_UDOI_95[1,2]
  mltl_UDOI_50_test_stat <- mltl_UDOI_50[1,2]
  mltl_PHR_95_test_stat <- mltl_PHR_95[1,2]
  mltl_PHR_50_test_stat <- mltl_PHR_50[1,2]
  mltl_95_BA_test_stat <- mltl_95_BA[1,2]
  mltl_50_BA_test_stat <- mltl_50_BA[1,2]
  
  # BFAL between islands
  mbtb_UDOI_95 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="UDOI", percent=95, conditional=T)
  mbtb_UDOI_50 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="UDOI", percent=50, conditional=T)
  mbtb_PHR_95 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="PHR", percent=95, conditional=T)
  mbtb_PHR_50 <- kerneloverlaphr(midBFAL_v_ternBFAL, method="PHR", percent=50, conditional=T)
  mbtb_95_BA <- kerneloverlaphr(midBFAL_v_ternBFAL, method="BA", percent=95, conditional = T) 
  mbtb_50_BA <- kerneloverlaphr(midBFAL_v_ternBFAL, method="BA", percent=50, conditional = T)
  mbtb_UDOI_95_test_stat <- mbtb_UDOI_95[1,2]
  mbtb_UDOI_50_test_stat <- mbtb_UDOI_50[1,2]
  mbtb_PHR_95_test_stat <- mbtb_PHR_95[1,2]
  mbtb_PHR_50_test_stat <- mbtb_PHR_50[1,2]
  mbtb_95_BA_test_stat <- mbtb_95_BA[1,2]
  mbtb_50_BA_test_stat <- mbtb_50_BA[1,2]
  
  # within Midway
  mlmb_UDOI_95 <- kerneloverlaphr(midLAAL_v_midBFAL, method="UDOI", percent=95, conditional=T)
  mlmb_UDOI_50 <- kerneloverlaphr(midLAAL_v_midBFAL, method="UDOI", percent=50, conditional=T)
  mlmb_PHR_95 <- kerneloverlaphr(midLAAL_v_midBFAL, method="PHR", percent=95, conditional=T)
  mlmb_PHR_50 <- kerneloverlaphr(midLAAL_v_midBFAL, method="PHR", percent=50, conditional=T)
  mlmb_95_BA <- kerneloverlaphr(midLAAL_v_midBFAL, method="BA", percent=95, conditional = T) 
  mlmb_50_BA <- kerneloverlaphr(midLAAL_v_midBFAL, method="BA", percent=50, conditional = T)
  mlmb_UDOI_95_test_stat <- mlmb_UDOI_95[1,2]
  mlmb_UDOI_50_test_stat <- mlmb_UDOI_50[1,2]
  mlmb_PHR_95_test_stat <- mlmb_PHR_95[1,2]
  mlmb_PHR_50_test_stat <- mlmb_PHR_50[1,2]
  mlmb_95_BA_test_stat <- mlmb_95_BA[1,2]
  mlmb_50_BA_test_stat <- mlmb_50_BA[1,2]
  
  # within Tern
  tltb_UDOI_95 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="UDOI", percent=95, conditional=T)
  tltb_UDOI_50 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="UDOI", percent=50, conditional=T)
  tltb_PHR_95 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="PHR", percent=95, conditional=T)
  tltb_PHR_50 <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="PHR", percent=50, conditional=T)
  tltb_95_BA <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="BA", percent=95, conditional = T) 
  tltb_50_BA <- kerneloverlaphr(ternLAAL_v_ternBFAL, method="BA", percent=50, conditional = T)
  tltb_UDOI_95_test_stat <- tltb_UDOI_95[1,2]
  tltb_UDOI_50_test_stat <- tltb_UDOI_50[1,2]
  tltb_PHR_95_test_stat <- tltb_PHR_95[1,2]
  tltb_PHR_50_test_stat <- tltb_PHR_50[1,2]
  tltb_95_BA_test_stat <- tltb_95_BA[1,2]
  tltb_50_BA_test_stat <- tltb_50_BA[1,2]
  
  c_comb <- c(j, alab_UDOI_95_test_stat, alab_UDOI_50_test_stat, mltl_UDOI_95_test_stat, mltl_UDOI_50_test_stat, mbtb_UDOI_95_test_stat, mbtb_UDOI_50_test_stat, mlmb_UDOI_95_test_stat, mlmb_UDOI_50_test_stat, tltb_UDOI_95_test_stat, tltb_UDOI_50_test_stat, simtype)
  
  allLAAL <- all_LAAL_averaged_estUD
  allBFAL <- all_BFAL_averaged_estUD
  midLAAL <- lm_averaged_estUD
  midBFAL <- bm_averaged_estUD
  ternLAAL <- lt_averaged_estUD
  ternBFAL <- bt_averaged_estUD
  
  allLAAL_UDOI_95 <- getverticeshr(all_LAAL_averaged_estUD, 95)
  allBFAL_UDOI_95 <- getverticeshr(all_BFAL_averaged_estUD, 95)
  allLAAL_UDOI_50 <- getverticeshr(all_LAAL_averaged_estUD, 50)
  allBFAL_UDOI_50 <- getverticeshr(all_BFAL_averaged_estUD, 50)
  
  allLAAL_UDOI_95$name <- "allLAAL_UDOI_95"
  allBFAL_UDOI_95$name <- "allBFAL_UDOI_95"
  allLAAL_UDOI_50$name <- "allLAAL_UDOI_50"
  allBFAL_UDOI_50$name <- "allBFAL_UDOI_50"
  
  all_simpolys <- rbind(allLAAL_UDOI_95, allLAAL_UDOI_50)
  all_simpolys <- rbind(all_simpolys, allBFAL_UDOI_95)
  all_simpolys <- rbind(all_simpolys, allBFAL_UDOI_50)
  
  midLAAL_UDOI_95 <- getverticeshr(lm_averaged_estUD, 95)
  midBFAL_UDOI_95 <- getverticeshr(bm_averaged_estUD, 95)
  midLAAL_UDOI_50 <- getverticeshr(lm_averaged_estUD, 50)
  midBFAL_UDOI_50 <- getverticeshr(bm_averaged_estUD, 50)
  
  midLAAL_UDOI_95$name <- "midLAAL_UDOI_95"
  midBFAL_UDOI_95$name <- "midBFAL_UDOI_95"
  midLAAL_UDOI_50$name <- "midLAAL_UDOI_50"
  midBFAL_UDOI_50$name <- "midBFAL_UDOI_50"
  
  mid_simpolys <- rbind(midLAAL_UDOI_95, midLAAL_UDOI_50)
  mid_simpolys <- rbind(mid_simpolys, midBFAL_UDOI_95)
  mid_simpolys <- rbind(mid_simpolys, midBFAL_UDOI_50)
  
  ternLAAL_UDOI_95 <- getverticeshr(lt_averaged_estUD, 95)
  ternBFAL_UDOI_95 <- getverticeshr(bt_averaged_estUD, 95)
  ternLAAL_UDOI_50 <- getverticeshr(lt_averaged_estUD, 50)
  ternBFAL_UDOI_50 <- getverticeshr(bt_averaged_estUD, 50)  
  
  ternLAAL_UDOI_95$name <- "ternLAAL_UDOI_95"
  ternBFAL_UDOI_95$name <- "ternBFAL_UDOI_95"
  ternLAAL_UDOI_50$name <- "ternLAAL_UDOI_50"
  ternBFAL_UDOI_50$name <- "ternBFAL_UDOI_50"
  
  tern_simpolys <- rbind(ternLAAL_UDOI_95, ternLAAL_UDOI_50)
  tern_simpolys <- rbind(tern_simpolys, ternBFAL_UDOI_95)
  tern_simpolys <- rbind(tern_simpolys, ternBFAL_UDOI_50)
  
  p_all <- ggplot(data=st_as_sf(allBFAL_UDOI_95)) + 
    geom_sf(data=st_as_sf(allLAAL_UDOI_95), fill="blue", alpha=0.2) +
    geom_sf(data=st_as_sf(allBFAL_UDOI_95), fill="red", alpha=0.2) +
    geom_sf(data=st_as_sf(allLAAL_UDOI_50), fill="blue", alpha=0.6) +
    geom_sf(data=st_as_sf(allBFAL_UDOI_50), fill="red", alpha=0.6) +
    ggtitle(paste("iter:", j, sep=""))
  
  p_mid <- ggplot(data=st_as_sf(midBFAL_UDOI_95)) + 
    geom_sf(data=st_as_sf(midLAAL_UDOI_95), fill="blue", alpha=0.2) +
    geom_sf(data=st_as_sf(midBFAL_UDOI_95), fill="red", alpha=0.2) +
    geom_sf(data=st_as_sf(midLAAL_UDOI_50), fill="blue", alpha=0.6) +
    geom_sf(data=st_as_sf(midBFAL_UDOI_50), fill="red", alpha=0.6) +
    ggtitle(paste("iter:", j, sep=""))
  
  p_tern <- ggplot(data=st_as_sf(ternBFAL_UDOI_95)) + 
    geom_sf(data=st_as_sf(ternLAAL_UDOI_95), fill="blue", alpha=0.2) +
    geom_sf(data=st_as_sf(ternBFAL_UDOI_95), fill="red", alpha=0.2) +
    geom_sf(data=st_as_sf(ternLAAL_UDOI_50), fill="blue", alpha=0.6) +
    geom_sf(data=st_as_sf(ternBFAL_UDOI_50), fill="red", alpha=0.6) +
    ggtitle(paste("iter:", j, sep=""))
  
  #p1 <- ggplot(data=st_as_sf(sim_A)) + 
  #  geom_sf(data=st_as_sf(ver1_95sim), fill="blue", alpha=0.2) +
  #  geom_sf(data=st_as_sf(ver2_95sim), fill="red", alpha=0.2) +
  #  geom_sf(data=st_as_sf(ver1_50sim), fill="blue", alpha=0.6) +
  #  geom_sf(data=st_as_sf(ver2_50sim), fill="red", alpha=0.6) +
  #  ggtitle(paste("300 pts/animal, N=100 each, disp = ", muOL1_s1, sep=""))
  #geom_sf(data=st_as_sf(ann_2_pts), color="red") + 
  #geom_sf(data=st_as_sf(ann_1_pts), color="blue") 
  
  rand_out_pts_birds[[j]] <- sim_CRW_sp
  rand_out_plot_all[[j]] <- p_all
  rand_out_plot_mid[[j]] <- p_mid
  rand_out_plot_tern[[j]] <- p_tern
  
  rand_out_poly_all[[j]] <- all_simpolys
  rand_out_poly_mid[[j]] <- mid_simpolys
  rand_out_poly_tern[[j]] <- tern_simpolys
  #rand_out_blobs[[j]] <- simpolys
  rand_out_table_birds <- rbind(rand_out_table_birds, c_comb)
  print(j)
}

rand_out_table_birds <- data.frame(rand_out_table_birds)
colnames(rand_out_table_birds) <- c("iter", "ALAB95", "ALAB50", "MLTL95", "MLTL50", "MBTB95", "MBTB50", "MLMB95", "MLMB50", "TLTB95", "TLTB50", "simtype")
beep(5)

write.csv(sim_CRW_pts,file="one_iteration_of_CRW_pts.csv")
write.csv(rand_out_table_birds,file="CRW_overlap_results_34to34.csv")
write.csv(big_storage,file="all_CRW_pts_34to34.csv")

#######
#Test stats
test_stats <- t(data.frame("comparison" = c("LAALxBFAL", "midLAALxternLAAL", "midBFALxternBFAL", "midLAALxmidBFAL", "ternLAALxternBFAL"), "95_udoi" = c(0.432, 0.779, 0.411, 0.635, 0.189), "50_udoi" = c(0.001, 0.057, 0.002, 0.004, 0)))
test_stats1 <- melt(test_stats)
test_stats1 <- test_stats1[test_stats1$Var1 != "comparison",]

pvals <- NULL
for(i in 1:nrow(test_stats1)){
  t <- test_stats1[i,3]
  pos <- i + 1
  iters <- as.numeric(as.character(rand_out_table_birds[,pos]))
  vall <- length(iters[iters > t]) / nrow(rand_out_table_birds)
  pvals <- c(pvals, vall)
}

###########################
#totOL_blobs <- rand_out_blobs
#totOL_plot <- rand_out_plot
#totOL_pts <- rand_out_pts
#totOL_table <- rand_out_table
#noOL_teststat <- c_comb

p2 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(ALAB95)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[1,3]))) + 
  ggtitle(paste("All LAALs vs BFALs 95", "p-val =", 1-pvals[1])) +
  xlab("UDOI")

p3 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(ALAB50)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[2,3]))) + 
  ggtitle(paste("All LAALs vs BFALs 50", "p-val =", 1-pvals[2])) +
  xlab("UDOI")

p4 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(MLTL95)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[3,3]))) + 
  ggtitle(paste("Mid vs Tern LAALs 95", "p-val =", 1-pvals[3])) +
  xlab("UDOI")

p5 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(MLTL50)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[4,3]))) + 
  ggtitle(paste("Mid vs Tern LAALs 50", "p-val =", 1-pvals[4])) +
  xlab("UDOI")

p6 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(MBTB95)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[5,3]))) + 
  ggtitle(paste("Mid vs Tern BFALs 95", "p-val =", 1-pvals[5])) +
  xlab("UDOI")

p7 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(MBTB50)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[6,3]))) + 
  ggtitle(paste("Mid vs Tern BFALs 50", "p-val =", 1-pvals[6])) +
  xlab("UDOI")

p8 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(MLMB95)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[7,3]))) + 
  ggtitle(paste("Mid LAALs vs BFALs 95", "p-val =", 1-pvals[7])) +
  xlab("UDOI")

p9 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(MLMB50)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[8,3]))) + 
  ggtitle(paste("Mid LAALs vs BFALs 50", "p-val =", 1-pvals[8])) +
  xlab("UDOI")

p10 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(TLTB95)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[9,3]))) + 
  ggtitle(paste("Tern LAALs vs BFALs 95", "p-val =", 1-pvals[9])) +
  xlab("UDOI")

p11 <- ggplot() + 
  geom_histogram(data=rand_out_table_birds, aes(x=as.numeric(as.character(TLTB50)))) + 
  geom_vline(xintercept = as.numeric(as.character(test_stats1[10,3]))) + 
  ggtitle(paste("Tern LAALs vs BFALs 50", "p-val =", 1-pvals[10])) +
  xlab("UDOI")


cowplot::plot_grid(p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol=2)


#### plot blobs on raster grid
#order goes left 95, right 95, left 50, right 50
totOL_poly <- totOL_blobs[[1]]
for (i in 1:100){
  boo <- totOL_blobs[[i]]
  totOL_poly <- rbind(totOL_poly, boo)
}
totOL_poly <- totOL_poly[-1,]
totOL_poly <- st_as_sf(totOL_poly)

rast <- raster(nrows=100, ncols=100, xmn = -90, xmx = 90, ymn= 5, ymx=160)
rast <- as(rast, "SpatialPolygonsDataFrame")
rast <- st_as_sf(rast)
rast$layer <- rownames(rast)

poly_50 <- totOL_poly[totOL_poly$area < 0.6,]
poly_95 <- totOL_poly[totOL_poly$area >0.8,]

dfL <- data.frame(lon=c(-60,-60,-50,-50,-60), 
                  lat=c(20,50,50,20,20))
polyL <- st_sf(st_sfc(st_polygon(list(as.matrix(dfL)))))
polyL$loc <- "left"

dfR <- data.frame(lon=c(60,60,50,50,60), 
                  lat=c(20,50,50,20,20))
polyR <- st_sf(st_sfc(st_polygon(list(as.matrix(dfR)))))
polyR$loc <- "right"

poly_50_L <- st_join(poly_50, polyL)
poly_50_R <- st_join(poly_50, polyR)
poly_95_L <- st_join(poly_95, polyL)
poly_95_R <- st_join(poly_95, polyR)

poly_50_L <- poly_50_L[!is.na(poly_50_L$loc),]
poly_50_R <- poly_50_R[!is.na(poly_50_R$loc),]
poly_95_L <- poly_95_L[!is.na(poly_95_L$loc),]
poly_95_R <- poly_95_R[!is.na(poly_95_R$loc),]

poly_50_L_rast <- st_intersection(poly_50_L, rast)
poly_50_R_rast <- st_intersection(poly_50_R, rast)
poly_95_L_rast <- st_intersection(poly_95_L, rast)
poly_95_R_rast <- st_intersection(poly_95_R, rast)

poly_50_L_agg <- aggregate(id ~ layer, data = poly_50_L_rast, FUN = length)
poly_50_R_agg <- aggregate(id ~ layer, data = poly_50_R_rast, FUN = length)
poly_95_L_agg <- aggregate(id ~ layer, data = poly_95_L_rast, FUN = length)
poly_95_R_agg <- aggregate(id ~ layer, data = poly_95_R_rast, FUN = length)

rast_50_L <- merge(rast, poly_50_L_agg, by.x="layer", by.y="layer")
rast_50_R <- merge(rast, poly_50_R_agg, by.x="layer", by.y="layer")
rast_95_L <- merge(rast, poly_95_L_agg, by.x="layer", by.y="layer")
rast_95_R <- merge(rast, poly_95_R_agg, by.x="layer", by.y="layer")


ggplot(data=st_as_sf(sim_A)) + 
  geom_sf(data=st_as_sf(ver1_95sim), fill="blue", alpha=0.2) +
  geom_sf(data=st_as_sf(ver2_95sim), fill="red", alpha=0.2) +
  geom_sf(data=st_as_sf(ver1_50sim), fill="blue", alpha=0.6) +
  geom_sf(data=st_as_sf(ver2_50sim), fill="red", alpha=0.6)



#############
soooo <- st_as_sf(sim_CRW_pts, coords = c("Longitude", "Latitude"))
soooo.line <- soooo %>% 
  group_by(Animal_ID) %>%
  arrange(num) %>%
  summarise(Animal_ID = first(Animal_ID), do_union = FALSE) %>%
  st_cast("LINESTRING")

w <- st_as_sf(first_pts, coords=c("Longitude", "Latitude"))

ggplot(data = soooo.line) +
  geom_sf(aes(color = Animal_ID)) +
  theme(legend.position ="") +
  geom_sf(data=w, aes(shape=tag_loc))
