# Randomization procedure for aggregated (all years both colonies by species) data
# Dallas Jordan
# Feb 23 2021, last edited March 29th 2021
# Changes March 2021: 50th and 95th isolpleths, permuting labels of trips rather than aggregated Lat/Longs, normalized UD
# ascgen() function lets you make own raster grid to pass to kernelUD, https://lists.faunalia.it/pipermail/animov/2006-May/000132.html

library(dplyr)
library(maptools)
library(rgdal)
library(GeoLocTools)
setupGeolocation()
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)
library(sf)

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

pdc_mercator_proj<-"+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
laea <- "+proj=laea +lat_0=20 +lon_0=180 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# load aggregated data

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_midway_withTrackID.Rdata")
LAALmid <- LAAL
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_tern_withTrackID.Rdata")
LAALtern <- LAAL

LAAL <- rbind(LAALmid, LAALtern)
LAAL$id <- "LAAL"

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_midway_withTrackID.Rdata")
BFALmid <- BFAL
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_tern_withTrackID.Rdata")
BFALtern <- BFAL

BFAL <- rbind(BFALmid, BFALtern)
BFAL$id <- "BFAL"

both <- rbind(LAAL,BFAL)

    ##### Need to calculate a single grid for both:
    both_1<-both[,1:3]
    sp::coordinates(both_1) <- c("x", "y")
    
    proj4string(both_1) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs") # placeholder
    both_1 <- spTransform(both_1,laea)
    grid_input <- calculate_sp_obj_extent(both_1)
    results <- kernelUD(both_1, grid=grid_input,same4all=T,extent=0.1)
    image(results)
    vud_points<-getvolumeUD(results) # for diagnosing extent parameter
    image(vud_points)
    original_results_95 <- kerneloverlaphr(results, method="UDOI", percent=95)
    original_results_50 <- kerneloverlaphr(results, method="UDOI", percent=50)
    original_BA_95 <- kerneloverlaphr(results, method="BA", percent=95)
    original_BA_50 <- kerneloverlaphr(results, method="BA", percent=50)
    BA_95 <- original_BA_95 # 	0.6943703
    BA_50 <- original_BA_50 #   0.6943703
    test_stat_95 <- 0.4673984 # manually putting this in because its easier then calling from matrix
    test_stat_50 <- 0.008576117 # manually putting this in because its easier then calling from matrix
    
    ################### Export Raster and contour lines SPDF
    names(results)
    r<- raster(as(results$LAAL,"SpatialPixelsDataFrame"))
    r2<- raster(as(results$BFAL,"SpatialPixelsDataFrame"))
    writeRaster(r, filename = "Rast_LAAL_aggregated.tif", options=c('TFW=YES'), overwrite=TRUE)
    
    
    x<-getverticeshr(results$BFAL, percent=95)
    y<-getverticeshr(results$LAAL, percent=95)
    x1<-getverticeshr(results$BFAL, percent=50)
    y1<-getverticeshr(results$LAAL, percent=50)
    rgdal::writeOGR(y1, getwd(), "LAALaggregated50" ,driver="ESRI Shapefile", verbose=TRUE, overwrite=TRUE) 
    
    ################### Normalize UDs
    s <- sum(values(r), na.rm=TRUE)
    values(r) <- values(r)/s
    df <- data.frame(v=values(r), id=seq(1,ncell(r)))
    df.2 <- df[rev(order(df$v)),]
    df.2$cs <- cumsum(df.2[,'v'])
    
    #CREATE BREAKS
    df.2$c<-cut(df.2$cs,breaks=seq(0,1,by=0.05),labels=seq(100,5,by=-5))
    df.3 <- df.2[order(df.2$id),]
    values(r) <- as.numeric(df.3$c)
    plot(r)
    
############ Step 2: randomization
both_2 <- both 
resample_this <- both_2
resample_this_sp <- split(resample_this, resample_this$id)
resample_this_BFAL_track <- split(resample_this_sp[[1]],resample_this_sp[[1]]$track)
resample_this_LAAL_track <- split(resample_this_sp[[2]],resample_this_sp[[2]]$track)
resample_all_tracks <- c(resample_this_BFAL_track,resample_this_LAAL_track)
available_BFAL_tracks <- unique(resample_this_sp[[1]]$track)
available_LAAL_tracks <- unique(resample_this_sp[[2]]$track)
n_abt <- length(available_BFAL_tracks) # you are calculating these to preserve sample size in randomizations
n_alt <- length(available_LAAL_tracks)
n_all <- n_abt+n_alt
counts <- both %>% count(id,track)

#95th isopleth

iter <- 100
significance_tally <- 0
iter_tally <- 0
results_storage <- vector()
for (i in 1:iter){
  pool <- 1:n_all
  BFAL_nums <- sample(pool, n_abt, replace=F)
  LAAL_nums <- setdiff(pool,BFAL_nums)
  BFAL_iter <- NULL
  LAAL_iter <- NULL
  for (j in 1:length(BFAL_nums)){
    BFAL_loop <- resample_all_tracks[[BFAL_nums[j]]]
    BFAL_iter <- rbind(BFAL_iter,BFAL_loop)
  }
  BFAL_iter <- BFAL_iter[,c("id","x","y")]
  BFAL_iter$id <- "BFALiter"
  for (j in 1:length(LAAL_nums)){
    LAAL_loop <- resample_all_tracks[[LAAL_nums[j]]]
    LAAL_iter <- rbind(LAAL_iter,LAAL_loop)
  }
  LAAL_iter <- LAAL_iter[,c("id","x","y")]
  LAAL_iter$id <- "LAALiter"
  final_iter <- rbind(BFAL_iter, LAAL_iter)
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  coordinates(final_iter) <- c("x","y")
  proj4string(final_iter) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  final_iter <- spTransform(final_iter,laea)
  
  kernel.boot <- kernelUD(final_iter, grid=grid_input, same4all = T, extent=0.1)
  
  results<-kerneloverlaphr(kernel.boot, method="UDOI", percent=95)
  results_storage[i] <- results[1,2]
  print(results_storage[i])
  if (results[1,2]<test_stat_95){
    print(results)
    print("+1!")
    significance_tally <- significance_tally+1
  }
}

p_value <- significance_tally/iter
mean_value <- mean(results_storage)
sd_value <- sd(results_storage)

p_value # p<0.001
mean_value # 1.257934
sd_value # 0.04705332

#50th isopleth

iter <- 100
significance_tally <- 0
iter_tally <- 0
results_storage <- vector()
for (i in 1:iter){
  pool <- 1:n_all
  BFAL_nums <- sample(pool, n_abt, replace=F)
  LAAL_nums <- setdiff(pool,BFAL_nums)
  BFAL_iter <- NULL
  LAAL_iter <- NULL
  for (j in 1:length(BFAL_nums)){
    BFAL_loop <- resample_all_tracks[[BFAL_nums[j]]]
    BFAL_iter <- rbind(BFAL_iter,BFAL_loop)
  }
  BFAL_iter <- BFAL_iter[,c("id","x","y")]
  BFAL_iter$id <- "BFALiter"
  for (j in 1:length(LAAL_nums)){
    LAAL_loop <- resample_all_tracks[[LAAL_nums[j]]]
    LAAL_iter <- rbind(LAAL_iter,LAAL_loop)
  }
  LAAL_iter <- LAAL_iter[,c("id","x","y")]
  LAAL_iter$id <- "LAALiter"
  final_iter <- rbind(BFAL_iter, LAAL_iter)
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  coordinates(final_iter) <- c("x","y")
  proj4string(final_iter) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  final_iter <- spTransform(final_iter,laea)
  
  kernel.boot <- kernelUD(final_iter, grid=grid_input, same4all = T, extent=0.1)
  
  results<-kerneloverlaphr(kernel.boot, method="UDOI", percent=50)
  results_storage[i] <- results[1,2]
  print(results_storage[i])
  if (results[1,2]<test_stat_50){
    print(results)
    print("+1!")
    significance_tally <- significance_tally+1
  }
}

p_value <- significance_tally/iter
mean_value <- mean(results_storage)
sd_value <- sd(results_storage)

p_value # p<0.001
mean_value # 0.2756724
sd_value # 0.02158613












# old - permutes lats and longs instead of just the ID associated with each point. Don't use!
significance_tally <- 0
iter_tally <- 0
results_storage <- vector()
for (i in 1:iter){
  resample_this <- both
  resample_this$x <- sample(resample_this$x,length(resample_this$x), replace=F)
  resample_this$y <- sample(resample_this$y,length(resample_this$y), replace=F)
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  coordinates(resample_this) <- c("x","y")
  proj4string(resample_this) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
  resample_this <- spTransform(resample_this,pdc_mercator_proj)
  
  kernel.boot <- kernelUD(resample_this, grid=grid_input, same4all = T)
  
  results<-kerneloverlaphr(kernel.boot, method="UDOI", percent=95)
  results_storage[i] <- results[1,2]
  if (results[1,2]<test_stat){
    print(results)
    print("+1!")
    significance_tally <- significance_tally+1
  }
}

p_value <- significance_tally/iter
mean_value <- mean(results_storage)
sd_value <- sd(results_storage)

year_im_on
p_value
mean_value
sd_value



# Ultimately didn't use this, leaving for future reference just in case!
    
# coordinates(LAAL) <- c("x","y")
# LAAL <- sp::SpatialPoints(LAAL, proj4string = CRS("+proj=longlat +datum=WGS84"))
# LAAL <- sp::spTransform(LAAL,pdc_mercator_proj)
# coordinates(BFAL) <- c("x","y")
# BFAL <- sp::SpatialPoints(BFAL, proj4string = CRS("+proj=longlat +datum=WGS84"))
# BFAL <- sp::spTransform(BFAL,pdc_mercator_proj)
# 
# ########### Calculate kernel and overlap of original data: 
# 
# LAALkernel <- kernelUD(LAAL, grid=grid_input, extent = 1)
# 
# BFALkernel <- kernelUD(BFAL, grid=grid_input, extent = 1)
# 
# combined <- list(LAALkernel, BFALkernel)
# class(combined) <- "estUDm"
# 
# original_results <- kerneloverlaphr(combined, method="PHR", percent=95)

########### Original script, used to adapt what is above:

year_im_on <- 2008
test_stat <- 0.488580284433701 # manually putting this in because its easier then writing code to get it from other script
iter <- 1000

significance_tally <- 0
iter_tally <- 0
results_storage <- vector()
for (i in 1:iter){
  resample_this <- resample_2012
  resample_this$x <- sample(resample_this$x,length(resample_this$x), replace=F)
  resample_this$y <- sample(resample_this$y,length(resample_this$y), replace=F)
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  coordinates(resample_this) <- c("x","y")
  proj4string(resample_this) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
  resample_this <- spTransform(resample_this,CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  kernel.boot <- kernelUD(resample_this, same4all = T)
  
  results<-kerneloverlaphr(kernel.boot, method="PHR", percent=95)
  results_storage[i] <- results[1,2]
  if (results[1,2]<test_stat){
    print(results)
    print("+1!")
    significance_tally <- significance_tally+1
  }
}

p_value <- significance_tally/iter
mean_value <- mean(results_storage)
sd_value <- sd(results_storage)

year_im_on
p_value
mean_value
sd_value