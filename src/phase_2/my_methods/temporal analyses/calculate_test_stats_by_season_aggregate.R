# Calculate test stats by season
# Derived from master_script.R, taking the portion where I calculate UDs and the overlap stats
# of each category. 
# May 1 2022
# Dallas Jordan

# Setup -------------------------------------------------------------------

# packages
library(dplyr)
library(lubridate)
library(maptools)
library(rgdal)
# library(GeoLocTools)
# setupGeolocation()
# Not available for R>=4.0
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)
library(sf)
library(ggplot2)

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
lcea <- "+proj=cea +lat_0=0 +lat_ts=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs" 
years <- c("2008","2009","2010","2011","2012")


# Load data ---------------------------------------------------------------
# USE RELATIVE PATHWAYS - setwd ONCE and leave it! 
setwd("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/")
# list Tern IS. BFAL FILES
d_spps_bfal_tern <- list.files(path = "./tern_postbreeding_exports/import_for_CRW/BFAL/", pattern = ".csv$", recursive=T, full.names = T)
d_spps_bfal_tern <- lapply(d_spps_bfal_tern, FUN = read.csv, header = T)

# list Tern IS. LAAL Files
d_spps_laal_tern <- list.files(path = "./tern_postbreeding_exports/import_for_CRW/LAAL/", pattern = ".csv$", recursive=T, full.names = T)
d_spps_laal_tern <- lapply(d_spps_laal_tern, FUN = read.csv, header=T)

d_spps_bfal_mid <- list.files(path = "./midway_postbreeding_exports/BFAL/",pattern = ".csv$", recursive=T, full.names = T)
d_spps_bfal_mid <- lapply(d_spps_bfal_mid, FUN = read.csv, header = T)
d_spps_laal_mid <- list.files(path = "./midway_postbreeding_exports/LAAL/", pattern = ".csv$", recursive=T, full.names = T)
d_spps_laal_mid <- lapply(d_spps_laal_mid, function(x) read.csv(x, header=T))

#BFAL TERN
d_spps_bfal_ternID <- c()
for(i in 1:length(d_spps_bfal_tern)){ # never hardcode for loop lengths, your code won't be applicable to any other data
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
for(i in 1:length(d_spps_laal_tern)){
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
for(i in 1:length(d_spps_bfal_mid)){
  w <- d_spps_bfal_mid[[i]]
  w$ID <- paste("BFAL_MID_", i, sep="")
  w$dtime <- ymd_hms(w$dtime)
  ee <- w %>% dplyr::select(ID, dtime, Longitude, Latitude)
  colnames(ee) <- colnames(d_spps_bfal_ternID)
  d_spps_bfal_midID <- rbind(d_spps_bfal_midID, ee)
}

#LAAL MID
d_spps_laal_midID <- c()
for(i in 1:length(d_spps_laal_mid)){
  w <- d_spps_laal_mid[[i]]
  w$ID <- paste("LAAL_MID_", i, sep="")
  w$dtime <- ymd_hms(w$dtime)
  ee <- w %>% dplyr::select(ID, dtime, Longitude, Latitude)
  colnames(ee) <- colnames(d_spps_bfal_ternID)
  d_spps_laal_midID <- rbind(d_spps_laal_midID, ee)
}

# make dataframe of all categories with IDs
comb <- rbind(d_spps_bfal_midID, d_spps_bfal_ternID, d_spps_laal_midID, d_spps_laal_ternID)

### HERE: BREAK INTO SEASONS ### 

comb$month <- month(comb$GMT)
comb <- comb %>% mutate(season=case_when(
  month==1 ~ 1,
  month==2 ~ 1,
  month==3 ~ 2,
  month==4 ~ 2,
  month==5 ~ 2,
  month==6 ~ 3,
  month==7 ~ 3,
  month==8 ~ 3,
  month==9 ~ 4,
  month==10 ~ 4,
  month==11 ~ 4,
  month==12 ~ 1
))
comb <- comb %>% mutate(season_name=case_when(
  season==1 ~ "Winter",
  season==2 ~ "Spring",
  season==3 ~ "Summer",
  season==4 ~ "Fall"
))
comb <- comb %>% mutate(month = month.name[month])

comb_spring <- comb %>% filter(month== "June" |  month=="July")
comb_summer <- comb %>% filter(month== "August" |  month=="September")
comb_fall <- comb %>% filter(month== "October" |  month=="November")
  
### SELECT WHAT YOU WANT GOING FORWARD ###

# ALL 
all_data <- comb
lm <- all_data[grep("LAAL_MID", all_data$ID), ]
lt <- all_data[grep("LAAL_TERN", all_data$ID), ]
bm <- all_data[grep("BFAL_MID", all_data$ID), ]
bt <- all_data[grep("BFAL_TERN", all_data$ID), ]

# SPRING
all_data <- comb_spring
lm <- all_data[grep("LAAL_MID", all_data$ID), ]
lt <- all_data[grep("LAAL_TERN", all_data$ID), ]
bm <- all_data[grep("BFAL_MID", all_data$ID), ]
bt <- all_data[grep("BFAL_TERN", all_data$ID), ]

# SUMMER
all_data <- comb_summer
lm <- all_data[grep("LAAL_MID", all_data$ID), ]
lt <- all_data[grep("LAAL_TERN", all_data$ID), ]
bm <- all_data[grep("BFAL_MID", all_data$ID), ]
bt <- all_data[grep("BFAL_TERN", all_data$ID), ]

# FALL 
all_data <- comb_fall
lm <- all_data[grep("LAAL_MID", all_data$ID), ]
lt <- all_data[grep("LAAL_TERN", all_data$ID), ]
bm <- all_data[grep("BFAL_MID", all_data$ID), ]
bt <- all_data[grep("BFAL_TERN", all_data$ID), ]

# Calculate grid value ----------------------------------------------------

all_data <- comb_fall
all_data <- all_data %>% group_by(ID) %>% mutate(n=n()) %>% filter(n>5)

all_data_1<-all_data[,c(1,3,4)]
sp::coordinates(all_data_1) <- c("Longitude", "Latitude")
proj4string(all_data_1) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # placeholder
all_data_1 <- spTransform(all_data_1,lcea)
grid_input <- calculate_sp_obj_extent(all_data_1,0.1)

# Effect of grid size, however, is negligible...
# http://r-sig-geo.2731867.n2.nabble.com/adehabitatHR-kerneloverlap-and-kernelUD-etc-td7586144.html

# https://trac.faunalia.it/animove/wiki/AnimoveHowto#Whatunitisthegridcellsizein
# some more good notes on parameters: 
# https://lists.faunalia.it/pipermail/animov/2006-May/000137.html


# Generate individual KDE -------------------------------------------------

# kernelUD output currently is probability density: absolute density (events per unit area) divided by total number of events
results <- kernelUD(all_data_1, grid=grid_input,same4all=T,extent=0.1,h=150000)# grid-input guarantees AT LEAST 300km x 300km - in actuality, it is very slightly more than 300km x 300km.
image(results[[1]])

names <- names(results)
lm_kde <- results[grep("LAAL_MID",names)]
lt_kde <- results[grep("LAAL_TERN",names)]
if (nrow(bm)>0){
  bm_kde <- results[grep("BFAL_MID",names)]
}
bt_kde <- results[grep("BFAL_TERN",names)]
all_LAAL <- append(lm_kde, lt_kde)
if (exists("bm_kde")){
  all_BFAL <- append(bm_kde, bt_kde)
} else {
  all_BFAL <- bt_kde
}


# Averaging -> new KDEs, account for track length ------------------------

# "0" values are preserved here, not converted to NA
# lt_holder length will change here if your cell size changes (calc_sp_obj_extent)

len <- nrow(results[[1]]@data)

# average lm #
lm_holder <- numeric(length = len)
for (i in 1:length(lm_kde)) {
  lm_kde[[i]]@data$ud[is.na(lm_kde[[i]]@data$ud)] <- 0
  add <- lm_kde[[i]]@data$ud
  lm_holder <- lm_holder+add
}
lm_holder <- lm_holder/length(lm_kde)
lm_holder

##### modify existing estUD object with averaged values, then rename
lm_averaged_estUD <- lm_kde[[1]]
lm_averaged_estUD@data$ud <- lm_holder
lm_averaged_estUD@data$ud[is.na(lm_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(lm_kde[[1]])
image(lm_averaged_estUD) 

# average lt #
lt_holder <- numeric(length = len)
for (i in 1:length(lt_kde)) {
  lt_kde[[i]]@data$ud[is.na(lt_kde[[i]]@data$ud)] <- 0
  add <- lt_kde[[i]]@data$ud
  lt_holder <- lt_holder+add
}
lt_holder <- lt_holder/length(lt_kde)
lt_holder

##### modify existing estUD object with averaged values, then rename
lt_averaged_estUD <- lt_kde[[1]]
lt_averaged_estUD@data$ud <- lt_holder
lt_averaged_estUD@data$ud[is.na(lt_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(lt_kde[[1]])
image(lt_averaged_estUD)         

# average bm #
if (exists("bm_kde")){
  bm_holder <- numeric(length = len)
  for (i in 1:length(bm_kde)) {
    bm_kde[[i]]@data$ud[is.na(bm_kde[[i]]@data$ud)] <- 0
    add <- bm_kde[[i]]@data$ud
    bm_holder <- bm_holder+add
  }
  bm_holder <- bm_holder/length(bm_kde)
  bm_holder
  
  ##### modify existing estUD object with averaged values, then rename
  bm_averaged_estUD <- bm_kde[[1]]
  bm_averaged_estUD@data$ud <- bm_holder
  bm_averaged_estUD@data$ud[is.na(bm_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
  image(bm_kde[[1]])
  image(bm_averaged_estUD)
}

# average bt #
bt_holder <- numeric(length = len)
for (i in 1:length(bt_kde)) {
  bt_kde[[i]]@data$ud[is.na(bt_kde[[i]]@data$ud)] <- 0
  add <- bt_kde[[i]]@data$ud
  bt_holder <- bt_holder+add
}
bt_holder <- bt_holder/length(bt_kde)
bt_holder

##### modify existing estUD object with averaged values, then rename
bt_averaged_estUD <- bt_kde[[1]]
bt_averaged_estUD@data$ud <- bt_holder
bt_averaged_estUD@data$ud[is.na(bt_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(bt_kde[[1]])
image(bt_averaged_estUD)

# average all_LAAL #
all_LAAL_holder <- numeric(length = len)
for (i in 1:length(all_LAAL)) {
  all_LAAL[[i]]@data$ud[is.na(all_LAAL[[i]]@data$ud)] <- 0
  add <- all_LAAL[[i]]@data$ud
  all_LAAL_holder <- all_LAAL_holder+add
}
all_LAAL_holder <- all_LAAL_holder/length(all_LAAL)
all_LAAL_holder

##### modify existing estUD object with averaged values, then rename
all_LAAL_averaged_estUD <- all_LAAL[[1]]
all_LAAL_averaged_estUD@data$ud <- all_LAAL_holder
all_LAAL_averaged_estUD@data$ud[is.na(all_LAAL_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(all_LAAL[[1]])
image(all_LAAL_averaged_estUD)

# average all_BFAL #
all_BFAL_holder <- numeric(length = len)
for (i in 1:length(all_BFAL)) {
  all_BFAL[[i]]@data$ud[is.na(all_BFAL[[i]]@data$ud)] <- 0
  add <- all_BFAL[[i]]@data$ud
  all_BFAL_holder <- all_BFAL_holder+add
}
all_BFAL_holder <- all_BFAL_holder/length(all_BFAL)
all_BFAL_holder

##### modify existing estUD object with averaged values, then rename
all_BFAL_averaged_estUD <- all_BFAL[[1]]
all_BFAL_averaged_estUD@data$ud <- all_BFAL_holder
all_BFAL_averaged_estUD@data$ud[is.na(all_BFAL_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(all_BFAL[[1]])
image(all_BFAL_averaged_estUD)    


# Generate class KDE ------------------------------------------------------

# reminder for converting into a spatial pixel dataframe if needed: 
# udspdf <- estUDm2spixdf(lm_kde)

allLAAL <- all_LAAL_averaged_estUD
allBFAL <- all_BFAL_averaged_estUD
midLAAL <- lm_averaged_estUD
if (exists("bm_averaged_estUD")){
  midBFAL <- bm_averaged_estUD
}
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
if (exists("midBFAL")){
  midBFAL_v_ternBFAL <- list(midBFAL,ternBFAL)
  class(midBFAL_v_ternBFAL)<-"estUDm"
  names(midBFAL_v_ternBFAL)<-c("midBFAL","ternBFAL")
}

# within Midway
if (exists("midBFAL")){
  midLAAL_v_midBFAL <- list(midLAAL,midBFAL)
  class(midLAAL_v_midBFAL)<-"estUDm"
  names(midLAAL_v_midBFAL)<-c("midLAAL","midBFAL")
}

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

path <- "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/pre_defense/permutation_tests/test_stats/"
setwd(path)
# just save everything out and load into other script (calculate_temporal_pvalues) to visualize the test stats easier
save(list = ls(all.names = TRUE), file = "all_test_stats_octobernovember.RData")

