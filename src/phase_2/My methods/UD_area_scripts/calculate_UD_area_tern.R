library(dplyr)
library(ggplot2)
library(maptools)
library(rgdal)
library(GeoLocTools)
setupGeolocation()
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
years <- c("2008","2009","2010","2011","2012")

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  colnames(data_loop)= c("dtime","x","y")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","x","y")]


# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  colnames(data_loop)= c("dtime","x","y")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","x","y")]


# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)

data.spL <- split(data.sp, data.sp$id)

resample_2008 <- rbind(data.spL[["LAAL_2008"]],data.spL[["BFAL_2008"]])
resample_2009 <- rbind(data.spL[["LAAL_2009"]],data.spL[["BFAL_2009"]])
resample_2010 <- rbind(data.spL[["LAAL_2010"]],data.spL[["BFAL_2010"]])
resample_2011 <- rbind(data.spL[["LAAL_2011"]],data.spL[["BFAL_2011"]])
resample_2012 <- rbind(data.spL[["LAAL_2012"]],data.spL[["BFAL_2012"]])

# Replace stuff

resample_list_2012 <- split(resample_2012, resample_2012$id)
current_LAAL <- resample_list_2012[["LAAL_2012"]]
current_BFAL <- resample_list_2012[["BFAL_2012"]]

# don't replace stuff

coordinates(current_LAAL) <- c("x","y")
proj4string(current_LAAL) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
current_LAAL <- spTransform(current_LAAL,CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

coordinates(current_BFAL) <- c("x","y")
proj4string(current_BFAL) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs")
current_BFAL<- spTransform(current_BFAL,CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

kernel.ref <- kernelUD(current_LAAL, same4all = T)
data.kernel.poly <- getverticeshr(kernel.ref, percent = 90, unin = "m", unout = "km2") 
data.kernel.poly$id
data.kernel.poly$area

kernel.ref <- kernelUD(current_BFAL, same4all = T)
data.kernel.poly <- getverticeshr(kernel.ref, percent = 90, unin = "m", unout = "km2") 
data.kernel.poly$id
data.kernel.poly$area

#########################################################
# PRESERVE INDIVIDUAL ID:

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

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  files_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(files_list)){
    files_list[[j]]$iid <- j
  }
  data_loop <- do.call(rbind, files_list)
  colnames(data_loop)= c("dtime","x","y","delete1","delete2","iid")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("dtime","x","y", "iid", "id")]

# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  files_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(files_list)){
    files_list[[j]]$iid <- j
  }
  data_loop <- do.call(rbind, files_list)
  colnames(data_loop)= c("dtime","delete1","delete2","x","y","iid")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("dtime","x","y", "iid", "id")]

# combine into 1 dataframe, preserve year and individual id

data.sp <- rbind(data1.sp,data2.sp)

data.spL <- split(data.sp, data.sp$id)

resample_2008 <- rbind(data.spL[["LAAL_2008"]],data.spL[["BFAL_2008"]])
resample_2009 <- rbind(data.spL[["LAAL_2009"]],data.spL[["BFAL_2009"]])
resample_2010 <- rbind(data.spL[["LAAL_2010"]],data.spL[["BFAL_2010"]])
resample_2011 <- rbind(data.spL[["LAAL_2011"]],data.spL[["BFAL_2011"]])
resample_2012 <- rbind(data.spL[["LAAL_2012"]],data.spL[["BFAL_2012"]])

LAAL_area <- data.frame()
BFAL_area <- data.frame()

# Replace stuff

# replace years when needed: 
current_year = 2012
resample_list_2012 <- split(resample_2012, resample_2012$id)
current_LAAL <- resample_list_2012[["LAAL_2012"]]
current_BFAL <- resample_list_2012[["BFAL_2012"]]

# replace for each iteraiton
iid=7
current_individual_LAAL <- current_LAAL %>% filter(iid==6) # change the number here to get different individuals
current_individual_BFAL <- current_BFAL %>% filter(iid==7) # change the number here to get different individuals

# don't replace stuff

pdc_mercator_proj<-sf::st_crs(3832)

current_individual_LAAL <- current_individual_LAAL[,c("x","y")]
current_individual_LAAL <- sp::SpatialPoints(current_individual_LAAL, proj4string = CRS("+proj=longlat +datum=WGS84"))
current_individual_LAAL <- sp::spTransform(current_individual_LAAL,pdc_mercator_proj$proj4string)

current_individual_BFAL <- current_individual_BFAL[,c("x","y")]
current_individual_BFAL <- sp::SpatialPoints(current_individual_BFAL, proj4string = CRS("+proj=longlat +datum=WGS84"))
current_individual_BFAL <- sp::spTransform(current_individual_BFAL,pdc_mercator_proj$proj4string)

# The grid parameter here needs to be the extent in meters/'grid' = 300km
grid_input <- calculate_sp_obj_extent(current_individual_LAAL)
kernel.ref <- kernelUD(current_individual_LAAL, grid=grid_input, extent = 1)
data.kernel.poly <- getverticeshr(kernel.ref, percent = 95, unin = "m", unout = "km2") 
addLAAL_area<-data.kernel.poly$area

# image(kernel.ref[1])
# title("Output of getvolumeUD")
# xyzv <- as.image.SpatialGridDataFrame(kernel.ref[1])
# contour(xyzv, add=TRUE)

LAAL_area <- rbind(LAAL_area,c(iid,current_year,addLAAL_area))

# The grid parameter here needs to be the extent in meters/'grid' = 300km 
grid_input <- calculate_sp_obj_extent(current_individual_BFAL)
kernel.ref <- kernelUD(current_individual_BFAL, grid=grid_input, extent = 1)
data.kernel.poly <- getverticeshr(kernel.ref, percent = 95, unin = "m", unout = "km2") 
addBFAL_area<-data.kernel.poly$area

BFAL_area <- rbind(BFAL_area,c(iid,current_year,addBFAL_area))

colnames(LAAL_area) <- c("iid","year","area_km2")
colnames(BFAL_area) <- c("iid","year","area_km2")

# when you are done iterating/calculating for each individual

write.csv(LAAL_area, "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/LAAL_area_tern.csv", row.names=F)
write.csv(BFAL_area, "/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/BFAL_area_tern.csv", row.names=F)

# Descriptive plot for these areas: 

## Violin/Boxplots (distribution of covariates by manatee/no manatee)
library(wesanderson)
pal <- wes_palette("Zissou1", 5, type="discrete") # fun colors
LAAL_area$spp = "LAAL"
BFAL_area$spp = "BFAL"

area <- rbind(LAAL_area,BFAL_area)

## All covariates in the same plot
manatee_long <- manateedata %>%  # gather covariates into long format
  pivot_longer(
    cols = depth:mud,
    names_to = "env_variable")

ggplot(area, aes(x=factor(year), y=area_km2, fill=factor(spp))) +
  geom_boxplot() +
  theme_classic() +
  labs(x="year", title ="95th UD area by year and species at Tern", fill="species") + 
  theme(legend.position = "right", legend.text = element_text(size=16), legend.title = element_text(size=16)) +
  scale_fill_manual(values=c(pal[1],pal[3])) +
  theme(plot.title = element_text(size=24,hjust = 0.5), axis.title = element_text(size=20), axis.text = element_text(size=16))
    
  # Make sure title is big!


##########################################
# Trip length box plots

# create trip length data:

trips <- read.csv("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/trip_length_by_individual.csv", header=F)
colnames(trips) <- c("days","year","species","island")
midway_trips <- trips %>% filter(island=="Midway")
tern_trips <- trips %>% filter(island=="Tern")

ggplot(trips, aes(x=factor(year), y=days, fill=factor(species))) +
  geom_boxplot() +
  facet_wrap(~island) +
  theme_classic() +
  labs(x="year", title ="Postbreeding trip duration (days)", fill="species") + 
  theme(legend.position = "right",legend.text = element_text(size=16), legend.title = element_text(size=16)) +
  scale_fill_manual(values=c(pal[1],pal[3])) +
  theme(plot.title = element_text(size=24,hjust = 0.5), axis.title = element_text(size=20), axis.text = element_text(size=16), strip.text = element_text(size=16))

