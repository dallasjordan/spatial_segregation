# Generating data structures for presentation/manuscript figures

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
###################################################
###################################################

# MIDWAY 

###################################################
###################################################
###################################################

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

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/")
save(LAAL, file="LAALdata_midway.Rdata")
save(BFAL, file="BFALdata_midway.Rdata")

###################################################
###################################################
###################################################

# TERN

###################################################
###################################################
###################################################

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
  data_loop <- data_loop[,c("tripid","zm.lon","zm.lat")]
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
  data_loop <- data_loop[,c("tripid","Longitude","Latitude")]
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

LAAL <-rbind(data.spL[["LAAL_2008"]], data.spL[["LAAL_2009"]], data.spL[["LAAL_2010"]], data.spL[["LAAL_2011"]], data.spL[["LAAL_2012"]])
BFAL <-rbind(data.spL[["BFAL_2008"]], data.spL[["BFAL_2009"]], data.spL[["BFAL_2010"]], data.spL[["BFAL_2011"]], data.spL[["BFAL_2012"]])

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/")
save(LAAL, file="LAALdata_tern.Rdata")
save(BFAL, file="BFALdata_tern.Rdata")







###################################################
###################################################

# Repeat the above, but preserving trip ID

###################################################
###################################################

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
years <- c("2008","2009","2010","2011","2012")



###################################################
###################################################
###################################################

# MIDWAY

###################################################
###################################################
###################################################

###################################################
#COMBINE ALL INDIVIDUALS IN A YEAR, PRESERVE INDIVIDUAL ID:

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  file_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(file_list)){
      file_list[[j]]$id <- paste0(years[i],"_",j)
  }
  data_loop <- do.call(rbind, file_list)
  colnames(data_loop)= c("dtime","x","y","track")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","x","y","track")]
data1.xyt <- data1


# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  file_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(file_list)){
    file_list[[j]]$id <- paste0(years[i],"_",j)
  }
  data_loop <- do.call(rbind, file_list)
  colnames(data_loop)= c("dtime","x","y","track")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","x","y","track")]
data2.xyt <- data2


# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)

data.spL <- split(data.sp, data.sp$id)

LAAL <-rbind(data.spL[["LAAL_2008"]], data.spL[["LAAL_2009"]], data.spL[["LAAL_2010"]], data.spL[["LAAL_2011"]], data.spL[["LAAL_2012"]])
BFAL <-rbind(data.spL[["BFAL_2008"]], data.spL[["BFAL_2009"]], data.spL[["BFAL_2010"]], data.spL[["BFAL_2011"]], data.spL[["BFAL_2012"]])

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/")
save(LAAL, file="LAALdata_midway_withTrackID.Rdata")
save(BFAL, file="BFALdata_midway_withTrackID.Rdata")

###################################################
###################################################
###################################################

# TERN

###################################################
###################################################
###################################################

years <- c("2008","2009","2010","2011","2012")

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  file_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(file_list)){
    file_list[[j]]$id <- paste0(years[i],"_",j)
  }
  data_loop <- do.call(rbind, file_list)
  colnames(data_loop)= c("tripid","year","day_gmt","time_gmt","x","y","track")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","x","y","track")]
data1.xyt <- data1

# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  file_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(file_list)){
    file_list[[j]]$id <- paste0(years[i],"_",j)
  }
  data_loop <- do.call(rbind, file_list)
  colnames(data_loop)= c("tripid","year","day_gmt","time_gmt","x","y","track")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","x","y","track")]
data2.xyt <- data2

# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)

data.spL <- split(data.sp, data.sp$id)

LAAL <-rbind(data.spL[["LAAL_2008"]], data.spL[["LAAL_2009"]], data.spL[["LAAL_2010"]], data.spL[["LAAL_2011"]], data.spL[["LAAL_2012"]])
BFAL <-rbind(data.spL[["BFAL_2008"]], data.spL[["BFAL_2009"]], data.spL[["BFAL_2010"]], data.spL[["BFAL_2011"]], data.spL[["BFAL_2012"]])

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/")
save(LAAL, file="LAALdata_tern_withTrackID.Rdata")
save(BFAL, file="BFALdata_tern_withTrackID.Rdata")


###################################################
###################################################

# Repeat the above, but preserving trip ID AND DATE

###################################################
###################################################

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
years <- c("2008","2009","2010","2011","2012")



###################################################
###################################################
###################################################

# MIDWAY

###################################################
###################################################
###################################################

###################################################
#COMBINE ALL INDIVIDUALS IN A YEAR, PRESERVE INDIVIDUAL ID:

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  file_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(file_list)){
    file_list[[j]]$id <- paste0(years[i],"_",j)
  }
  data_loop <- do.call(rbind, file_list)
  colnames(data_loop)= c("dtime","x","y","track")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","dtime","x","y","track")]
data1.xyt <- data1


# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  file_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(file_list)){
    file_list[[j]]$id <- paste0(years[i],"_",j)
  }
  data_loop <- do.call(rbind, file_list)
  colnames(data_loop)= c("dtime","x","y","track")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","dtime","x","y","track")]
data2.xyt <- data2


# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)

data.spL <- split(data.sp, data.sp$id)

LAAL <-rbind(data.spL[["LAAL_2008"]], data.spL[["LAAL_2009"]], data.spL[["LAAL_2010"]], data.spL[["LAAL_2011"]], data.spL[["LAAL_2012"]])
BFAL <-rbind(data.spL[["BFAL_2008"]], data.spL[["BFAL_2009"]], data.spL[["BFAL_2010"]], data.spL[["BFAL_2011"]], data.spL[["BFAL_2012"]])

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/")
save(LAAL, file="LAALdata_midway_withTrackID_dtime.Rdata")
save(BFAL, file="BFALdata_midway_withTrackID_dtime.Rdata")

###################################################
###################################################
###################################################

# TERN

###################################################
###################################################
###################################################

years <- c("2008","2009","2010","2011","2012")

# Make LAAL data object

SPECIES = "LAAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data1<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  file_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(file_list)){
    file_list[[j]]$id <- paste0(years[i],"_",j)
  }
  data_loop <- do.call(rbind, file_list)
  colnames(data_loop)= c("tripid","year","day_gmt","time_gmt","x","y","track")
  data_loop$id <- paste0("LAAL_",years[i])
  data1 <- rbind(data1,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
data1.sp <- data1[,c("id","day_gmt","time_gmt","x","y","track")]
data1.xyt <- data1

# Make BFAL data object

SPECIES = "BFAL"

setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
             SPECIES,"/"))
data2<-data.frame()
for (i in 1:5){
  files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
  setwd(paste0(getwd(),"/",years[i]))
  file_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
  for (j in 1:length(file_list)){
    file_list[[j]]$id <- paste0(years[i],"_",j)
  }
  data_loop <- do.call(rbind, file_list)
  colnames(data_loop)= c("tripid","year","day_gmt","time_gmt","x","y","track")
  data_loop$id <- paste0("BFAL_",years[i])
  data2 <- rbind(data2,data_loop)
  setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
               SPECIES,"/"))
}

data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
data2.sp <- data2[,c("id","day_gmt","time_gmt","x","y","track")]
data2.xyt <- data2

# combine into 1 dataframe, preserve id

data.sp <- rbind(data1.sp,data2.sp)

data.spL <- split(data.sp, data.sp$id)

LAAL <-rbind(data.spL[["LAAL_2008"]], data.spL[["LAAL_2009"]], data.spL[["LAAL_2010"]], data.spL[["LAAL_2011"]], data.spL[["LAAL_2012"]])
BFAL <-rbind(data.spL[["BFAL_2008"]], data.spL[["BFAL_2009"]], data.spL[["BFAL_2010"]], data.spL[["BFAL_2011"]], data.spL[["BFAL_2012"]])

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/")
save(LAAL, file="LAALdata_tern_withTrackID_dtime.Rdata")
save(BFAL, file="BFALdata_tern_withTrackID_dtime.Rdata")

