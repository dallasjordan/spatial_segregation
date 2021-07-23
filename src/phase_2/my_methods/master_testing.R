# testing June 11 2021

# functions and objects
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
years <- c("2008","2009","2010","2011","2012")

#########################
###### LOAD DATA ########
#########################

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_midway_withTrackID.Rdata")
LAALmid <- LAAL
LAALmid$id <- paste0("lm",LAALmid$track)
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_tern_withTrackID.Rdata")
LAALtern <- LAAL
LAALtern$id <- paste0("lt",LAALtern$track)

LAAL <- rbind(LAALmid, LAALtern)

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_midway_withTrackID.Rdata")
BFALmid <- BFAL
BFALmid$id <- paste0("bm",BFALmid$track)
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_tern_withTrackID.Rdata")
BFALtern <- BFAL
BFALtern$id <- paste0("bt",BFALtern$track)

BFAL <- rbind(BFALmid, BFALtern)

# comparisons: allLAAL v all BFAL
#              ternLAAL v midwayLAAL
#              ternBFAL v midwayBFAL
#              ternLAAL v ternBFAL
#              midwayLAAL v midwayBFAL
all_data <- rbind(LAAL,BFAL)
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

#########################
#### CALC GRID VALUE ####
#########################

all_data_1<-all_data[,1:3]
sp::coordinates(all_data_1) <- c("x", "y")
proj4string(all_data_1) <- CRS("+proj=longlat +zone=2 +datum=WGS84 +no_defs") # placeholder
all_data_1 <- spTransform(all_data_1,laea)
grid_input <- calculate_sp_obj_extent(all_data_1)
# https://trac.faunalia.it/animove/wiki/AnimoveHowto#Whatunitisthegridcellsizein
# some more good notes about parameters: 

###########################
# GENERATE INDIVIDUAL KDE #
###########################

# kernelUD output currently is probability density: absolute density (event per unit area) divided by total number of events
results <- kernelUD(all_data_1, grid=grid_input,same4all=T,extent=0.1) # grid-input guarantees AT LEAST 300km x 300km - in actuality, it is more than 300km x 300km.
# image(results)

kerneloverlaphr(results)

names <- names(results)
lm_kde <- results[grep("lm",names)]
lt_kde <- results[grep("lt",names)]
bm_kde <- results[grep("bm",names)]
bt_kde <- results[grep("bt",names)]

