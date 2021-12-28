# Plotting by island/species/season
# Dallas Jordan
# Nov 15 2021

# Purpose of this script is to make 2 4-panel plots, one for each species, of modelled locations
# in the North Pacific. Uses Ptolemy, ggarrange. Each one of 4 panels is a different season. 
# First load in all the tracking data, add months, add season counter, then plot


# Setup -------------------------------------------------------------------

library(mapproj)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgdal)
library(ggspatial)
library(ggmap)
library(ggthemes)
library(stars)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(dplyr)
library(ggnewscale)
library(ggsn)

# CRITICAL: 
sf::sf_use_s2(FALSE)
# The above line is necessary as of April 17 2021, when 'sf' switched to using S2 spheroid 
# geometry instead of a GEOS routine that assumed projected coordinates. Until 'mapping' 
# updates its back-end, switching S2 geometry off is necessary for my script to work. 
# See: https://github.com/r-spatial/sf/issues/1649
# Great reference here: https://r-spatial.github.io/sf/articles/sf7.html


# Load in tracking data ---------------------------------------------------

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

mid_data <- rbind(lm,bm)
tern_data <-rbind(lt,bt)

# fix some dates for later
mid_data$dtime <- as.Date(mid_data$dtime)
tern_data$dtime <- as.Date(tern_data$dtime)
lm$dtime <- as.Date(lm$dtime)
bm$dtime <- as.Date(bm$dtime)
lt$dtime <- as.Date(lt$dtime)
bt$dtime <- as.Date(bt$dtime)

lm$spp <- "LAAL"
lm$island <- "Midway"
bm$spp <- "BFAL"
bm$island <- "Midway"
lt$spp <- "LAAL"
lt$island <- "Tern"
bt$spp <- "BFAL"
bt$island <- "Tern"

LAAL_data <- rbind(lm,lt)
BFAL_data <- rbind(bm,bt)

LAAL_data$month <- month(LAAL_data$dtime)
LAAL_data <- LAAL_data %>% mutate(season=case_when(
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
LAAL_data <- LAAL_data %>% mutate(season_name=case_when(
  season==1 ~ "Winter",
  season==2 ~ "Spring",
  season==3 ~ "Summer",
  season==4 ~ "Fall"
))
LAAL_data <- LAAL_data %>% mutate(month = month.name[month])

BFAL_data$month <- month(BFAL_data$dtime)
BFAL_data <- BFAL_data %>% mutate(season=case_when(
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
BFAL_data <- BFAL_data %>% mutate(season_name=case_when(
  season==1 ~ "Winter",
  season==2 ~ "Spring",
  season==3 ~ "Summer",
  season==4 ~ "Fall"
))
BFAL_data <- BFAL_data %>% mutate(month = month.name[month])

# filtering out to just peak postbreeding months, other months don't have enough data:
LAAL_points_sufficient_data <- LAAL_data %>% filter(month %in% c("June","July","August","September","October","November"))
BFAL_points_sufficient_data <- BFAL_data %>% filter(month %in% c("June","July","August","September","October","November"))


# LAAL_points_sufficient_data <- LAAL_data %>% filter(month %in% c(6:11))
# BFAL_points_sufficient_data <- BFAL_data %>% filter(month %in% c(6:11))

LAAL_points_sufficient_data <- LAAL_points_sufficient_data %>% 
  mutate(month_factor = factor(month, levels=c("June","July","August","September","October","November")))
BFAL_points_sufficient_data <- BFAL_points_sufficient_data %>% 
  mutate(month_factor = factor(month, levels=c("June","July","August","September","October","November")))
# check to see if Ptolemy needs 0 to 360 or not
# convlon <- function(lon){
#   ifelse(lon > 180, lon - 360, lon)
# }
# 
# lm$x <- convlon(lm$x)
# lt$x <- convlon(lt$x)
# bm$x <- convlon(bm$x)
# bt$x <- convlon(bt$x)

# Basemap -----------------------------------------------------------------


# Generate basemap 

### TWO METHODS ###
# Ptolemy package, formerly "nPacMaps". Very amazing for North Pacific Mapping

# First method, generate a basemap based on the extent and crs of data
# create sf object to get extent: 
LAAL_points = st_as_sf(LAAL_points_sufficient_data, coords = c("x","y"), remove = FALSE, crs=4326)
BFAL_points = st_as_sf(BFAL_points_sufficient_data, coords = c("x","y"), remove = FALSE, crs=4326)
# npac_base<- ptolemy::extract_gshhg(points, resolution = "h", epsg = 3349, buffer = 5000,
#                                    simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS
# 
# # for faster test plotting
# npac_base_res <- ptolemy::extract_gshhg(allBFAL.rast.sf, resolution = "l", epsg = NULL, buffer = 5000,
#                                         simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS

npac_base_i <- ptolemy::extract_gshhg(BFAL_points, resolution = "i", epsg = 3832, buffer = 0,
                                      simplify = FALSE) # using an sf object to set the appropriate boundary box and CRS


LAAL_points <- st_transform(LAAL_points, crs=3832)

season_names <- c(
  `1` = "Winter",
  `2` = "Spring",
  `3` = "Summer",
  `4` = "Fall"
)
figure_LAAL <- ggplot() + 
  # base map and other parameters
  geom_sf(data=LAAL_points, aes(color=island), size=0.005, alpha = 0.5)+
  geom_sf(data=npac_base_i, fill="grey60")+
  scale_color_manual(values=c("orange","firebrick"))+
  theme(legend.key.size = unit(3,"line"))+
  guides(color = guide_legend(override.aes = list(size=10)))+
  facet_wrap(~month_factor)
  #facet_wrap(~season, labeller = as_labeller(season_names))+
  #ggtitle("LAAL locations by month and breeding colony")
figure_LAAL

gp <- ggplotGrob(figure_LAAL)
gtable::gtable_show_layout(gp)



# now do again for BFAL 

BFAL_points <- st_transform(BFAL_points, crs=3832)

season_names <- c(
  `1` = "Winter",
  `2` = "Spring",
  `3` = "Summer",
  `4` = "Fall"
)
figure_BFAL <- ggplot() + 
  # base map and other parameters
  geom_sf(data=BFAL_points, aes(color=island), size=0.005, alpha = 0.5)+
  geom_sf(data=npac_base_i, fill="grey60")+
  scale_color_manual(values=c("turquoise2","royalblue4"))+
  theme(legend.key.size = unit(3,"line"))+
  guides(color = guide_legend(override.aes = list(size=10)))+
  facet_wrap(~month_factor)
  #facet_wrap(~season, labeller = as_labeller(season_names))+
  #ggtitle("BFAL locations by month and breeding colony")
figure_BFAL

library(ggpubr)
combined_figure <- ggarrange(figure_LAAL, figure_BFAL, 
                             labels = c("a", "b"),
                             ncol = 1, nrow = 2)
combined_figure


