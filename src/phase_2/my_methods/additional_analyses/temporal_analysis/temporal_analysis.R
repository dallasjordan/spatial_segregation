# Temporal analysis
# Dallas Jordan
# March 17 2022

# Purpose of this script is to conduct my overlap analysis (calculate UDOI/BA/PRS metrics) for overall LAAL x BFAL 
# and the 4 island x spp combinations on a monthly basis. The results will go in a very large table.


# Setup -------------------------------------------------------------------

library(dplyr)
library(sf)
library(adehabitatHR)


lcea <- "+proj=cea +lat_0=0 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

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


counts_LAAL <- LAAL_points_sufficient_data %>% count(island,month_factor, id)
counts_LAAL
write.csv(counts_LAAL,"/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/counts_LAAL.csv")

counts_BFAL <- BFAL_points_sufficient_data %>% count(island,month_factor)
counts_BFAL
write.csv(counts_BFAL,"/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/spatial_segregation_git/counts_BFAL.csv")
