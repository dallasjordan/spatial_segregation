# simple moveVis track for one of my birds for my Master's Seminar

library(moveVis)
library(move)
library(lubridate)
library(ggplot2)
library(dplyr)

x <- read.csv("E:/project_data/spatial_segregation/data/midway_postbreeding_exports/LAAL/2011/pb_1111_01.csv")
x$name = "Laysan Albatross"

x$dtime <- strptime(x$dtime,"%Y-%m-%d %H:%M:%S")
x$dtime <- as.POSIXct(x$dtime)

convlon <- function(lon){
  ifelse(lon > 180, lon - 360, lon)
}

x$Longitude <- convlon(x$Longitude)

y <- df2move(x, proj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), x = 'Longitude', y = 'Latitude', time = 'dtime', track_id = NULL)

m <- align_move(y, res = 'min', unit = "hours")

m <- sp::spTransform(m, sp::CRS("+proj=cea +lon_0=150 +x_0=0 +y_0=0 +R=6371228 +units=m +no_defs +type=crs"))

frames <- frames_spatial(m, path_colours = c("red"),
                         map_service = "osm", map_type = "terrain", alpha = 0.5, cross_dateline = TRUE) %>%
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_timestamps(m, type = "label") %>%
  add_gg(gg = expr(theme(aspect.ratio = 1))) %>%
  add_gg(gg = expr(geom_point(aes(x=463347,y=3120213), shape=19)))


frames[[100]] # preview one of the frames, e.g. the 100th frame

# animate frames
animate_frames(frames, out_file = "moveVis.gif")







library(moveVis)
library(move)
library(lubridate)
library(ggplot2)
library(dplyr)
library(sf)

BFAL <- rbind(BFALmid, BFALtern)
"E:/project_data/spatial_segregation/data/LAALdata_midway_withTrackID.Rdata"



x <- read.csv("E:/project_data/spatial_segregation/data/midway_postbreeding_exports/LAAL/2009/pb_1109_03.csv")
x$name = "Laysan Albatross"

x$dtime <- strptime(x$dtime,"%Y-%m-%d %H:%M:%S")
x$dtime <- as.POSIXct(x$dtime)

convlon <- function(lon){
  ifelse(lon > 180, lon - 360, lon)
}

x$Longitude <- convlon(x$Longitude)

y <- df2move(x, proj = "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs", x = "Longitude", y = "Latitude", time = "dtime", track_id = "name")
m <- align_move(y, res = 'max', unit = "hours")

stopifnot(st_crs(m)==st_crs(4326))
frames <- frames_spatial(m, cross_dateline = TRUE,
                         trace_show = TRUE, map_service = "osm", map_type = "terrain") %>%
          add_gg(gg = expr(theme(aspect.ratio = 0.3)))%>%
          add_gg(gg = expr(ggplot2::coord_sf(crs = sf::st_crs(3975), datum = sf::st_crs(3975), clip = "on", expand = F)))# stretching the y axis a biteme(aspect.ratio = 0.1)))




  add_gg(gg = expr(ggplot2::coord_sf(crs = sf::st_crs(4326), datum = sf::st_crs(4326), clip = "on", expand = F)))
  

  add_gg(gg = expr(theme(aspect.ratio = 0.3))) # stretching the y axis a biteme(aspect.ratio = 0.1)))
  
  
  add_gg(gg = expr(ggplot2::coord_sf(crs = sf::st_crs(3975), datum = sf::st_crs(3975), clip = "on", expand = F))) %>%
  add_gg(gg = expr(theme(panel.background = element_rect(fill = NA), panel.ontop = TRUE)))


frames[[1]] # preview one of the frames, e.g. the 100th frame

# animate frames
animate_frames(frames, out_file = "moveVis2.gif")
