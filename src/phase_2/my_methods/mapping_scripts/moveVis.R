# simple moveVis track for one of my birds for my Master's Seminar

library(moveVis)
library(move)
library(lubridate)

BFAL <- rbind(BFALmid, BFALtern)
"E:/project_data/spatial_segregation/data/LAALdata_midway_withTrackID.Rdata"
x <- read.csv("E:/project_data/spatial_segregation/data/midway_postbreeding_exports/LAAL/2011/pb_1111_01.csv")

?df2move

x$dtime <- strptime(x$dtime,"%Y-%m-%d %H:%M:%S")
x$dtime <- as.POSIXct(x$dtime)

convlon <- function(lon){
  ifelse(lon > 180, lon - 360, lon)
}

x$Longitude <- convlon(x$Longitude)

y <- df2move(x, proj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), x = 'Longitude', y = 'Latitude', time = 'dtime', track_id = NULL)

m <- align_move(y, res = 'min', unit = "hours")

#m <- sp::spTransform(m, sp::CRS("+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

frames <- frames_spatial(m, path_colours = c("red"),
                         map_service = "osm", map_type = "grayscale", alpha = 0.5, cross_dateline = FALSE) %>%
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_timestamps(m, type = "label") %>%
  add_gg(gg = expr(theme(aspect.ratio = 0.5)))
  add_progress()


frames[[100]] # preview one of the frames, e.g. the 100th frame

# animate frames
animate_frames(frames, out_file = "moveVis.gif")
