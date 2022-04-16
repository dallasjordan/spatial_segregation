library(sf)
library(tidyverse)

load("~/Downloads/vert95_allLAAL.Rdata")

class(vert95_allLAAL)

eqc <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

al95c <- sf::st_as_sf(vert95_allLAAL) %>% 
  st_transform(4326) %>%
  st_shift_longitude() %>%
  st_transform(eqc)


plot(al95c)

load("~/Downloads/npac_base_res.Rdata")

figure1 <- ggplot() + 
  # LAAL raster
  geom_sf(data=al95c, color=("#E39191"), size=1, fill=alpha("#E39191",0.2)) +
  geom_sf(data=npac_base_res, fill='grey60') +
  theme_bw() +
  coord_sf(expand=F)

figure1
