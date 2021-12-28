# Area calculations of KDEs by species and colony
# You have an old version of this script that calculated by year as well, but the scripts are quite clunky and you can do it 
# in way fewer lines. Just re-writing it here: 
# Dec 3 2021
# Dallas Jordan


# Setup -------------------------------------------------------------------

library(adehabitatHR)
library(dplyr)
library(ggplot2)

setwd("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/final_push/final_ud/")
path <- getwd()
file_list <- list.files(path=path)

for (i in 1:length(file_list)){
  load(file_list[i])
}

# should have allBFAL/allLAAL/midBFAL/midLAAL/ternBFAL/ternLAAL data objects


# Calculate vertices, get area --------------------------------------------

midLAAL.poly95 <- getverticeshr(midLAAL, percent = 95, unin = "m", unout = "km2") 
midLAAL.poly95$id
midLAAL.poly95$area

midLAAL.poly50 <- getverticeshr(midLAAL, percent = 50, unin = "m", unout = "km2") 
midLAAL.poly50$id
midLAAL.poly50$area

midBFAL.poly95 <- getverticeshr(midBFAL, percent = 95, unin = "m", unout = "km2") 
midBFAL.poly95$id
midBFAL.poly95$area

midBFAL.poly50 <- getverticeshr(midBFAL, percent = 50, unin = "m", unout = "km2") 
midBFAL.poly50$id
midBFAL.poly50$area

ternLAAL.poly95 <- getverticeshr(ternLAAL, percent = 95, unin = "m", unout = "km2") 
ternLAAL.poly95$id
ternLAAL.poly95$area

ternLAAL.poly50 <- getverticeshr(ternLAAL, percent = 50, unin = "m", unout = "km2") 
ternLAAL.poly50$id
ternLAAL.poly50$area

ternBFAL.poly95 <- getverticeshr(ternBFAL, percent = 95, unin = "m", unout = "km2") 
ternBFAL.poly95$id
ternBFAL.poly95$area

ternBFAL.poly50 <- getverticeshr(ternBFAL, percent = 50, unin = "m", unout = "km2") 
ternBFAL.poly50$id
ternBFAL.poly50$area

