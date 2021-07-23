# implementing 'amt' to calculate HR overlap 
# March 5 2021 Dallas Jordan
# references: 
  # https://cran.r-project.org/web/packages/amt/vignettes/p2_hr.html#fn2
  # https://www.biorxiv.org/content/10.1101/2020.08.19.256859v2
  # you have several helpful resources in your workona Statistical Resources folder/tab


library(amt)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(sf)
library(rgdal)

      leroy <- amt_fisher %>% filter(name == "Leroy")
      lupe <- amt_fisher %>% filter(name == "Lupe")
      
      # get a CRS assigned
      amt_fisher<- transform_coords(amt_fisher, sp::CRS("+init=epsg:4326"))
      
      # trast = template raster
      trast <- make_trast(amt_fisher %>% filter(name %in% c("Leroy", "Lupe")), res = 50)
      trast <- make_trast(lupe, res = 50)
      hr_leroy <- hr_kde(leroy, trast = trast, levels = c(0.5, 0.9))
      hr_lupe <- hr_kde(lupe, trast = trast, levels = c(0.5, 0.9))
      

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_midway.Rdata")
LAALmid <- LAAL
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/LAALdata_tern.Rdata")
LAALtern <- LAAL

LAAL <- rbind(LAALmid, LAALtern)
LAAL$id <- "LAAL"

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_midway.Rdata")
BFALmid <- BFAL
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_tern.Rdata")
BFALtern <- BFAL

BFAL <- rbind(BFALmid, BFALtern)
BFAL$id <- "BFAL"

both <- rbind(LAAL,BFAL)
df1 <-tibble(both)

laea <- sp::CRS("+proj=laea +lon_0=180")
tr1<-make_track(df1,x,y,id=id, crs = laea)

LAALtrack <- tr1 %>% filter(id=="LAAL")
BFALtrack <- tr1 %>% filter(id=="BFAL")

trast <- make_trast(tr1 %>% filter(id %in% c("LAAL", "BFAL")), res = 300)

hr_LAAL <- hr_kde(LAALtrack, trast=trast, levels = 0.95)
hr_BFAL <- hr_kde(BFALtrack, trast=trast, levels = 0.95)

# hr and phr are directional, this means the order matters. For all other overlap measures the order does not matter.

test<-hr_overlap(hr_BFAL, hr_LAAL,type='udoi', conditional=T) 
test


######################
# BFAL between Midway and Tern

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_midway.Rdata")
BFALmid <- BFAL
load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/publication_figures/BFALdata_tern.Rdata")
BFALtern <- BFAL

BFALmid$id <- "BFALmid"
BFALtern$id <- "BFALtern"

both <- rbind(BFALmid,BFALtern)
df1 <-tibble(both)

tr1<-make_track(df1,x,y,id=id)

BFALmidtrack <- tr1 %>% filter(id=="BFALmid")
BFALterntrack <- tr1 %>% filter(id=="BFALtern")

trast <- make_trast(tr1 %>% filter(id %in% c("BFALmid", "BFALtern")), res = 300)

hr_BFALm <- hr_kde(BFALmidtrack,  levels = c(0.5, 0.95))
hr_BFALt <- hr_kde(BFALterntrack, levels = c(0.5, 0.95))

# hr and phr are directional, this means the order matters. For all other overlap measures the order does not matter.

hr_overlap(hr_BFALm, hr_BFALt,type='udoi') 








