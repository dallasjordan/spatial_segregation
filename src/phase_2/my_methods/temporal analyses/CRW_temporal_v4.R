# Revised CRW code: 
# determines distributions of data
# pulls step length and turning angle from fit distributions rather than non-parametrically
# Re run CRW on aggregate and in temporal chunks, spring, summer, and fall
# Thoughts to speed code up: 
# In CRW code, trim it to not print all the points? Do you need to do that? No, just one for an example of a figure.
# Dallas Jordan April 21 2022

library(dplyr)
library(lubridate)
library(reproj)
library(adehabitatLT)
library(fitdistrplus)
library(gamlss)
library(SuppDists)

library(geosphere)
library(beepr)
library(lubridate)
library(maptools)
library(reproj)
library(dplyr)
library(tidyr)
library(ggplot2)
library(adehabitatLT)
library(adehabitatHR)
library(sf)
library(reshape2)

# Step 1.1 - Fit distributions -----------------------------------------------

# Grab turning angles and step lengths

# USE RELATIVE PATHWAYS - setwd ONCE and leave it! 
setwd("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/")
# list Tern IS. BFAL FILES
d_spps_bfal_tern <- list.files(path = "./tern_postbreeding_exports/import_for_CRW/BFAL/", pattern = ".csv$", recursive=T, full.names = T)
d_spps_bfal_tern <- lapply(d_spps_bfal_tern, FUN = read.csv, header = T)

# list Tern IS. LAAL Files
d_spps_laal_tern <- list.files(path = "./tern_postbreeding_exports/import_for_CRW/LAAL/", pattern = ".csv$", recursive=T, full.names = T)
d_spps_laal_tern <- lapply(d_spps_laal_tern, FUN = read.csv, header=T)

d_spps_bfal_mid <- list.files(path = "./midway_postbreeding_exports/BFAL/",pattern = ".csv$", recursive=T, full.names = T)
d_spps_bfal_mid <- lapply(d_spps_bfal_mid, FUN = read.csv, header = T)
d_spps_laal_mid <- list.files(path = "./midway_postbreeding_exports/LAAL/", pattern = ".csv$", recursive=T, full.names = T)
d_spps_laal_mid <- lapply(d_spps_laal_mid, function(x) read.csv(x, header=T))

#BFAL TERN
d_spps_bfal_ternID <- c()
for(i in 1:length(d_spps_bfal_tern)){ # never hardcode for loop lengths, your code won't be applicable to any other data
  w <- d_spps_bfal_tern[[i]]
  w$ID <- paste("BFAL_TERN_", i, sep="")
  if(is.na(ymd(w$GMT_YYYY.MM.DD)) == TRUE){
    w$GMT <- mdy_hm(paste(w$GMT_YYYY.MM.DD, w$GMT_HH.mm.ss, sep=" "))
  }else{
    w$GMT <- ymd_hm(paste(w$GMT_YYYY.MM.DD, w$GMT_HH.mm.ss, sep=" "))
  }
  d_spps_bfal_ternID <- rbind(d_spps_bfal_ternID, w)
}
d_spps_bfal_ternID <- d_spps_bfal_ternID %>% dplyr::select(ID, GMT, Longitude, Latitude )

#LAAL TERN
d_spps_laal_ternID <- c()
for(i in 1:length(d_spps_laal_tern)){
  w <- d_spps_laal_tern[[i]]
  w$ID <- paste("LAAL_TERN_", i, sep="")
  if(is.na(ymd(w$GMT_YYYY.MM.DD)) == TRUE){
    w$GMT <- mdy_hm(paste(w$GMT_YYYY.MM.DD, w$GMT_HH.mm.ss, sep=" "))
  }else{
    w$GMT <- ymd_hm(paste(w$GMT_YYYY.MM.DD, w$GMT_HH.mm.ss, sep=" "))
  }
  d_spps_laal_ternID <- rbind(d_spps_laal_ternID, w)
}
d_spps_laal_ternID <- d_spps_laal_ternID %>% dplyr::select(ID, GMT, zm.lon, zm.lat )
names(d_spps_laal_ternID) <- names(d_spps_bfal_ternID)

#BFAL MID
d_spps_bfal_midID <- c()
for(i in 1:length(d_spps_bfal_mid)){
  w <- d_spps_bfal_mid[[i]]
  w$ID <- paste("BFAL_MID_", i, sep="")
  w$dtime <- ymd_hms(w$dtime)
  ee <- w %>% dplyr::select(ID, dtime, Longitude, Latitude)
  colnames(ee) <- colnames(d_spps_bfal_ternID)
  d_spps_bfal_midID <- rbind(d_spps_bfal_midID, ee)
}

#LAAL MID
d_spps_laal_midID <- c()
for(i in 1:length(d_spps_laal_mid)){
  w <- d_spps_laal_mid[[i]]
  w$ID <- paste("LAAL_MID_", i, sep="")
  w$dtime <- ymd_hms(w$dtime)
  ee <- w %>% dplyr::select(ID, dtime, Longitude, Latitude)
  colnames(ee) <- colnames(d_spps_bfal_ternID)
  d_spps_laal_midID <- rbind(d_spps_laal_midID, ee)
}

# make dataframe of all categories with IDs
d_birds <- rbind(d_spps_bfal_midID, d_spps_bfal_ternID, d_spps_laal_midID, d_spps_laal_ternID)

# reproject to make everything in PDC mercator, doing this because we need an equidistant projection
newLL <- reproj(cbind(d_birds$Longitude, d_birds$Latitude), target = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", source = "+proj=longlat +datum=WGS84 +no_defs")

# append projected lat/lons (in meters) to original big dataframe of all categories with IDs
d_birds$Longitude2 <- newLL[,1]
d_birds$Latitude2 <- newLL[,2]

# add column of ID and time per location
d_birds$match <- paste(d_birds$ID, d_birds$GMT, sep=".")

# Make traj object (necessary for step lengths and turning angles) - not necessary
# You can use a function from bayesmove "prep_data"
#d_birds$date <- d_birds$GMT
#bayesmove::prep_data(dat = d_birds, coord.names = c("Longitude2", "Latitude2"), id = "ID")
d_birds_traj <- as.ltraj(cbind(d_birds$Longitude2, d_birds$Latitude2), date=d_birds$GMT, id = d_birds$ID, typeII = T,
                         proj4string = CRS("+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
# convert to SPDF
d_birds_sp <- ltraj2spdf(d_birds_traj)

# convert to dataframe
d_birds_df <- data.frame(d_birds_sp)

# merge by time column in big dataframe of all categoeries with IDs
d_birds_df <- merge(d_birds_df, d_birds[,5:7], by.x="pkey", by.y="match")
colnames(d_birds_df)[13] <- "Longitude"
colnames(d_birds_df)[14] <- "Latitude"

# make new dataframe to fix names of key column
nms <- data.frame(do.call(rbind, strsplit(d_birds_df$pkey, split="_", fixed=T)))
nms2 <- data.frame(do.call(rbind, strsplit(nms$X3, split=".", fixed=T)))
nms <- data.frame(nms, nms2)
colnames(nms) <- c("SPP", "tag_loc", "ptID", "track", "trackpt")

# append to dataframe that was built from traj object
d_birds_df <- data.frame(d_birds_df, nms)
d_birds_df$Animal_ID <- paste(d_birds_df$SPP, d_birds_df$tag_loc, d_birds_df$track, sep="_")

# ALL ANGLES LAAL 

angles_LAAL <- d_birds_df$abs.angle[d_birds_df$SPP == "LAAL"]
angles_LAAL <- angles_LAAL[!is.na(angles_LAAL)]
angles_LAAL <- (angles_LAAL * 180) / pi # convert to degrees from radians

# PULL OUT SEASONAL ANGLES FOR LAAL

angles_LAAL_spring <- d_birds_df$abs.angle[d_birds_df$SPP == "LAAL" & (month(d_birds_df$date)==6 | month(d_birds_df$date)==7)]
angles_LAAL_spring <- angles_LAAL_spring[!is.na(angles_LAAL_spring)]
angles_LAAL_spring <- (angles_LAAL_spring * 180) / pi # convert to degrees from radians

angles_LAAL_summer <- d_birds_df$abs.angle[d_birds_df$SPP == "LAAL" & (month(d_birds_df$date)==8 | month(d_birds_df$date)==9)]
angles_LAAL_summer <- angles_LAAL_summer[!is.na(angles_LAAL_summer)]
angles_LAAL_summer <- (angles_LAAL_summer * 180) / pi # convert to degrees from radians

angles_LAAL_fall <- d_birds_df$abs.angle[d_birds_df$SPP == "LAAL" & (month(d_birds_df$date)==10 | month(d_birds_df$date)==11)]
angles_LAAL_fall <- angles_LAAL_fall[!is.na(angles_LAAL_fall)]
angles_LAAL_fall <- (angles_LAAL_fall * 180) / pi # convert to degrees from radians

# ALL ANGLES BFAL
angles_BFAL <- d_birds_df$abs.angle[d_birds_df$SPP == "BFAL"]
angles_BFAL <- angles_BFAL[!is.na(angles_BFAL)]
angles_BFAL <- (angles_BFAL * 180) / pi

# PULL OUT SEASONAL ANGLES FOR BFAL

angles_BFAL_spring <- d_birds_df$abs.angle[d_birds_df$SPP == "BFAL" & (month(d_birds_df$date)==6 | month(d_birds_df$date)==7)]
angles_BFAL_spring <- angles_BFAL_spring[!is.na(angles_BFAL_spring)]
angles_BFAL_spring <- (angles_BFAL_spring * 180) / pi # convert to degrees from radians

angles_BFAL_summer <- d_birds_df$abs.angle[d_birds_df$SPP == "BFAL" & (month(d_birds_df$date)==8 | month(d_birds_df$date)==9)]
angles_BFAL_summer <- angles_BFAL_summer[!is.na(angles_BFAL_summer)]
angles_BFAL_summer <- (angles_BFAL_summer * 180) / pi # convert to degrees from radians

angles_BFAL_fall <- d_birds_df$abs.angle[d_birds_df$SPP == "BFAL" & (month(d_birds_df$date)==10 | month(d_birds_df$date)==11)]
angles_BFAL_fall <- angles_BFAL_fall[!is.na(angles_BFAL_fall)]
angles_BFAL_fall <- (angles_BFAL_fall * 180) / pi # convert to degrees from radians

# ALL DISTANCES LAAL

dist_LAAL <- d_birds_df$dist[d_birds_df$SPP == "LAAL"]
dist_LAAL <- dist_LAAL[!is.na(dist_LAAL)]

# PULL OUT SEASONAL DISTANCES FOR LAAL

dist_LAAL_spring <- d_birds_df$dist[d_birds_df$SPP == "LAAL" & (month(d_birds_df$date)==6 | month(d_birds_df$date)==7)]
dist_LAAL_spring <- dist_LAAL_spring[!is.na(dist_LAAL_spring)]

dist_LAAL_summer <- d_birds_df$dist[d_birds_df$SPP == "LAAL" & (month(d_birds_df$date)==8 | month(d_birds_df$date)==9)]
dist_LAAL_summer <- dist_LAAL_summer[!is.na(dist_LAAL_summer)]

dist_LAAL_fall <- d_birds_df$dist[d_birds_df$SPP == "LAAL" & (month(d_birds_df$date)==10 | month(d_birds_df$date)==11)]
dist_LAAL_fall <- dist_LAAL_fall[!is.na(dist_LAAL_fall)]

# ALL DIST BFAL

dist_BFAL <- d_birds_df$dist[d_birds_df$SPP == "BFAL"]
dist_BFAL <- dist_BFAL[!is.na(dist_BFAL)]

# PULL OUT SEASONAL DISTANCES FOR BFAL

dist_BFAL_spring <- d_birds_df$dist[d_birds_df$SPP == "BFAL" & (month(d_birds_df$date)==6 | month(d_birds_df$date)==7)]
dist_BFAL_spring <- dist_BFAL_spring[!is.na(dist_BFAL_spring)]

dist_BFAL_summer <- d_birds_df$dist[d_birds_df$SPP == "BFAL" & (month(d_birds_df$date)==8 | month(d_birds_df$date)==9)]
dist_BFAL_summer <- dist_BFAL_summer[!is.na(dist_BFAL_summer)]

dist_BFAL_fall <- d_birds_df$dist[d_birds_df$SPP == "BFAL" & (month(d_birds_df$date)==10 | month(d_birds_df$date)==11)]
dist_BFAL_fall <- dist_BFAL_fall[!is.na(dist_BFAL_fall)]

comb <- d_birds_df

comb$month <- month(comb$date)
comb <- comb %>% mutate(season=case_when(
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
comb <- comb %>% mutate(season_name=case_when(
  season==1 ~ "Winter",
  season==2 ~ "Spring",
  season==3 ~ "Summer",
  season==4 ~ "Fall"
))
comb <- comb %>% mutate(month = month.name[month])

comb_spring <- comb %>% filter(month== "June" |  month=="July")
comb_summer <- comb %>% filter(month== "August" |  month=="September")
comb_fall <- comb %>% filter(month== "October" |  month=="November")


d_birds_df$month <- month(d_birds_df$date)
d_birds_df <- d_birds_df %>% mutate(season=case_when(
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
d_birds_df <- d_birds_df %>% mutate(season_name=case_when(
  season==1 ~ "Winter",
  season==2 ~ "Spring",
  season==3 ~ "Summer",
  season==4 ~ "Fall"
))
d_birds_df <- d_birds_df %>% mutate(month = month.name[month])

d_birds_df_spring <- d_birds_df %>% filter(month== "June" |  month=="July")
d_birds_df_summer <- d_birds_df %>% filter(month== "August" |  month=="September")
d_birds_df_fall <- d_birds_df %>% filter(month== "October" |  month=="November")

# you now have angle_XXX and dist_XXX, big dataframes of turning angles and distances by species. You are combining
# them 


# Step 1.2 Check distributions of combined and of seasonal angles/distances --------

# all LAAL distances
plotdist(dist_LAAL, histo=T,demp=T)
descdist(dist_LAAL, boot=1000)

fw <- fitdist(dist_LAAL, "weibull")
plot(fw)

fg <- fitdist(dist_LAAL/1000, "gamma")
plot(fg)

fln <- fitdist(dist_LAAL,"lnorm") # winner
plot(fln)
summary(fln)

fn <- fitdist(dist_LAAL,"norm")
plot(fn)

# seasonal LAAL distances
# spring
plotdist(dist_LAAL_spring, histo=T, demp=T)
descdist(dist_LAAL_spring, boot=1000)

fg <- fitdist(dist_LAAL_spring/1000, "gamma")
plot(fg)

fln <- fitdist(dist_LAAL_spring,"lnorm") # winner
plot(fln)
summary(fln)

# summer 
plotdist(dist_LAAL_summer, histo=T, demp=T)
descdist(dist_LAAL_summer, boot=1000)

fg <- fitdist(dist_LAAL_summer/1000, "gamma")
plot(fg)

fln <- fitdist(dist_LAAL_summer, "lnorm") # winner
plot(fln)

dist_LAAL_summer_normalized <- (dist_LAAL_summer-min(dist_LAAL_summer))/(max(dist_LAAL_summer)-min(dist_LAAL_summer))
dlsn.df <- data.frame(dist_LAAL_summer_normalized)
dlsn.df[5757,] <- .000001
dlsn.df[3383,] <- .999998
fb <- fitdist(dlsn.df$dist_LAAL_summer_normalized,"beta")
plot(fb)

fw <- fitdist(dist_LAAL_summer,"weibull")
plot(fw)

fe <- fitdist(dist_LAAL_summer/1000, "exp")
plot(fe)

fn <- fitdist(dist_LAAL_summer, "invGauss")

# neither beta, lnorm, nor weibull are great fits. Try fitting using gamlss
fit <- fitDist(dist_LAAL_summer_normalized, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)
summary(fit)
# suggests exponential Gaussian. This is a very uncommon dist, not going to use. 

# fall
plotdist(dist_LAAL_fall, histo=T, demp=T)
descdist(dist_LAAL_fall, boot=1000)

fln <- fitdist(dist_LAAL_fall,'lnorm') # winner
plot(fln)

fe <- fitdist(dist_LAAL_fall/1000,"exp")
plot(fe)

fw <- fitdist(dist_LAAL_fall,"weibull")
plot(fw)


# all BFAL distances
plotdist(dist_BFAL, histo=T,demp=T)
descdist(dist_BFAL, boot=1000)

fw <- fitdist(dist_BFAL, "weibull")
plot(fw)

fg <- fitdist(dist_BFAL/1000, "gamma")
plot(fg)

fln <- fitdist(dist_BFAL,"lnorm") # winner
plot(fln)
summary(fln)

fn <- fitdist(dist_BFAL,"norm")
plot(fn)


# par(mfrow = c(2, 2))
# plot.legend <- c("Weibull", "lognormal", "gamma")
# denscomp(list(fw, fln), legendtext = plot.legend)
# qqcomp(list(fw, fln), legendtext = plot.legend)
# cdfcomp(list(fw, fln), legendtext = plot.legend)
# ppcomp(list(fw, fln), legendtext = plot.legend)

# Seasonal BFAL distances
# Spring
plotdist(dist_BFAL_spring, histo=T, demp=T)
descdist(dist_BFAL_spring, boot=1000)

fw <- fitdist(dist_BFAL_spring, "weibull")
plot(fw)

fg <- fitdist(dist_BFAL_spring/1000, "gamma")
plot(fg)

fln <- fitdist(dist_BFAL_spring,"lnorm") # winner but not great
plot(fln)

fn <- fitdist(dist_BFAL_spring,"norm")
plot(fn)

# Summer
plotdist(dist_BFAL_summer, histo=T, demp=T)
descdist(dist_BFAL_summer, boot=1000)

fw <- fitdist(dist_BFAL_summer, "weibull")
plot(fw)

fg <- fitdist(dist_BFAL_summer/1000, "gamma")
plot(fg)

fln <- fitdist(dist_BFAL_summer,"lnorm") # winner
plot(fln)

fn <- fitdist(dist_BFAL_summer,"norm")
plot(fn)

# Fall
plotdist(dist_BFAL_fall, histo=T, demp=T)
descdist(dist_BFAL_fall, boot=1000)

fw <- fitdist(dist_BFAL_fall, "weibull")
plot(fw)

fg <- fitdist(dist_BFAL_fall/1000, "gamma")
plot(fg)

fln <- fitdist(dist_BFAL_fall,"lnorm") # winner
plot(fln)

fn <- fitdist(dist_BFAL_fall,"norm")
plot(fn)


# LAAL turning angles

# All
plotdist(angles_LAAL,histo=T,demp=T)
descdist(angles_LAAL,boot=1000)
fu <- fitdist(angles_LAAL_spring,"unif") # winner
plot(fu)

# Spring  
plotdist(angles_LAAL_spring, histo=T,demp=T)
descdist(angles_LAAL_spring, boot=1000)
fu <- fitdist(angles_LAAL_spring,"unif") # winner
plot(fu)

# Summer 
plotdist(angles_LAAL_summer, histo=T,demp=T)
descdist(angles_LAAL_summer, boot=1000)
fu <- fitdist(angles_LAAL_summer,"unif") # winner
plot(fu)

# Fall 
plotdist(angles_LAAL_fall, histo=T,demp=T)
descdist(angles_LAAL_fall, boot=1000)
fu <- fitdist(angles_LAAL_fall,"unif") # winner
plot(fu)

# BFAL turning angles

# All
plotdist(angles_BFAL,histo=T,demp=T)
descdist(angles_BFAL,boot=1000)
fu <- fitdist(angles_BFAL,"unif") # winner
plot(fu)

# Spring  
plotdist(angles_BFAL_spring, histo=T,demp=T)
descdist(angles_BFAL_spring, boot=1000)
fu <- fitdist(angles_BFAL_spring,"unif") # winner
plot(fu)

# Summer 
plotdist(angles_BFAL_summer, histo=T,demp=T)
descdist(angles_BFAL_summer, boot=1000)
fu <- fitdist(angles_BFAL_summer,"unif") # winner
plot(fu)

# Fall 
plotdist(angles_BFAL_fall, histo=T,demp=T)
descdist(angles_BFAL_fall, boot=1000)
fu <- fitdist(angles_BFAL_fall,"unif") # winner
plot(fu)

# Step 1.3 SUMMARY OF DISTRIBUTION RESULTS 4/15 ------------------------------------

# All distances were lognormal, whether aggregated or by season. There were a couple
# that were not great fits, but out of the typical dispersal distance distributions
# (lognormal, gamma, weibull, negative exponential; see Limpert et al. 2001 and Diefenbach et al. 2008)
# lognormal was the best fit. gamlss package suggested inverse Gaussian might be a good fit, 
# but I couldn't figure out how to work with invGauss in R. 

# All angles were uniform 0-359. Pretty straightforward. I thought there might be seasonal differences but 
# the fits were all close to uniform. 

# Step 2 - CRW drawn from distributions --------------------------------------------
# Going to conduct CRW with distributions that I determined. When I draw from 
# the lognormal distribution, I have to specify meanlog and sdlog. What values do I pick?
# Seasonal and aggregated LAAL/BFAL all are ~11.6 and ~.9. Maybe just pick the aggregated species values?

### pull first row of each, simulate from there

# IMPORTANT - CHANGE SEASON HERE FOR EACH TIME YOU RUN BY SEASON

d_birds_df <- d_birds_df_summer
comb <- comb_summer

first_pts <- NULL
for(i in unique(comb$Animal_ID)){
  p <- comb[comb$Animal_ID == i,]
  p <- p[which(p$trackpt == min(p$trackpt)),]
  first_pts <- rbind(first_pts, p)
}

setwd("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/overlap_sensitivity/")

# CRW STARTS HERE
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

# setup cluster in single computer world
# detect and define n.cores
n.cores <- parallel::detectCores()-1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()
## [1] TRUE
#how many workers are available? (optional)
foreach::getDoParWorkers()

# CREATE FUNCTION FOR CREATING AND TESTING NEW RANDOMLY GENERATED POINTS 
CRW_FUN_PAR <- function(temp_FP, land_mask, num_tracks, lenndf, x, IDs){
  
  bloopL <- list() # create empty list to store bloops 1:num_tracks
  temp_FP$new_dis <- 0
  temp_FP$new_ang <- 0
  #temp_FP$date <- ymd_hms(temp_FP$date)
  temp_FP$num <- 0
  temp_FP$iter <- 0
  
  for(n in 1:num_tracks){
    draws <- lenndf$nobs[n] # get this tracks length
    bloop <- temp_FP %>% filter(Animal_ID==IDs[n]) %>%
      dplyr::select(Animal_ID, Longitude, Latitude,#date, 
                    new_dis, new_ang, num)
    
    # set the new track sim counter to zero
    lenn_new_sim_df = 0
    stop=1
    lenn <-lenndf$nobs[n]
    while (lenn_new_sim_df<(lenn)){
      new_x <- c()
      new_y <- c()
      #timee <- .POSIXct(character())
      # DRAW RANDOM VARIALBES
      new_dis <- rlnorm(draws,meanlog = 11.5943286, sdlog=0.9036841) # values pulled from fitdistr for relevant time period
      new_ang <- runif(draws,min=0,max=359) # uniform
      new_x[1] <- bloop$Longitude[stop]
      new_y[1] <- bloop$Latitude[stop]
      #timee[1] <- ymd_hms(bloop$date[stop])
      if (draws <2){
        new_x[2] <- new_x[1] + new_dis[1]*cos(new_ang[1])
        new_y[2] <- new_y[1] + new_dis[1]*sin(new_ang[1])
        #timee[2] <- ymd_hms(timee[1]) + 5
      }else{
        for(d in 2:draws){
          new_x[d] <- new_x[d-1] + new_dis[d]*cos(new_ang[d])
          new_y[d] <- new_y[d-1] + new_dis[d]*sin(new_ang[d])
          #timee[d] <- ymd_hms(timee[d-1]) + 5
        }
      }
      
      newptsdf <- data.frame("Lon"=new_x, "Lat"=new_y)
      
      testpoints <- st_as_sf(x = newptsdf, coords = c("Lon", "Lat")) 
      testpoints <- testpoints%>% st_set_crs(st_crs(land_mask))
      
      in_land <- st_join(testpoints, land_mask, join = st_intersects)
      
      if (all(is.na(in_land$type)) & all(-1708182<new_x) & all(new_x<9844089) & all(1185132<new_y) & all(new_y<9310913)){
        e <- data.frame(Longitude = new_x, Latitude = new_y, 
                        #date = timee, 
                        new_dis = new_dis, new_ang = new_ang, num = 1:length(new_ang))
        e$Animal_ID <- IDs[n]
        bloop <- rbind(bloop, e)
        lenn_new_sim_df <- nrow(bloop)
        
      }else{
        # find the index just prior to the first fail point
        stop <- min(c(which((is.na(in_land$type))==FALSE), which((-1708182<new_x)==F), which((new_x<9844089)==F), which((1185132<new_y)==F), which((new_y<9310913)==F)))-1
        e <- data.frame(Longitude = new_x[1:stop], Latitude = new_y[1:stop], 
                        #date = timee[1:stop], 
                        new_dis = new_dis[1:stop], new_ang = new_ang[1:stop], num = 1:stop)
        e$Animal_ID <- IDs[n]
        bloop <- rbind(bloop,e)
        lenn_new_sim_df <- nrow(bloop)-1
        draws = lenn - lenn_new_sim_df
      } # end of elseif
    } # end of while
    
    bloop$track_num <- n # this is effectively assigning the sim number of 1:num_tracks
    bloop <- bloop[!(bloop$num==0),] # remove initial point which for some reason duplicates in some iterations
    bloopL[[n]] <- bloop
  } # end of num_tracks for loop
  return(bind_rows(bloopL))
} # end of function


# set up landmask
library(sf)
sf::sf_use_s2(FALSE)
land_mask_bbox<-st_bbox(c(xmin = -3008182, xmax = 10844089, ymax = 18027535, ymin = 185132), crs = st_crs(3832))
land_mask_sfc <- st_as_sfc(land_mask_bbox)
land_mask <- ptolemy::extract_gshhg(data = land_mask_sfc, resolution = "i", epsg = NULL, buffer = 5000,
                                    simplify = FALSE) 
land_mask$type = "land"
simtype <- "birds"

# define spp and colony
spp <- sapply(strsplit(first_pts$Animal_ID, "_"), "[[", 1)
col <- sapply(strsplit(first_pts$Animal_ID, "_"), "[[", 2)
sppcol <- paste(col, spp, sep = "_")
first_pts$sppcol <- sppcol

d_birds_df$spp <- sapply(strsplit(d_birds_df$Animal_ID, "_"), "[[", 1)
d_birds_df$col <- sapply(strsplit(d_birds_df$Animal_ID, "_"), "[[", 2)
d_birds_df$sppcol <- paste(d_birds_df$col, d_birds_df$spp, sep = "_")

# BEGIN CRW
start <- Sys.time()
save_spp_col_list <- list()
for(sc in unique(sppcol)){
  #get dataframe of first track points for each animal ID in species_colony[sc]
  temp_FP <- first_pts[first_pts$sppcol == sc,]
  tmp_d_birds_df <- d_birds_df[d_birds_df$sppcol == sc,]
  
  # Calculate the number of tracks for this species / colony combo
  num_tracks <- length(unique(temp_FP$Animal_ID))
  # calculate length of tracks of all animal IDs 
  lenndf <- tmp_d_birds_df %>% group_by(Animal_ID) %>% summarise(nobs = n())
    print(sc)
  # THIS IS WHERE WE MUST APPLY FUNCTION - create
  my.seed <- 76513
  set.seed( my.seed, kind = "L'Ecuyer-CMRG" )
  sim_CRW_pts_L <- parallel::mclapply(X=c(1:1), FUN = CRW_FUN_PAR, temp_FP = temp_FP, land_mask = land_mask, num_tracks = num_tracks, lenndf = lenndf, IDs = unique(temp_FP$Animal_ID))
  save_spp_col_list <- append(save_spp_col_list,sim_CRW_pts_L) # SAVE THIS OUT
} # end of species / colony loop

save(save_spp_col_list,file="CRW_v4_augustseptember.Rdata")
end <- Sys.time()

test_output <- save_spp_col_list %>% bind_rows() %>% group_by(Animal_ID) %>% summarize(nobs=n())
test_original <-d_birds%>%group_by(ID)%>%summarize(nobs=n())
