# Lotek twilight times with SGAT
# Michael Sumner

# 2014-09-01

# Combining Michael Sumner's codes: ternisland03.R and lotek_riseset.R into one that pulls in lotek ldata and BAS prs data and then does
# the same estimation

# 2015-07-23 - Melinda Conners - run with modifications


library(SGAT)
library(BAStag)
library(chron)
library(pracma)
library(raster)
library(maptools)
library(rgeos)
library(maptools)


# Pre Define Stuff :
##-----------------------------------------------------------------------------

# Record some details that we know about this tag
lon.home <- -166.29 + 360
lat.home <- 23.87

## maximum range in lon/lat (for land mask plus maximum bound of track)
lon.range <- c(140, 245)
lat.range <- c(23, 62)

# Land mask data

data(wrld_simpl)

land.mask <- function(xlim, ylim, n = 4, land = TRUE) {
  r <- raster(nrows = n*diff(ylim), ncols = n*diff(xlim),
              xmn=xlim[1],xmx=xlim[2],
              ymn=ylim[1],ymx=ylim[2],
              crs=proj4string(wrld_simpl))
  r <- cover(rasterize(elide(wrld_simpl,shift=c(-360,0)),r,1,silent=TRUE),
             rasterize(wrld_simpl,r,1,silent=TRUE),
             rasterize(elide(wrld_simpl,shift=c(360,0)),r,1,silent=TRUE))
  r <- as.matrix(is.na(r))[nrow(r):1,]
  if(land) r <- !r
  xbin <- seq(xlim[1],xlim[2],length=ncol(r)+1)
  ybin <- seq(ylim[1],ylim[2],length=nrow(r)+1)
  
  function(p) {
    r[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
  }
}


is.sea <- land.mask(xlim = lon.range, ylim = lat.range, n = 4, land = FALSE)


log.prior <- function(p)  {
  f <- is.sea(p)
  ifelse(f & !is.na(f),0,-1000)
}


nlocation <- function(s) {
  dim(if(is.list(s)) s[[1]] else s)[1]
}

##-----------------------------------------------------------------------------------
# Add Data of Dates when birds are at colony

## filenames of ColonyDate Files (files named by bird)

colony_dates<-list.files("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/sst_fix_ugh/ugh/")
colony_dates_sub<- substr(colony_dates,1,16)

## filenames of start and end index
bookend<-list.files("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/FinalEst_Run01_07272015/fit_sst_truncated/startfix/locations") 
bookend.sub<-substr(bookend,1,16)
#    

######################################################################################
##-------------- BAS -----------------------------------------------------------------
######################################################################################

dp.bas <- "/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/prs_rename/"

setwd(dp.bas)
filenames<- list.files()
zenith = 93.5


i = i+1 # This isn't in a loop, because it requires manual intervention for the twilight fixes

load(filenames[i])

mgcbird<-substr(filenames[i],1,16)
year<-substr(mgcbird,5,8)

# Location estimation- naive threshold-------------------------------------------------
path1 <- threshold.path(prs$Twilight, prs$Rise, zenith = zenith, unfold = TRUE)


######################################################################################
##-------------- LAT -----------------------------------------------------------------
######################################################################################

# dp.lotek <- "/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/ldata_rename/"
# 
# setwd(dp.lotek)
# filenames<- list.files()
# zenith = 93.5
# 
# 
# i = i+1 # This isn't in a loop, because it requires manual intervention for the twilight fixes
# 
# load(filenames[i])
# 
# mgcbird<-substr(filenames[i],1,16)
# year<-substr(mgcbird,5,8)
# 
# # Location estimation- naive threshold-------------------------------------------------
# twl <- data.frame(Rise = as.POSIXct(ldata$RiseMinute,tz="HST"), 
#                   Set = as.POSIXct(ldata$SetMinute,tz="HST"))
# 
# # interleave rise and set
# interleave <- function(x, y) {
#   .POSIXct(as.vector(t(cbind(x, y))), tz = "GMT")
# }
# 
# twl1 <- data.frame(Twilight = interleave(twl$Rise, twl$Set), Rise = c(TRUE, FALSE))
# 
# #####################################
# twl1<-twl1[45:539,] # Check Twilights
# #####################################
# 
# path1 <- threshold.path(twl1$Twilight, twl1$Rise, zenith = zenith, unfold=FALSE)






##----- equinox fix------------------------------------------------------------

# this clunky tool allows you to interpolate latitudes between the start and end of equinox drift.

fake<- cbind(1:length(path1$x[,2]),path1$time,path1$x[,2])
plot(fake[,1],fake[,3]) # plot latitude timeseries

latcoor<-locator() # identify start and end points of equinox within timeseries. Below allows for up to three equinox chunks. Can add more for longer trips.

lat1  <- latcoor$y[1:2]
lix1  <- round(latcoor$x[1:2])
lat2  <- latcoor$y[3:4]
lix2  <- round(latcoor$x[3:4])
lat3  <- latcoor$y[5:6]
lix3  <- round(latcoor$x[5:6])

lat.index1<-seq(lix1[1],lix1[2],by=1)
lat.index2<-seq(lix2[1],lix2[2],by=1)
lat.index3<-seq(lix3[1],lix3[2],by=1)

latfix1<-approx(lix1,lat1,method="linear",xout=lat.index1) #linear interpolation between
latfix2<-approx(lix2,lat2,method="linear",xout=lat.index2)
latfix3<-approx(lix3,lat3,method="linear",xout=lat.index3)

fake[lat.index1,3]<-latfix1$y
fake[lat.index2,3]<-latfix2$y
fake[lat.index3,3]<-latfix3$y

plot(fake[,1],fake[,3]) # Does this look better?


# put equinox fixes into path1
path1$x[,2]<-fake[,3]

# The estelle.metropolis model uses light data in each iteration so even though you fix these lats here, the model will continue to add equinox drift. 
# So these lats are used in the initial location estimation, from which then you can identify a max latitude to constrain the model (in the automated MCMC estimations part of the code)


##-- COLONY FIX LOCATIONS ---------------------------------------------------------------

# bookend is a 3 col matrix of : bird, startday, endday (start and end dates of post breeding trip)

# colony dates


##-----------------------------------------------------------------------------------
## ---- Find Start and End index (times)

  ix2<-which(bookend.sub==mgcbird)
  
  load(paste("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/FinalEst_Run01_07272015/fit_sst_truncated/startfix/locations/",bookend[ix2],sep=""))

  startday<-time.trunc[1]
  endday<-time.trunc[length(time.trunc)]


## ---- Find Other Colony Times

    colony_dates<-list.files("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/sst_fix_ugh/ugh/")
    colony_dates_sub<- substr(colony_dates,1,16)

    ix3<-which(colony_dates_sub==mgcbird)
    alldates<-read.csv(paste("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/sst_fix_ugh/ugh/",colony_dates[ix3],sep=""))
  
    if (year=="2011" | year =="2012"){
      dates<-as.POSIXct(as.character(alldates$x),"%d/%m/%Y",tz="GMT") # BAS data format
    }else{
      dates<-as.POSIXct(as.character(alldates$x),"%m/%d/%Y",tz="GMT")} # LAT data format


   # ---- Delete any Colony Times that are after Startx and before Endx
   colonyfixdates<-as.POSIXct(c(as.character(dates[which(dates<startday)]),substr(as.character(startday),1,10),substr(as.character(endday),1,10)),"%Y-%m-%d",tz="GMT")


## ---- find indices for colony positions:

#closest :

    time<-path1$time
    colix<-vector()
    for (j in 1:length(colonyfixdates)){colix[j]<-which(abs(time-colonyfixdates[j])==min(abs(time-colonyfixdates[j])))}
    colx<-unique(colix)

##------------ Build log priors and fixed locations -----------------------------------------x0 <- path1$x

    # Check initial locations
    x0 <- path1$x
    x0[x0[,1] < 0, 1] <- x0[x0[,1] < 0, 1] + 360
    plot(x0)
    

    ### #Delete then na approx locations that are on land mask
    bad <- log.prior(x0) < 0
    x0[bad, 1] <- NA
    x0[bad, 2] <- NA
    x0 <- zoo::na.approx(x0)
    plot(x0)

    # Check start / end NB trip colony locations
    abline(h = lat.home, v = lon.home)
    points(x0[colx[1],,drop = FALSE], col = "dodgerblue", pch = 16)  # check colony location start
    points(x0[nrow(colx),,drop = FALSE], col = "firebrick", pch = 16)# check colony location end
    

    ## now fix our known colony locations
    fixedx <- rep(FALSE, nlocation(x0))
    fixedx[colx] <- TRUE # fix colony locations
    x0[fixedx,1] <- lon.home
    x0[fixedx,2] <- lat.home
    z0 <- trackMidpts(x0)

    # Plot - - - Looking reasonable? If not ... use the "na zapper" (see code chunk below)
    data(wrld_simpl)
    plot(x0)
    points(z0,pch = 16, cex = 2, col = rgb(30L, 144L, 255L, 80L, max = 255))
    lines(x0)
    plot(wrld_simpl, add = TRUE)
    plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)
    
    ## --- If have issue with some on land locations ... use "na zapper"
    nacoors<-identify(x0)
    x0[nacoors, 1] <- NA
    x0[nacoors, 2] <- NA
    x0 <- na.approx(x0)
    z0 <- trackMidpts(x0)

    plot(x0)
    points(z0,pch = 16, cex = 2, col = rgb(30L, 144L, 255L, 80L, max = 255))
    lines(x0)
    plot(wrld_simpl, add = TRUE)
    plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)


##-----------------------------------------------------------------
# Model

# To see the Light Parameters (for alpha)
m <- seq(0,25,length=80)
plot(m,dlnorm(m,1.2,0.8),type="l",xlab="Twilight Error (minutes)")


# And the Speed parameters (for beta) # I used the same beta parameters and alpha parameters for all birds.
km <- trackDist(x0)
hr <- diff(as.numeric(twl1$Twilight))/3600
spd <- km/hr
mean(spd)
# hist(spd,15,freq=FALSE)
# lines(0:80,dgamma(0:80,2.27,0.13))
parameters.gamma(mean(spd[spd<20]),sd(spd[spd<20])) 
##  shape   rate 
##   2.2    0.35



# Now we build the model with these components, and run it.

# Both Estelle and Stella variants of the model assume that the speed of travel between successive (x) locations is gamma distributed with shape beta[1] and rate beta[2]. 

model <- threshold.model(twl1$Twilight, twl1$Rise,
                         twilight.model = "ModifiedLogNormal",
                         alpha = c(1.2,0.8), beta=c(2.2, 0.35),
                         logp.x=log.prior,logp.z=log.prior,
                         x0 = x0, z0 = z0, zenith = zenith, fixedx = fixedx)
   

x.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
z.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))

# Run initial fit !
fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=100,thin=20,chains=1)

# Check initial fit
xm <- location.mean(fit$x)
zm <- location.mean(fit$z)
plot(zm, type = "l")
points(xm, pch = "+")
points(zm, col = rgb(30L, 144L, 255L, 80L, max = 255))
plot(wrld_simpl, add = TRUE)
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)

# save initial fit
dp3<-"/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/"
save(fit,file=paste(dp3,"fit01_2015_",filenames[i],sep=""))

rm(list=ls()[! ls() %in% c("lon.home","lat.home","filenames","i","is.sea","land.mask","log.prior","nlocation","colony_dates","colony_dates_sub","bookend","bookend.sub","zenith")])







