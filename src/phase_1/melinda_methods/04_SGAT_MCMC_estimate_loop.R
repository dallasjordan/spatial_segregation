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

# Record some details that we know about this tag, the deployment date will be used to ignore data prior to this time.
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


## ---------------------------------------------------------------------------------------------------------
dp.fit<-"/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fit01/"

setwd(dp.fit)

filenames.fit<- list.files()

for (i in 1:length(filenames.fit)) {
  
  load(filenames.fit[i])
  
  fixedx<-fit$model$fixedx
  x0<-fit$model$x0
  z0<-fit$model$z0
  model<-fit$model
  
  
  x.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
  z.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))
  
  fit <- estelle.metropolis(model,x.proposal,z.proposal,
                            iters=300,thin=20,chains=2)
  
  data(wrld_simpl)
  plot(chain.last(chain.collapse(fit$x)))
  lines(chain.last(chain.collapse(fit$x)))
  plot(wrld_simpl, add = TRUE)
  plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)
  
  
  # need to set log prior of maximum latitude
  xm <- apply(chain.collapse(fit$x), 1:2, mean)
  zm <- apply(chain.collapse(fit$z), 1:2, mean)
  max.lat<-max(xm[,2])
  min.lat<-min(xm[,2])
  min.lon <- 140
  max.lon <- 238
  
  data(wrld_simpl)
  is.sea <- land.mask(xlim = c(min.lon, max.lon), ylim = c(min.lat, max.lat), n = 4, land = FALSE) # range limit is adjusted for each bird.
  
  ## Define the log prior for x and z
  log.prior <- function(p)  {
    f <- is.sea(p)
    
    ifelse(f & !is.na(f),0,-1000)
  }
  
  
  
  
  for(k in 1:3) {
    x.proposal <- mvnorm(chain.cov(fit$x),s=0.25)
    z.proposal <- mvnorm(chain.cov(fit$z),s=0.25)
    fit <- estelle.metropolis(model,x.proposal,z.proposal,
                              x0=chain.last(fit$x),
                              z0=chain.last(fit$z),
                              iters=300,thin=20,chains=2)
    
    
    ## hope to see no obvious patterns in these chains
    ## once converged the values of each are constant within a range
    #     opar <- par(mfrow=c(length(fit$x),2),mar=c(3,5,2,1)+0.1)
    #     for(k in 1:length(fit$x)) {
    # #       matplot(t(fit$x[[k]][!fixedx,1,]),type="l",lty=1,col="dodgerblue",ylab="Lon")
    # #       matplot(t(fit$x[[k]][!fixedx,2,]),type="l",lty=1,col="firebrick",ylab="Lat")
    #     }
    #     par(opar) 
  }
  
  xm <- apply(chain.collapse(fit$x), 1:2, mean)
  zm <- apply(chain.collapse(fit$z), 1:2, mean)
  plot(zm, type = "l")
  points(xm, pch = "+")
  points(zm, col = rgb(30L, 144L, 255L, 80L, max = 255))
  plot(wrld_simpl, add = TRUE)
  plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)
  
  
  x.proposal <- mvnorm(chain.cov(fit$x),s=0.25)
  z.proposal <- mvnorm(chain.cov(fit$z),s=0.25)
  fit <- estelle.metropolis(model,x.proposal,z.proposal,
                            x0=chain.last(fit$x),
                            z0=chain.last(fit$z),
                            iters=1000,thin=20,chains=2)
  
  
  
  # save image of plot
  bird<-strsplit(filenames.fit[i],"\\.")[[1]][1]
  
  xm <- apply(chain.collapse(fit$x), 1:2, mean)
  zm <- apply(chain.collapse(fit$z), 1:2, mean)
  
  cairo_pdf(file=paste("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fitPDF/",bird,".pdf",sep=""), height=6.95, width=8.25,pointsize=12)
  plot(zm, type = "l")
  points(xm, pch = "+")
  points(zm, col = rgb(30L, 144L, 255L, 80L, max = 270))
  points(zm, pch = 16, cex = 2, col = rgb(30L, 144L, 255L, 80L, max = 255))
  plot(wrld_simpl, add = TRUE)
  plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)
  title(bird)
  dev.off()
  
  #save Pimage Plot
  library(maps)
  m2 <- map("world2", plot = FALSE)
  pz <- Pimage(fit, type = "intermediate")
  sz <- cut(pz, "1 months")
  cairo_pdf(file=paste("C:/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fitPimagePDF/",bird,".pdf",sep=""), height=9, width=12)
  plot(sz, col = terrain.colors(25), breaks = seq(80, 2500, length = 26), legend = FALSE, addfun = function() map(m2, add = TRUE))
  dev.off()
  
  #save Fit and Pimage (pz)
  
  save(fit,file=paste("C:/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fitRun/",bird,"_FitRun.Rdata",sep=""))
  save(pz,file=paste("C:/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fitPimage/",bird,"_FitPimage.Rdata",sep=""))
  
  
  rm(list=ls()[! ls() %in% c("filenames.fit","i","land.mask","is.sea","log.prior","nlocation","zenith","w1")])
  
  
}

  
  
  
  
  
  