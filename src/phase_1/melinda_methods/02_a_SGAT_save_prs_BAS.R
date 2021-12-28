
# Clean up BAS light data. Save final twilight hours as 'prs'


library(BAStag)
library(SGAT)

dp <- "/Users/Melinda/Dropbox/Chapter_03_Finalize/data/Tern_GLS_allrawfiles/BAS/Tern_BAS_1112_LIG/"
dp2<- "/Users/Melinda/Dropbox/Chapter_03_Finalize/data/Tern_GLS_allrawfiles/BAS/Tern_BAS_1112_SST/"


setwd(dp)

filenames<- list.files()


for (i in 1:length(filenames)) {
  
  f1 <- filenames[i]
  filenames[i]
  
  lf <- file.path(paste(dp,f1,sep=""))
  lf2 <- file.path(paste(dp2,f1,sep=""))
  
  rawdata <- read.lig(lf)
  rawtemp <- read.table(gsub("lig$", "tem", lf2), sep = ",", header = FALSE, colClasses = c("NULL", "character", "NULL", "numeric"), col.names = c("notok", "gmt", "excel", "temp"))
  
  rawtemp$Date <- as.POSIXct(strptime(rawtemp$gmt, "%d/%m/%y %H:%M:%S"))
  
  
  ## ------------------------------------------------------------------------
  start <- c(-166.286012 + 360, 23.869071)
  end <- start
  
  ## sensible simple maximum bounds for the track
  lon.min <- 110 
  lon.max <- 255
  lat.min <- 23
  lat.max <- 62
  
  
  ## ------------------------------------------------------------------------
  
  #   if don't want to process all light data in file, here's where you truncate it. 
  #   I processed all of it (minus the very begninning and end)
  #   And then truncated the tracks to just the PB season after the full estimation of trip , so that if I wanted to look at locations during the breeding season later, I'd have the locations already estimated.

 
  # if wanted to use start/end dates from SST timeseries code, to only process light from postbreeding trip, can use those dates here
  cutoff0 <- as.POSIXct("2012-1-11 00:00:00")
  cutoff1 <-  as.POSIXct("2012-11-1 00:00:00")
  
  d <- subset(rawdata, Date >= cutoff0 & Date <= cutoff1)
  

  ## ------------------------------------------------------------------------
  offset <- 0
  threshold <- 4
  im <- light.image(d, offset = offset, zmax = 64)
  # Click on night segments, 5 times!
  seeds <- tsimage.locator(im, n = 5)
  tsimage.points(seeds,offset=offset,pch=16,col="red")
  
  ## ------------------------------------------------------------------------
  ## twilight pairs
  prs <- find.twilights(d, threshold, include=seeds)
  im <- light.image(d, offset = offset)
  tsimage.points(prs$Twilight,offset=offset,pch=16,cex=0.5,
                 col=ifelse(prs$Rise,"dodgerblue","firebrick"))
  

  ## ------------------------------------------------------------------------
  # EDit twilights using interactive tool
  prs <- twilight.editW(d, prs, offset = offset, threshold = threshold)
  

  ##----save prs object------------------------------------------------
  
  dp3<-"/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/prs/"
  prs.name<-paste(dp3,strsplit(filenames[i],"\\.")[[1]][1],"prs.Rdata",sep="")
  save(prs,file=prs.name)
  
  
  rm(list=ls()[! ls() %in% c("filenames","dp","dp2")])
  