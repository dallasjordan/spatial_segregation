library(SGAT)
library(BAStag)
library(chron)
library(pracma)
library(raster)
library(maptools)
library(rgeos)
library(maptools)
## enough tracks are incomplete enough to where Im gonna nee to go through one at a time.
#shouldnt take more than a minute for each bird

sst <- read.csv("C:/Users/Melinda/Dropbox/Chapter_03_Finalize/data/metadata/gls_sst_truncate_metadata.csv")

## import metadata file with matching ids (mgc_id and fit_file_id)
fitid <- read.csv("C:/Users/Melinda/Dropbox/Chapter_03_Finalize/data/metadata/2015_gls_filenames.csv")

# Loop through each SGAT fit file to 
# 1. Truncate trip based on SST dates and add object to Rdatafile
# 2. save Rdata file with new MGC_id
# 3. Export truncated trip as a .csv with MGC_id

setwd("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fitRun/")
filenames<- list.files()
no.sst.file<-matrix(NA,length(filenames),1)

load("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fit_sst_trunc/locations/rdata/gls_2012_2912002_fit01trunc.Rdata")


i = i+1

load(filenames[i])

# parameters of interest (time and lat/lon)
time<-fit$model$time

#who is this?
n <- nchar(filenames[i])
n <- n-13
mgc_name<-substr(filenames[i],12,27)

# now match mgc_id with SST file for start and end dates

sstix<-which(sst$MGC_id==mgc_name)

year<-substr(mgc_name,5,8)


# if (year=="2011" | year =="2012"){
#   dates<-as.POSIXct(as.character(alldates$x),"%d/%m/%Y",tz="GMT")
# }else{
#   dates<-as.POSIXct(as.character(alldates$x),"%m/%d/%Y",tz="GMT")}



start.nb<-as.POSIXct(as.character(sst$start_date_sst_flag[sstix]),"%m/%d/%Y",tz="GMT")
end.nb<-as.POSIXct(as.character(sst$end_date_sst_flag[sstix]),"%m/%d/%Y",tz="GMT")

#   ifelse(is.na(start.nb) || is.na(end.nb))

xm <- location.mean(fit$x)
zm <- location.mean(fit$z)

startx<-which(abs(time-start.nb)==min(abs(time-start.nb)))
endx<- which(abs(time-end.nb)==min(abs(time-end.nb)))

xm.trunc<-xm[startx:endx,]
zm.trunc<-zm[startx:endx-1,]
time.trunc<-time[startx:endx]

plot(xm.trunc, pch = "+",xlab="",ylab="", ylim= c(13,63), xlim=c(138, 235),cex.axis=1.5)
lines(zm.trunc,type = "l")
points(zm.trunc, col = rgb(30L, 144L, 255L, 80L, max = 255))
data(wrld_simpl)
title(paste(substr(filenames[i],1,n)))
points(-166.29 + 360,23.87,pch=23, col="red", bg="yellow",cex=1)
points(xm.trunc[1,1],xm.trunc[1,2],pch=21, bg = "green")
points(xm.trunc[length(xm.trunc[,1]),1],xm.trunc[length(xm.trunc[,1]),2],pch=21, bg = "red")
plot(wrld_simpl, add = TRUE, col="grey" )
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE,col="grey")




nacoors<-identify(xm)
xm[nacoors, 1] <- NA
xm[nacoors, 2] <- NA
xm <- zoo::na.approx(xm)
zm <- trackMidpts(xm)

startx<-startx-10

startx<-10
startx<-startx-10
endx<-endx+5

endx<-dim(xm)[1]

startx<-1195
endx<-724



#save time,xm,zm Object
save(time.trunc,xm.trunc,zm.trunc,file=paste("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fit_sst_trunc/",mgc_name,"_fit01trunc.Rdata",sep=""))

#save time,xm,zm csv file
csv.m<-as.data.frame(as.character(time.trunc))
csv.m<-cbind(csv.m,xm.trunc,zm.trunc)
colnames(csv.m)<-c("GMT","xm.lon","xm.lat","zm.lon","zm.lat")
write.csv(csv.m,file=paste("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fit_sst_trunc/",mgc_name,"_fit01trunc.csv",sep=""),row.names=FALSE)

#save pdf

cairo_pdf(file=paste("/Users/Melinda/Dropbox/Chapter_03_Finalize/analysis/SGAT_2015/fit_02_08042015/fit_sst_trunc/pdf/",mgc_name,"_fit01trunc.pdf",sep=""), height=6.95, width=8.25,pointsize=12)
plot(xm.trunc, type="n",xlab="",ylab="", ylim= c(23,63), xlim=c(138, 235),cex.axis=1.5)
plot(wrld_simpl, add = TRUE, col="grey" )
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE,col="grey")
points(xm.trunc, pch = 16, cex = 3, col = rgb(30L, 144L, 255L, 80L, max = 255),xlab="",ylab="", ylim= c(23,63), xlim=c(138, 235),cex.axis=1.5)
lines(zm.trunc)
data(wrld_simpl)
title(paste(substr(filenames[i],1,n)))
points(-166.29 + 360,23.87,pch=23, col="red", bg="yellow",cex=3)
dev.off()


rm(list=ls()[! ls() %in% c("filenames","i","sst","fitid","no.sst.file")])

