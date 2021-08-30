# Script to do permutations for "side project" with Julia and Lesley, showing how scale impacts overlap calculations and 
# that previous methods don't necessarily work for our data set
# Dallas Jordan
# August 12

# This script 1. Loads in data,
#             2. Makes kernelUD (estUD) objects,
#             3. Permutes following Clay et al.,
#             4. Calculates Williamson's SOI,
#             5. Permutes for WSOI

# Setup -------------------------------------------------------------------

library(adehabitatHR)
library(dplyr)

library(dplyr)
library(maptools)
library(rgdal)
# library(GeoLocTools)
# setupGeolocation()
# Not available for R>=4.0
library(mapdata)
library(rgeos)
library(raster)
library(adehabitatHR)
library(sp)
library(sf)

# Load in tracks for species A and species B

sppA_read <- file.choose()
sppA <- readRDS(sppA_read)

sppB_read <- file.choose()
sppB <- readRDS(sppB_read)

all_data <- rbind(sppA, sppB)

all_data_1<-all_data[,1:3]
sp::coordinates(all_data_1) <- c("lon", "lat")

# Calculate all estUD for each animalID: href bandwidth, arbitrary grid=100 -> nicer images but longer calculation time
results <- kernelUD(all_data_1, grid=100,same4all=T,extent=0.1)
image(results[[100]])

# seperate by species
names <- names(results)
a_kde <- results[grep("A",names)]
b_kde <- results[grep("B",names)]


# Average KDEs by species  ------------------------------------------------

# average A #
a_holder <- numeric(length = 6700)
for (i in 1:length(a_kde)) {
  a_kde[[i]]@data$ud[is.na(a_kde[[i]]@data$ud)] <- 0
  add <- a_kde[[i]]@data$ud
  a_holder <- a_holder+add
}
a_holder <- a_holder/length(a_kde)
a_holder

##### modify existing estUD object with averaged values, then rename
a_averaged_estUD <- a_kde[[1]]
a_averaged_estUD@data$ud <- a_holder
a_averaged_estUD@data$ud[is.na(a_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(a_kde[[1]])
image(a_averaged_estUD) 

# average B #
b_holder <- numeric(length = 6700)
for (i in 1:length(b_kde)) {
  b_kde[[i]]@data$ud[is.na(b_kde[[i]]@data$ud)] <- 0
  add <- b_kde[[i]]@data$ud
  b_holder <- b_holder+add
}
b_holder <- b_holder/length(b_kde)
b_holder

##### modify existing estUD object with averaged values, then rename
b_averaged_estUD <- b_kde[[1]]
b_averaged_estUD@data$ud <- b_holder
b_averaged_estUD@data$ud[is.na(b_averaged_estUD@data$ud)] <- 0 # ok if it sums to >1!
image(b_kde[[1]])
image(b_averaged_estUD) 

# re-assign names
sppA_UD <- a_averaged_estUD
sppB_UD <- b_averaged_estUD

# save so you don't gotta run all this again
path <- paste0(getwd(),"/overlap_sensitivity/")
save(sppA_UD, file=paste0(path,"sppA_averaged_UD_case1.Rdata"))
save(sppB_UD, file=paste0(path,"sppB_averaged_UD_case1.Rdata"))

# load back in later
  # for PC
  setwd("E:/project_data/spatial_segregation/data/overlap_sensitivity/")
  path <- getwd()
  file_list <- list.files(path=path)
  
  load(file_list[3])
  load(file_list[4])
  # for Mac
  # FILL THIS IN WHEN YOU TRANSFER TO MAC, 8/27/2021 setwd(
  path <- getwd()
  file_list <- list.files(path=path)
  
  for (i in 1:length(file_list)){
    load(file_list[i])
  }

### done loading, now convert to estUDm class ###
  
A_v_B <- list(sppA_UD,sppB_UD)
class(A_v_B)<-"estUDm"
names(A_v_B)<-c("A","B")
# can run this to test things
# image(A_v_B)
# vert95_sppA <- getverticeshr(sppA_UD,percent=95)
# vert95_sppB <- getverticeshr(sppB_UD,percent=95)
# plot(vert95_sppA, col=NA, border="red", main="spp. A in red, spp. B in blue, 95th UD")
# plot(vert95_sppB, col=NA, border="blue",add=T)

ab_UDOI_95 <- kerneloverlaphr(A_v_B, method="UDOI", percent=95, conditional=F)
ab_UDOI_50 <- kerneloverlaphr(A_v_B, method="UDOI", percent=50, conditional=T)
ab_PHR_95 <- kerneloverlaphr(A_v_B, method="PHR", percent=95, conditional=T)
ab_PHR_50 <- kerneloverlaphr(A_v_B, method="PHR", percent=50, conditional=T)
ab_95_BA <- kerneloverlaphr(A_v_B, method="BA", percent=95, conditional = T) 
ab_50_BA <- kerneloverlaphr(A_v_B, method="BA", percent=50, conditional = T)
ab_UDOI_95_test_stat <- ab_UDOI_95[1,2]
ab_UDOI_50_test_stat <- ab_UDOI_50[1,2]
ab_PHR_95_test_stat <- ab_PHR_95[1,2]
ab_PHR_50_test_stat <- ab_PHR_50[1,2]
ab_95_BA_test_stat <- ab_95_BA[1,2]
ab_50_BA_test_stat <- ab_50_BA[1,2]  

path <- "E:/project_data/spatial_segregation/data/overlap_sensitivity/test_stats"
setwd(path)
save(list = ls(all.names = TRUE), file = "all_test_stats.RData")

# Rasters/Contours and GGplot visualization -------------------------------

### allLAAL v allBFAL ###

#95% contours and 100% UD RASTERS
# LAAL 
image(sppA_UD)
plot(getverticeshr(sppA_UD,percent=95), add=T)
vert95_sppA <- getverticeshr(sppA_UD,percent=95)
gg95_allLAAL <- fortify(vert95_sppA)
# path <- 
# save(vert95_sppA,file=path)

# for plotting UD raster in overallLAALBFAL_mapping: 
sppA.ud.vol <- getvolumeUD(sppA_UD, standardize=TRUE)
plot(sppA.ud.vol)
sppA.ud.vol.raster <- raster(sppA.ud.vol)
plot(sppA.ud.vol.raster)
# path <- 
# save(sppA.ud.vol.raster,file=path)



# Permutations of UDOI ----------------------------------------------------

# can run independent of the above
# prep and bookkeeping

sppA_read <- file.choose()
sppA <- readRDS(sppA_read)

sppB_read <- file.choose()
sppB <- readRDS(sppB_read)

all_data <- rbind(sppA, sppB)

all_data_1<-all_data[,1:3]
sp::coordinates(all_data_1) <- c("lon", "lat")

# Calculate all estUD for each animalID: href bandwidth, arbitrary grid=100 -> nicer images but longer calculation time
results <- kernelUD(all_data_1, grid=100,same4all=T,extent=0.1)
image(results[[100]])
all_ud <- results

path <- "E:/project_data/spatial_segregation/data/overlap_sensitivity/test_stats"
setwd(path)
load("all_test_stats.RData") # generated in earlier part of this script

add_a <- sppA
add_b <- sppB
add_a$spp <- "A"
add_b$spp <- "B"
resample_this <- rbind(add_a,add_b)
resample_this_spp <- split(resample_this, resample_this$spp)
resample_this_a_track <- split(resample_this_spp[[1]],resample_this_spp[[1]]$Animal_ID)
resample_this_b_track <- split(resample_this_spp[[2]],resample_this_spp[[2]]$Animal_ID)
resample_all_tracks <- c(resample_this_a_track,resample_this_b_track)
available_a_tracks <- unique(resample_this_spp[[1]]$Animal_ID)
available_b_tracks <- unique(resample_this_spp[[2]]$Animal_ID)
n_aat <- length(available_a_tracks) # you are calculating these to preserve sample size in randomizations
n_abt <- length(available_b_tracks)
n_all <- n_aat+n_abt
counts <- resample_this %>% count(spp,Animal_ID)

# 95th/50th isopleth and BA
iter <- 1000
significance_tallyUDOI_95 <- 0
significance_tallyUDOI_50 <- 0
significance_tallyBA_95 <- 0
significance_tallyBA_50 <- 0
iter_tally <- 0
resultsUDOI_95_storage <- vector(length=iter, mode="numeric")
resultsUDOI_50_storage <- vector(length=iter, mode="numeric")
resultsBA_95_storage <- vector(length=iter, mode="numeric")
resultsBA_50_storage <- vector(length=iter, mode="numeric")

for (i in 1:iter){
  pool <- 1:n_all
  a_nums <- sample(pool, n_aat, replace=F)
  b_nums <- setdiff(pool,a_nums)
  a_iter <- list()
  b_iter <- list()
  
  for (j in 1:length(a_nums)){
    a_loop <- all_ud[[a_nums[j]]]
    a_iter <- append(a_iter,a_loop)
  }
  
  for (w in 1:length(b_nums)){
    b_loop <- all_ud[[b_nums[w]]]
    b_iter <- append(b_iter,b_loop)
  }
  
  # average a_iter #
  a_iter_holder <- numeric(length = 6700)
  for (d in 1:length(a_iter)) {
    a_iter[[d]]@data$ud[is.na(a_iter[[d]]@data$ud)] <- 0
    add <- a_iter[[d]]@data$ud
    a_iter_holder <- a_iter_holder+add
  }
  a_iter_holder <- a_iter_holder/length(a_iter)
  a_iter_holder
  
  ##### modify existing estUD object with averaged values, then rename
  a_iter_avg <- a_iter[[1]]
  a_iter_avg@data$ud <- a_iter_holder
  a_iter_avg@data$ud[is.na(a_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
  #image(a_iter[[1]])
  #image(a_iter_avg) 
  
  # average b_iter #
  b_iter_holder <- numeric(length = 6700)
  for (f in 1:length(b_iter)) {
    b_iter[[f]]@data$ud[is.na(b_iter[[f]]@data$ud)] <- 0
    add <- b_iter[[f]]@data$ud
    b_iter_holder <- b_iter_holder+add
  }
  b_iter_holder <- b_iter_holder/length(b_iter)
  b_iter_holder
  
  ##### modify existing estUD object with averaged values, then rename
  b_iter_avg <- b_iter[[1]]
  b_iter_avg@data$ud <- b_iter_holder
  b_iter_avg@data$ud[is.na(b_iter_avg@data$ud)] <- 0 # ok if it sums to >1!
  #image(b_iter[[1]])
  #image(b_iter_avg) 
  
  # Create estUDm 
  i_sppA <-a_iter_avg
  i_sppB <-b_iter_avg
  A_v_B <- list(i_sppA,i_sppB)
  class(A_v_B)<-"estUDm"
  names(A_v_B)<-c("i_sppA","i_sppB")
  #image(A_v_B)
  
  # can run this to test things
  # image(A_v_B)
  # plot(getverticeshr(i_sppA,percent=95), add=T)
  # vert95_i_sppA <- getverticeshr(i_sppA,percent=95)
  # vert95_i_sppB <- getverticeshr(i_sppB,percent=95)
  # plot(vert95_i_sppA, col=NA, border="red", main="spp. A in red, spp. B in blue, 95th UD")
  # plot(vert95_i_sppB, col=NA, border="blue",add=T)
  
  iter_ab_UDOI_95 <- kerneloverlaphr(A_v_B, method="UDOI", percent=95, conditional=T)
  iter_ab_UDOI_50 <- kerneloverlaphr(A_v_B, method="UDOI", percent=50, conditional=T)
  iter_ab_BA_95 <- kerneloverlaphr(A_v_B, method="BA", percent=95, conditional=T) 
  iter_ab_BA_50 <- kerneloverlaphr(A_v_B, method="BA", percent=50, conditional=T)
  
  iter_tally <- iter_tally+1
  print(iter_tally)
  
  # proportion of randomized overlaps that are less than observed
  # so if all randomized overlaps are greater than the observed, p=0, there is significant segregation
  # if the observed overlap was less than all randomizations, then P ≤ 0.001, reject the null hypothesis that 
  # there is no difference in the spatial distributions of the two groups
  
  resultsUDOI_95_storage[i] <- iter_ab_UDOI_95[1,2]
  print(resultsUDOI_95_storage[i])
  if (iter_ab_UDOI_95[1,2]<ab_UDOI_95_test_stat){
    print(iter_ab_UDOI_95)
    print("+1!")
    significance_tallyUDOI_95 <- significance_tallyUDOI_95+1
  }
  
  resultsUDOI_50_storage[i] <- iter_ab_UDOI_50[1,2]
  print(resultsUDOI_50_storage[i])
  if (iter_ab_UDOI_50[1,2]<ab_UDOI_50_test_stat){
    print(iter_ab_UDOI_50)
    print("+1!")
    significance_tallyUDOI_50 <- significance_tallyUDOI_50+1
  }
  
  resultsBA_95_storage[i] <- iter_ab_BA_95[1,2]
  print(resultsBA_95_storage[i])
  if (iter_ab_BA_95[1,2]<ab_95_BA_test_stat){
    print(iter_ab_BA_95)
    print("+1!")
    significance_tallyBA_95 <- significance_tallyBA_95+1
  }
  
  resultsBA_50_storage[i] <- iter_ab_BA_50[1,2]
  print(resultsBA_50_storage[i])
  if (iter_ab_BA_50[1,2]<ab_50_BA_test_stat){
    print(iter_ab_BA_50)
    print("+1!")
    significance_tallyBA_50 <- significance_tallyBA_50+1
  }
}

p_valueUDOI_95 <- significance_tallyUDOI_95/iter
mean_valueUDOI_95 <- mean(resultsUDOI_95_storage)
sd_valueUDOI_95 <- sd(resultsUDOI_95_storage)

p_valueUDOI_95 # 0
mean_valueUDOI_95 # 1.286
sd_valueUDOI_95 # 0.055

p_valueUDOI_50 <- significance_tallyUDOI_50/iter
mean_valueUDOI_50 <- mean(resultsUDOI_50_storage)
sd_valueUDOI_50 <- sd(resultsUDOI_50_storage)

p_valueUDOI_50 # 0
mean_valueUDOI_50 # 0.149
sd_valueUDOI_50 # 0.025

p_valueBA_95 <- significance_tallyBA_95/iter
mean_valueBA_95 <- mean(resultsBA_95_storage)
sd_valueBA_95 <- sd(resultsBA_95_storage)

p_valueBA_95 # 0
mean_valueBA_95 # 0.887
sd_valueBA_95 # 0.013

p_valueBA_50 <- significance_tallyBA_50/iter
mean_valueBA_50 <- mean(resultsBA_50_storage)
sd_valueBA_50 <- sd(resultsBA_50_storage)

p_valueBA_50 # 0
mean_valueBA_50 # 0.374
sd_valueBA_50 # 0.033




