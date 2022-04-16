# Wind Boxplots
# Dallas Jordan
# Last updated Nov 8   2021 - some of the pathways will be messed up, you added one folder layer for better organization

# 2M and 10M are in this script

# Setup -------------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(readr)

years <- c("2007","2008","2009","2010","2011","2012")


# Load in data ------------------------------------------------------------

# merging files and loading everything in by species and island - these files have id, dtime, x, y, track, year, month, day, daily_avg u per height, daily_avg v per height, wind speed, wind direction 
# Tern 
Species <- "BFAL"
Island <- "tern"
years <- c(2008,2009,2010,2011,2012)
setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/2M/tern',Species,"/"))

file.list <- list.files(getwd())

big.df <- data.frame()
for (i in 1:length(file.list)){
  add.file <- read_csv(file.list[i])
  add.file$island <- "Tern"
  add.file$spp <- Species
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

#write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/2M/allBFAL_wind2M_Tern.csv"),row.names = FALSE)

# Midway
Species <- "BFAL"
years <- c(2008,2009,2010,2011,2012)
setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/2M/midway',Species,"/"))

file.list <- list.files(getwd())

big.df <- data.frame()
for (i in 1:length(file.list)){
  add.file <- read_csv(file.list[i])
  add.file$island <- "Midway"
  add.file$spp <- Species
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

#write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/2M/allBFAL_wind2M_Midway.csv"),row.names = FALSE)

# load in collated files -------------------------------------

allLAAL_wind2M_Midway <- read_csv('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/2M/allLAAL_wind2M_Midway.csv')
allBFAL_wind2M_Midway <- read_csv('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/2M/allBFAL_wind2M_Midway.csv')
allLAAL_wind2M_Tern <- read_csv('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/2M/allLAAL_wind2M_Tern.csv')
allBFAL_wind2M_Tern <- read_csv('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/2M/allBFAL_wind2M_Tern.csv')

all_wind2M<- rbind(allLAAL_wind2M_Midway,allBFAL_wind2M_Midway,allLAAL_wind2M_Tern ,allBFAL_wind2M_Tern)

# can drop Record number
all_wind2M<- all_wind2M[,-1]
write.csv(all_wind2M,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/all_wind2M.csv"),row.names = FALSE)


# START HERE - load in master wind with everything -------------------------

all_wind2M <- read_csv("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/all_wind2M.csv")
all_wind2M$mths<-as.numeric(all_wind2M$mths)
all_wind2M <- all_wind2M %>% mutate(MonthName = month.name[mths])
counts_wind2M <- all_wind2M %>% count(MonthName, spp, island)

# very basic 
library(wesanderson)
pal <- wes_palette("Zissou1", 5, type="discrete")
ggplot(all_wind2M, aes(island,wind_speed2M, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "daily averaged wind2M by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))

# by month - this is the one you probably want
ggplot(all_wind2M, aes(MonthName,wind_speed2M, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "daily averaged wind2M by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))+
  facet_wrap(.~island)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(limits = month.name) +
  xlab("Month")

# wind direction by month
ggplot(all_wind2M, aes(MonthName,wind_direction2M, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "daily averaged wind2M direction by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))+
  facet_wrap(.~island)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(limits = month.name) +
  xlab("Month")




#==============================================
# REPEAT AGAIN FOR 10M
#==============================================

# Load in data ------------------------------------------------------------

# merging files and loading everything in by species and island - these files have id, dtime, x, y, track, year, month, day, daily_avg u per height, daily_avg v per height, wind speed, wind direction 
# Tern 
Species <- "BFAL"
Island <- "tern"
years <- c(2008,2009,2010,2011,2012)
setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/10M/tern',Species,"/"))

file.list <- list.files(getwd())

big.df <- data.frame()
for (i in 1:length(file.list)){
  add.file <- read_csv(file.list[i])
  add.file$island <- "Tern"
  add.file$spp <- Species
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

#write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/10M/allBFAL_wind10M_Tern.csv"),row.names = FALSE)

# Midway
Species <- "BFAL"
years <- c(2008,2009,2010,2011,2012)
setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/10M/midway',Species,"/"))

file.list <- list.files(getwd())

big.df <- data.frame()
for (i in 1:length(file.list)){
  add.file <- read_csv(file.list[i])
  add.file$island <- "Midway"
  add.file$spp <- Species
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/10M/allBFAL_wind10M_Midway.csv"),row.names = FALSE)

# load in collated files -------------------------------------

allLAAL_wind10M_Midway <- read_csv('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/10M/allLAAL_wind10M_Midway.csv')
allBFAL_wind10M_Midway <- read_csv('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/10M/allBFAL_wind10M_Midway.csv')
allLAAL_wind10M_Tern <- read_csv('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/10M/allLAAL_wind10M_Tern.csv')
allBFAL_wind10M_Tern <- read_csv('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/wind/wind_extracted/10M/allBFAL_wind10M_Tern.csv')

# fix the two incorrectly named columns in allLAAL_midway
allLAAL_wind10M_Midway <- rename(allLAAL_wind10M_Midway,daily_avg_u10m=daily_avg_u2m)
allLAAL_wind10M_Midway <- rename(allLAAL_wind10M_Midway,daily_avg_v10m=daily_avg_v2m)

all_wind10M<- rbind(allLAAL_wind10M_Midway,allBFAL_wind10M_Midway,allLAAL_wind10M_Tern ,allBFAL_wind10M_Tern)

# can drop Record number
all_wind10M<- all_wind10M[,-1]
#write.csv(all_wind10M,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/all_wind10M.csv"),row.names = FALSE)


# START HERE - load in master wind with everything -------------------------

all_wind10M <- read_csv("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/all_wind10M.csv")
all_wind10M$mths<-as.numeric(all_wind10M$mths)
all_wind10M <- all_wind10M %>% mutate(MonthName = month.name[mths])
counts_wind10M <- all_wind10M %>% count(MonthName, spp, island)

# very basic 
library(wesanderson)
pal <- wes_palette("Zissou1", 5, type="discrete")
ggplot(all_wind10M, aes(island,wind_speed10M, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "daily averaged wind10M speed by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))

# by month - this is the one you probably want
ggplot(all_wind10M, aes(MonthName,wind_speed10M, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "daily averaged wind10M speed by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))+
  facet_wrap(.~island)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(limits = month.name) +
  xlab("Month")

# broken down by species/island
# don't need these!
# all_wind10M_lm <- all_wind10M %>% filter(island=="Midway") %>% 
#   filter(spp=="LAAL")
# all_wind10M_bm <- all_wind10M %>% filter(island=="Midway") %>% 
#   filter(spp=="BFAL")
# all_wind10M_lt <- all_wind10M %>% filter(island=="Tern") %>% 
#   filter(spp=="LAAL")
# all_wind10M_bt <- all_wind10M %>% filter(island=="Tern") %>% 
#   filter(spp=="BFAL")

all_wind10M = all_wind10M %>% 
  unite(class, island, spp, sep = " ", remove = FALSE)

# excluding months with little data: 
all_wind10M <- all_wind10M %>% filter(MonthName %in% c("June","July","August","September","October","November"))

all_wind10M <- all_wind10M %>% mutate(class_factor = factor(class, levels=c("Midway LAAL","Tern LAAL","Midway BFAL","Tern BFAL")))

pal2 <- c("firebrick","orangered","dodgerblue","royalblue4")

ggplot() +
  geom_boxplot(all_wind10M, mapping=aes(MonthName,wind_speed10M, fill=class_factor), outlier.shape = NA,position = position_dodge(preserve = "single"))+
  theme_classic() + 
  labs(title = "daily averaged wind10M speed by island and species") +
  scale_fill_manual(name="species",values=c(pal2[2],pal2[1],pal2[3],pal2[4]))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(limits = month.name[6:11]) +
  xlab("Month")












# wind direction by month
ggplot(all_wind10M, aes(MonthName,wind_direction10M, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "daily averaged wind10M direction by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))+
  facet_wrap(.~island)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(limits = month.name) +
  xlab("Month")

