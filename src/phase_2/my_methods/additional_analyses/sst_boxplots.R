# SST Boxplots
# Dallas Jordan
# Last updated Nov 8 2021


# Setup -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)


# load in collated files -------------------------------------

setwd("/Users/dallasjordan/projects/spatial_segregation/files/data/oceanographic_data/SST/daily_averaged_SST_data")
file.list <- list.files(getwd())
load(file.list[1])
load(file.list[2])
load(file.list[3])
load(file.list[4])
load(file.list[5])

# drop dtime columns - just makes it harder, can recover that later if you need it
# allLAAL_SST_Midway <- allLAAL_SST_Midway[,-c(2,3)]
# allBFAL_SST_Midway <- allBFAL_SST_Midway[,-c(2,3)]
# allLAAL_SST_Tern <- allLAAL_SST_Tern[,-2]
# allBFAL_SST_Tern <- allBFAL_SST_Tern[,-2]

# # for monthly comparison: 
dailySST_midwayLAAL$month <- month(dailySST_midwayLAAL$dtime)
dailySST_midwayLAAL<- dailySST_midwayLAAL %>% mutate(month = month.name[month])
dailySST_midwayLAAL$year <- year(dailySST_midwayLAAL$dtime)
dailySST_midwayLAAL$island <- "Midway"
dailySST_midwayLAAL$spp <- "LAAL"
# 
dailySST_midwayBFAL$month <- month(dailySST_midwayBFAL$dtime)
dailySST_midwayBFAL<- dailySST_midwayBFAL %>% mutate(month = month.name[month])
dailySST_midwayBFAL$year <- year(dailySST_midwayBFAL$dtime)
dailySST_midwayBFAL$island <- "Midway"
dailySST_midwayBFAL$spp <- "BFAL"
# 
dailySST_ternLAAL$month <- month(dailySST_ternLAAL$dtime)
dailySST_ternLAAL<- dailySST_ternLAAL %>% mutate(month = month.name[month])
dailySST_ternLAAL$year <- year(dailySST_ternLAAL$dtime)
dailySST_ternLAAL$island <- "Tern"
dailySST_ternLAAL$spp <- "LAAL"
# 
dailySST_ternBFAL$month <- month(dailySST_ternBFAL$dtime)
dailySST_ternBFAL<- dailySST_ternBFAL %>% mutate(month = month.name[month])
dailySST_ternBFAL$year <- year(dailySST_ternBFAL$dtime)
dailySST_ternBFAL$island <- "Tern"
dailySST_ternBFAL$spp <- "BFAL"


all_SST <- rbind(dailySST_midwayBFAL,dailySST_midwayLAAL,dailySST_ternBFAL,dailySST_ternLAAL)

save(all_SST,file="/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/all_DailyAverageSST.Rdata")


# START HERE - load in master SST with everything -------------------------

load("/Users/dallasjordan/projects/spatial_segregation/files/data/oceanographic_data/SST/daily_averaged_sst_data/all_DailyAverageSST.Rdata")

# IMPORTANT NOTE - TEMP COLLECTION FOR 2010 LAAL AT TERN FAILED, all files did not record correctly. Exclude from yearly analysis

library(wesanderson)
pal <- wes_palette("Zissou1", 5, type="discrete")
ggplot(all_SST, aes(month,DailyMeanSST, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "SST by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))+
  facet_wrap(~island)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(limits = month.name)

# broken down by species/island
all_SST = all_SST %>% 
  unite(class, island, spp, sep = " ", remove = FALSE)

# excluding months with little data: 
all_SST <- all_SST %>% filter(month %in% c("June","July","August","September","October","November"))
all_SST <- all_SST %>% mutate(class_factor = factor(class, levels=c("Midway LAAL","Tern LAAL","Midway BFAL","Tern BFAL")))

pal2 <- c("firebrick","orangered","dodgerblue","royalblue4")
ggplot() +
  geom_boxplot(all_SST, mapping=aes(month,DailyMeanSST, fill=class_factor), outlier.shape = NA,position = position_dodge(preserve = "single"))+
  theme_classic() + 
  labs(title = "Daily mean SST by island and species") +
  scale_fill_manual(name="species",values=c(pal2[2],pal2[1],pal2[3],pal2[4]))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(limits = month.name[6:11]) +
  scale_y_discrete(limits=1:30,breaks=c(5,10,15,20,25)) +
  xlab("Month")




counts <- all_SST %>% count(Month,island,spp)


# Calculate averages ------------------------------------------------------

Midway_avg <- all_SST %>% filter(island=="Midway") %>% filter(spp=="BFAL") %>% group_by(month)
Midway_avg <- Midway_avg %>% summarize(Mean=mean(DailyMeanSST))
Midway_avg

Tern_avg <- all_SST %>% filter(island=="Tern") %>% filter(spp=="BFAL") %>% group_by(month)
Tern_avg <- Tern_avg %>% summarize(Mean=mean(DailyMeanSST))
Tern_avg
