# chla Boxplots
# Dallas Jordan
# Last updated Nov 1 2021 - some of the pathways will be messed up, you added one folder layer for better organization


# Setup -------------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(readr)

years <- c("2007","2008","2009","2010","2011","2012")


# Load in data ------------------------------------------------------------

# merging files and loading everything in by species and island - these files have ____
# Tern 
Species <- "BFAL"
Island <- "Tern"
years <- c(2008,2009,2010,2011,2012)
setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla_extracted/tern/',Species,"/"))

file.list <- list.files(getwd())

big.df <- data.frame()
for (i in 1:length(file.list)){
  add.file <- read_csv(file.list[i])
  add.file$island <- "Tern"
  add.file$spp <- Species
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

# write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/allBFAL_chla_Tern.csv"),row.names = FALSE)

# Midway
Species <- "BFAL"
years <- c(2008,2009,2010,2011,2012)
setwd(paste0('~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/chla_extracted/midway/',Species,"/"))

file.list <- list.files(getwd())

big.df <- data.frame()
for (i in 1:length(file.list)){
  add.file <- read_csv(file.list[i])
  add.file$island <- "Midway"
  add.file$spp <- Species
  if (i==1){big.df <- add.file
  } else {big.df <- rbind(big.df, add.file)}
}

write.csv(big.df,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/allBFAL_chla_Midway.csv"),row.names = FALSE)


# load in collated files -------------------------------------

allLAAL_chla_Midway <- read_csv("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/allLAAL_chla_Midway.csv")
allBFAL_chla_Midway <- read_csv("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/allBFAL_chla_Midway.csv")
allLAAL_chla_Tern <- read_csv("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/allLAAL_chla_Tern.csv")
allBFAL_chla_Tern <- read_csv("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/allBFAL_chla_Tern.csv")

all_chla <- rbind(allLAAL_chla_Midway,allLAAL_chla_Tern,allBFAL_chla_Midway,allBFAL_chla_Tern)

# can drop Record number
all_chla <- all_chla[,-1]
write.csv(all_chla,paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/all_chla.csv"),row.names = FALSE)


# START HERE - load in master chla with everything -------------------------

all_chla <- read_csv("~/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/oceanographic_data/all_chla.csv")
all_chla$mths<-as.numeric(all_chla$mths)
all_chla <- all_chla %>% mutate(MonthName = month.name[mths])
counts_chla <- all_chla %>% count(MonthName, spp, island)
# converting numbers in month to month names so scale_x_discrete(limits=month.name works)
all_chla$mths<-as.numeric(all_chla$mths)
all_chla <- all_chla %>% mutate(MonthName = month.name[mths])

library(wesanderson)
pal <- wes_palette("Zissou1", 5, type="discrete")
ggplot(all_chla, aes(island,chla, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "chla by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))

ggplot(all_chla, aes(MonthName,chla, fill=factor(spp))) +
  geom_boxplot()+
  theme_classic() + 
  labs(title = "chlA by island and species") +
  scale_fill_manual(name="species",values=c(pal[1],pal[3]))+
  facet_wrap(.~island)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(limits = month.name)

# broken down by species/island
all_chla = all_chla %>% 
  unite(class, island, spp, sep = " ", remove = FALSE)

# excluding months with little data: 
all_chla <- all_chla %>% filter(MonthName %in% c("June","July","August","September","October","November"))
all_chla <- all_chla %>% mutate(class_factor = factor(class, levels=c("Midway LAAL","Tern LAAL","Midway BFAL","Tern BFAL")))

pal2 <- c("firebrick","orangered","dodgerblue","royalblue4")
ggplot() +
  geom_boxplot(all_chla, mapping=aes(MonthName,chla, fill=class_factor),outlier.shape = NA,position = position_dodge(preserve = "single"))+
  theme_classic() + 
  labs(title = "chlA by island and species") +
  scale_fill_manual(name="species",values=c(pal2[2],pal2[1],pal2[3],pal2[4]))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(limits = month.name[6:11]) +
  scale_y_continuous(limits=c(0,1.5))+
  ylab("ChlA (mg m^-3)")+
  xlab("Month")
