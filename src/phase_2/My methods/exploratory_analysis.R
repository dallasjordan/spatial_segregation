# Exploratory visualization of segregation
# Dallas Jordan
# Dec 18 2020

library(tidyverse)
library(readr)
library(readxl)
library(gridExtra)
library(wesanderson)
library(lattice)
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggpubr)

#########################################################################################

convlon <- function(lon){
  ifelse(lon < 0, lon + 360, lon)
}
years <- c("2008","2009","2010","2011","2012")

# combine postbreeding data into 4 datasets, by island/species

  ## Midway 
  
    SPECIES = "LAAL"
    
    setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
                 SPECIES,"/"))
    data1<-data.frame()
    for (i in 1:5){
      files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
      setwd(paste0(getwd(),"/",years[i]))
      data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
      colnames(data_loop)= c("dtime","x","y")
      data_loop$id <- "midwayLAAL"
      data_loop$year <- years[i]
      data1 <- rbind(data1,data_loop)
      setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
                   SPECIES,"/"))
    }
    
    data1 <- data1[!is.na(data1$x) & !is.na(data1$y),]
    data1.sp <- data1[,c("id","x","y")]
    
    SPECIES = "BFAL"
    
    setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
                 SPECIES,"/"))
    data2<-data.frame()
    for (i in 1:5){
      files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
      setwd(paste0(getwd(),"/",years[i]))
      data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
      colnames(data_loop)= c("dtime","x","y")
      data_loop$id <- paste0("midwayBFAL_",years[i])
      data2 <- rbind(data2,data_loop)
      setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/midway_postbreeding_exports/",
                   SPECIES,"/"))
    }
    
    data2 <- data2[!is.na(data2$x) & !is.na(data2$y),]
    data2.sp <- data2[,c("id","x","y")]
    
  ## Tern
    
    SPECIES = "LAAL"
    
    setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
                 SPECIES,"/"))
    data3<-data.frame()
    for (i in 1:5){
      files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
      setwd(paste0(getwd(),"/",years[i]))
      data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
      colnames(data_loop)= c("dtime","x","y")
      data_loop$id <- paste0("ternLAAL_",years[i])
      data3 <- rbind(data3,data_loop)
      setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
                   SPECIES,"/"))
    }
    
    data3 <- data3[!is.na(data3$x) & !is.na(data3$y),]
    data3.sp <- data3[,c("id","x","y")]
    
    SPECIES = "BFAL"
    
    setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
                 SPECIES,"/"))
    data4<-data.frame()
    for (i in 1:5){
      files <- list.files(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",SPECIES,"/",years[i]), pattern = "csv")
      setwd(paste0(getwd(),"/",years[i]))
      data_loop <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
      colnames(data_loop)= c("dtime","x","y")
      data_loop$id <- paste0("ternBFAL_",years[i])
      data4 <- rbind(data4,data_loop)
      setwd(paste0("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/tern_postbreeding_exports/",
                   SPECIES,"/"))
    }
    
    data4 <- data4[!is.na(data4$x) & !is.na(data4$y),]
    data4.sp <- data4[,c("id","x","y")]

# combine into one giant dataset

  all.data <- rbind(data1.sp,data2.sp,data3.sp,data4.sp)
  
#############################################

# Longitude boxplot by spp. and by year
  
  mL2008 <- data1[data1$year==2008,]
  mL2009 <- data1[data1$year==2009,]
  mL2010 <- data1[data1$year==2010,]
  mL2011 <- data1[data1$year==2011,]
  mL2012 <- data1[data1$year==2012,]
  
  mB2008 <- data2[data2$id=="midwayBFAL_2008",]
  mB2009 <- data2[data2$id=="midwayBFAL_2009",]
  mB2010 <- data2[data2$id=="midwayBFAL_2010",]
  mB2011 <- data2[data2$id=="midwayBFAL_2011",]
  mB2012 <- data2[data2$id=="midwayBFAL_2012",]
  
  tL2008 <- data3[data3$id=="ternLAAL_2008",]
  tL2009 <- data3[data3$id=="ternLAAL_2009",]
  tL2010 <- data3[data3$id=="ternLAAL_2010",]
  tL2011 <- data3[data3$id=="ternLAAL_2011",]
  tL2012 <- data3[data3$id=="ternLAAL_2012",]
  
  tB2008 <- data4[data4$id=="ternBFAL_2008",]
  tB2009 <- data4[data4$id=="ternBFAL_2009",]
  tB2010 <- data4[data4$id=="ternBFAL_2010",]
  tB2011 <- data4[data4$id=="ternBFAL_2011",]
  tB2012 <- data4[data4$id=="ternBFAL_2012",]
  
  pal <- wes_palette("Zissou1", 5, type="discrete") 
  lon_box_mL2008 <- ggplot(mL2008, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2008") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  lon_box_mL2009 <- ggplot(mL2009, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2009") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  lon_box_mL2010 <- ggplot(mL2010, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2010") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  lon_box_mL2011 <- ggplot(mL2011, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2011") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  lon_box_mL2012 <- ggplot(mL2012, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2012") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  longitude_figure_mL <- ggarrange(lon_box_mL2008, lon_box_mL2009, lon_box_mL2010,
                                    lon_box_mL2011,lon_box_mL2012,
                                ncol = 3, nrow = 2)
  longitude_figure_mL
  
  lon_box_mB2008 <- ggplot(mB2008, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2008") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  lon_box_mB2009 <- ggplot(mB2009, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2009") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  lon_box_mB2010 <- ggplot(mB2010, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2010") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  lon_box_mB2011 <- ggplot(mB2011, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2011") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  lon_box_mB2012 <- ggplot(mB2012, aes(id, x, fill=factor(id))) +
    geom_boxplot() +
    theme_classic() +
    labs(title ="2012") + 
    theme(legend.position = "none") +
    scale_fill_manual(name="presence", values=c(pal[1],pal[2], pal[3], pal[4]))
  
  
  longitude_figure_mB <- ggarrange(lon_box_mB2008, lon_box_mB2009, lon_box_mB2010,
                                lon_box_mB2011,lon_box_mB2012,
                                ncol = 3, nrow = 2)
  longitude_figure_mB

library(dplyr)

## Summary stats table for all your vars
my_vars <- your_data %>% select(., covar1:covarN) # select all the metric variables you want to include

my_summary <- my_vars %>% summarise_each(funs(
  n=sum(!is.na(.)),
  Mean = mean(.,na.rm = TRUE),
  Median = median(.,na.rm = TRUE),
  SD = sd(.,na.rm = TRUE),
  Min = min(.,na.rm = TRUE),
  Q25 = quantile(., 0.25, na.rm = T),
  Q75 = quantile(., 0.75, na.rm = T),
  Max = max(.,na.rm = TRUE)
))

my_table <- my_summary %>% tidyr::gather(stat, val) %>%
  tidyr::separate(stat, into = c("var", "stat"), sep = "_") %>%
  tidyr::spread(stat, val)


## Count observations
ggplot(your_data) +
  geom_bar(aes(factor(manatee)))


## Count observations separating by another covariate
ggplot(your_data) +
  geom_bar(aes(factor(manatee))) +
  facet_wrap(~covariate)


## Violin/Boxplots (distribution of covariates by manatee/no manatee)
ggplot(EMA, aes(factor(manatee), covariate)) +
  geom_violin() +
  geom_boxplot() +
  geom_point(alpha=.4)


## Multiple covariates (with similar scale) in the same plot
long_data <- your_data %>% gather(var,value, covar1:covarN) # gather covariates into long format

ggplot(long_data, aes(var, value, fill=factor(manatee))) +
  geom_violin() +
  geom_boxplot() +
  geom_point(alpha=.4)


#######################################################################

manateedata <- read_excel("assignment3/manateedata.xls")
colnames(manateedata) <- c("sample","season","month","manatee","depth","transp","temp","sal","gravel","crsSand","fineSand","mud")
manatee_long <- manateedata %>% 
  pivot_longer(
    cols = depth:mud,
    names_to = "env_variable")
pal <- wes_palette("Zissou1", 5, type="discrete") # fun colors

# Prepare summary stats
my_vars <- manateedata %>% select(., depth:mud)

my_summary <- my_vars %>% summarise_each(funs(
  n=sum(!is.na(.)),
  Mean = mean(.,na.rm = TRUE),
  Median = median(.,na.rm = TRUE),
  SD = sd(.,na.rm = TRUE),
  Min = min(.,na.rm = TRUE),
  Q25 = quantile(., 0.25, na.rm = T),
  Q75 = quantile(., 0.75, na.rm = T),
  Max = max(.,na.rm = TRUE)
))

# Summary stats in table

summary_stats_table <- my_summary %>% tidyr::gather(stat, val) %>%
  tidyr::separate(stat, into = c("var", "stat"), sep = "_") %>%
  tidyr::spread(stat, val)

kable(summary_stats_table,
      digits=3, booktabs=TRUE,
      caption='Summary statistics of collected data') %>%
  kable_styling(latex_options = "hold_position") 

# Total counts of presences absence
total_counts <- ggplot(manateedata) +
  geom_bar(aes(factor(manatee), fill=factor(manatee))) +
  labs(x="absence/presence") +
  theme_classic() +
  ggtitle("Total counts of absence/presence of manatees") +
  theme(plot.title = element_text(size=8, face="bold.italic"), legend.position="none") +
  scale_fill_manual(name="presence",values=c(pal[1],pal[3])) 

## Count observations separating by another covariate
counts_by_season <- ggplot(manateedata) +
  geom_bar(aes(factor(manatee), fill=factor(manatee))) +
  labs(x="absence/presence") +
  facet_wrap(~season) +
  theme_classic() +
  ggtitle("Counts of absence/presence of manatees by season") +
  theme(plot.title = element_text(size=8, face="bold.italic"),legend.position="none") +
  scale_fill_manual(name="presence",values=c(pal[1],pal[3])) 

## Count observations separating by month
counts_by_month<- ggplot(manateedata) +
  geom_bar(aes(factor(manatee), fill=factor(manatee))) +
  labs(x="absence/presence") +
  facet_wrap(~month) +
  theme_classic() +
  ggtitle("Counts of absence/presence of manatees by month") +
  theme(plot.title = element_text(size=8, face="bold.italic")) +
  scale_fill_manual(name="presence",values=c(pal[1],pal[3]))

lay <- rbind(c(1,2),
             c(3,2))
grid.arrange(total_counts,counts_by_month,counts_by_season, ncol=2, layout_matrix=lay)

## Violin/Boxplots (distribution of covariates by manatee/no manatee)
## All covariates in the same plot
manatee_long <- manateedata %>%  # gather covariates into long format
  pivot_longer(
    cols = depth:mud,
    names_to = "env_variable")

all_covar <- ggplot(manatee_long, aes(env_variable, value, fill=factor(manatee))) +
  geom_boxplot() +
  theme_classic() +
  labs(x="environmental variable category", title ="Measured environmental covariates by presence") + 
  theme(legend.position = "right") +
  scale_fill_manual(name="presence", values=c(pal[1],pal[3]))

all_covar


