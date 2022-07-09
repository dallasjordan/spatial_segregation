# Quick script to calculate active tracks and points by month and bimonthly periods
# May 22 2022
# Dallas Jordan

# Setup -------------------------------------------------------------------

library(sf)
library(dplyr)

eqc <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Load in tracking data ---------------------------------------------------

load("/Users/dallasjordan/Desktop/StonyBrook/SoMAS/Thesis/R/spatial_segregation/data/final_tracks/all_tracks_master.Rdata")
lm <- all_data[grep("lm", all_data$id), ]
lt <- all_data[grep("lt", all_data$id), ]
bm <- all_data[grep("bm", all_data$id), ]
bt <- all_data[grep("bt", all_data$id), ]

mid_data <- rbind(lm,bm)
tern_data <-rbind(lt,bt)

# fix some dates for later
mid_data$dtime <- as.Date(mid_data$dtime)
tern_data$dtime <- as.Date(tern_data$dtime)
lm$dtime <- as.Date(lm$dtime)
bm$dtime <- as.Date(bm$dtime)
lt$dtime <- as.Date(lt$dtime)
bt$dtime <- as.Date(bt$dtime)

lm$spp <- "LAAL"
lm$island <- "Midway"
bm$spp <- "BFAL"
bm$island <- "Midway"
lt$spp <- "LAAL"
lt$island <- "Tern"
bt$spp <- "BFAL"
bt$island <- "Tern"

LAAL_data <- rbind(lm,lt)
BFAL_data <- rbind(bm,bt)

LAAL_data$month <- month(LAAL_data$dtime)
LAAL_data <- LAAL_data %>% mutate(season=case_when(
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
LAAL_data <- LAAL_data %>% mutate(season_name=case_when(
  season==1 ~ "Winter",
  season==2 ~ "Spring",
  season==3 ~ "Summer",
  season==4 ~ "Fall"
))
LAAL_data <- LAAL_data %>% mutate(month = month.name[month])

BFAL_data$month <- month(BFAL_data$dtime)
BFAL_data <- BFAL_data %>% mutate(season=case_when(
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
BFAL_data <- BFAL_data %>% mutate(season_name=case_when(
  season==1 ~ "Winter",
  season==2 ~ "Spring",
  season==3 ~ "Summer",
  season==4 ~ "Fall"
))
BFAL_data <- BFAL_data %>% mutate(month = month.name[month])

# filtering out to just peak postbreeding months, other months don't have enough data:
LAAL_points_sufficient_data <- LAAL_data %>% filter(month %in% c("June","July","August","September","October","November"))
BFAL_points_sufficient_data <- BFAL_data %>% filter(month %in% c("June","July","August","September","October","November"))


# LAAL_points_sufficient_data <- LAAL_data %>% filter(month %in% c(6:11))
# BFAL_points_sufficient_data <- BFAL_data %>% filter(month %in% c(6:11))

LAAL_points_sufficient_data <- LAAL_points_sufficient_data %>% 
  mutate(month_factor = factor(month, levels=c("June","July","August","September","October","November")))
BFAL_points_sufficient_data <- BFAL_points_sufficient_data %>% 
  mutate(month_factor = factor(month, levels=c("June","July","August","September","October","November")))

LAAL_points_bimonthly <- LAAL_points_sufficient_data %>%
  mutate(bimonthly = case_when(
    grepl("June", month) ~ "Jun. - Jul.",
    grepl("July", month) ~ "Jun. - Jul.",
    grepl("August", month) ~ "Aug. - Sept.",
    grepl("September", month) ~ "Aug. - Sept.",
    grepl("October", month) ~ "Oct. - Nov.",
    grepl("November", month) ~ "Oct. - Nov."
  ))
LAAL_points_bimonthly <- LAAL_points_bimonthly %>%
  mutate(bimonthly_factor = factor(bimonthly, levels=c("Jun. - Jul.","Aug. - Sept.", "Oct. - Nov.")))

BFAL_points_bimonthly <- BFAL_points_sufficient_data %>%
  mutate(bimonthly = case_when(
    grepl("June", month) ~ "Jun. - Jul.",
    grepl("July", month) ~ "Jun. - Jul.",
    grepl("August", month) ~ "Aug. - Sept.",
    grepl("September", month) ~ "Aug. - Sept.",
    grepl("October", month) ~ "Oct. - Nov.",
    grepl("November", month) ~ "Oct. - Nov."
  ))
BFAL_points_bimonthly <- BFAL_points_bimonthly %>%
  mutate(bimonthly_factor = factor(bimonthly, levels=c("Jun. - Jul.","Aug. - Sept.", "Oct. - Nov.")))

summary_LAAL_jj <- LAAL_points_bimonthly %>% filter(bimonthly_factor=="Jun. - Jul.")
sum_midway <- summary_LAAL_jj %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_jj %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_jj %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_jj %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_jj %>% group_by(spp) %>% summarize(n=n())
sum_total_points


summary_LAAL_as <- LAAL_points_bimonthly %>% filter(bimonthly_factor=="Aug. - Sept.")
sum_midway <- summary_LAAL_as %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_as %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_as %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_as %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_as %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_LAAL_on <- LAAL_points_bimonthly %>% filter(bimonthly_factor=="Oct. - Nov.")
sum_midway <- summary_LAAL_on %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_on %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_on %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_on %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_on %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_jj <- BFAL_points_bimonthly %>% filter(bimonthly_factor=="Jun. - Jul.")
sum_midway <- summary_BFAL_jj %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_jj %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_jj %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_jj %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_jj %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_as <- BFAL_points_bimonthly %>% filter(bimonthly_factor=="Aug. - Sept.")
sum_midway <- summary_BFAL_as %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_as %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_as %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_as %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_as %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_on <- BFAL_points_bimonthly %>% filter(bimonthly_factor=="Oct. - Nov.")
sum_midway <- summary_BFAL_on %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_on %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_on %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_on %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_on %>% group_by(spp) %>% summarize(n=n())
sum_total_points


# Monthly -----------------------------------------------------------------

summary_LAAL_june <- LAAL_points_bimonthly %>% filter(month_factor=="June")
sum_midway <- summary_LAAL_june %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_june %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_june %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_june %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_june %>% group_by(spp) %>% summarize(n=n())
sum_total_points


summary_LAAL_july <- LAAL_points_bimonthly %>% filter(month_factor=="July")
sum_midway <- summary_LAAL_july %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_july %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_july %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_july %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_july %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_LAAL_august <- LAAL_points_bimonthly %>% filter(month_factor=="August")
sum_midway <- summary_LAAL_august %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_august %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_august %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_august %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_august %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_LAAL_september <- LAAL_points_bimonthly %>% filter(month_factor=="September")
sum_midway <- summary_LAAL_september %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_september %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_september %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_september %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_september %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_LAAL_october <- LAAL_points_bimonthly %>% filter(month_factor=="October")
sum_midway <- summary_LAAL_october %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_october %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_october %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_october %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_october %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_LAAL_november <- LAAL_points_bimonthly %>% filter(month_factor=="November")
sum_midway <- summary_LAAL_november %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_LAAL_november %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_LAAL_november %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_LAAL_november %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_LAAL_november %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_june <- BFAL_points_bimonthly %>% filter(month_factor=="June")
sum_midway <- summary_BFAL_june %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_june %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_june %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_june %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_june %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_july <- BFAL_points_bimonthly %>% filter(month_factor=="July")
sum_midway <- summary_BFAL_july %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_july %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_july %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_july %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_july %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_august <- BFAL_points_bimonthly %>% filter(month_factor=="August")
sum_midway <- summary_BFAL_august %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_august %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_august %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_august %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_august %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_september <- BFAL_points_bimonthly %>% filter(month_factor=="September")
sum_midway <- summary_BFAL_september %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_september %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_september %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_september %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_september %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_october <- BFAL_points_bimonthly %>% filter(month_factor=="October")
sum_midway <- summary_BFAL_october %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_october %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_october %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_october %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_october %>% group_by(spp) %>% summarize(n=n())
sum_total_points

summary_BFAL_november <- BFAL_points_bimonthly %>% filter(month_factor=="November")
sum_midway <- summary_BFAL_november %>% filter(island=="Midway") %>% group_by(id) %>% summarize(n=n())
nrow(sum_midway)
sum_midway_points <- summary_BFAL_november %>% filter(island=="Midway") %>% summarize(n=n())
sum_midway_points
sum_tern <- summary_BFAL_november %>% filter(island=="Tern") %>% group_by(id) %>% summarize(n=n())
nrow(sum_tern)
sum_tern_points <- summary_BFAL_november %>% filter(island=="Tern") %>% summarize(n=n())
sum_tern_points
sum_total_points <- summary_BFAL_november %>% group_by(spp) %>% summarize(n=n())
sum_total_points
