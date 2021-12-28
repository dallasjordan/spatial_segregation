# Demonstration of "random" under CRW vs "random" under track switching
# Steps: 
#   1. Go to totalOL script for calculating track permutation and grab random points
#   2. Load those in here. plot with diff colors
#   3. Load in the CRW random points Julia sent you
#   4. filter for the points you need plot with diff colors

library(dpylr)
library(data.table)


# Track permutation points and plot ---------------------------------------

all_data_read <- file.choose() #read in .rds that has all the originally simulated points (e.g. xxxx_fullOL-20211103.rds)
all_data <- readRDS(all_data_read)
sppA <- all_data %>% dplyr::filter(species=="A")
sppB <- all_data %>% dplyr::filter(species=="B")

df_list <- split(all_data, f=all_data$Animal_ID)

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

pool <- 1:n_all
a_nums <- sample(pool, n_aat, replace=F)
b_nums <- setdiff(pool,a_nums)
a_iter <- list()
b_iter <- list()

for (j in 1:length(a_nums)){
  a_iter[[j]] <- df_list[[a_nums[j]]]
}

for (w in 1:length(b_nums)){
  b_iter[[w]] <- df_list[[b_nums[w]]]
}


# now you have two lists of 100 dataframes each of randomized, id-permuted tracks. Combine into big dataframes and plot
library(data.table)
a.big_df <- rbindlist(a_iter)
b.big_df <- rbindlist(b_iter)

col_list = c("#983BADA3","#7CEBE7FA")
track_permutation_pts <- ggplot()+
  geom_point(data=b.big_df,aes(x=lon,y=lat,col=species),size=0.005)+
  geom_point(data=a.big_df,aes(x=lon,y=lat,col=species),size=0.005)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  labs(title="Randomized points through permuting track ID", color="Species")+
  theme(plot.title=element_text(size=11))+
  coord_cartesian(xlim=c(-200,200),ylim=c(-150,150))
track_permutation_pts

# CRW points and plot -----------------------------------------------------

CRW_allpoints<-readRDS(file.choose()) # totalOL_allpoints.rds
# just need one set of simulations
CRW_allpoints <- CRW_allpoints[[1]]
# this set has 100 tracks from each species, 
CRW_allpoints_df <- as.data.frame(CRW_allpoints)
CRW_allpoints_a <- CRW_allpoints_df %>% filter(spp=="A")
CRW_allpoints_b <- CRW_allpoints_df %>% filter(spp=="B")
  
col_list = c("#983BADA3","#7CEBE7FA")
CRW_pts <- ggplot()+
  geom_point(data=CRW_allpoints_b,aes(x=lon,y=lat,col=spp),size=0.005)+
  geom_point(data=CRW_allpoints_a,aes(x=lon,y=lat,col=spp),size=0.005)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  labs(title="CRW-generated random tracks", color="Species")+
  theme(plot.title=element_text(size=12))+
  coord_cartesian(xlim=c(-200,200),ylim=c(-150,150))

CRW_pts
  

# Each plot has 100 tracks from each species, each track has 350 points, for both panels
# There were 1000 iterations of these simulations, you just did 1 simulation with 100 tracks to 
# match your other method
# will need to explain which one was significant in methods


# Plot together -----------------------------------------------------------



library(ggpubr)

combined_figure <- ggarrange(track_permutation_pts, CRW_pts+rremove("ylab"), 
                             common.legend=T, legend="bottom",
                             ncol=2, nrow=1)
combined_figure

