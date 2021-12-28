# Making first overlap analysis figure - showing the "blobs" from simulated data to make clear
# what our "overlap scenarios" are. 
# Should be easy - load in blobs, plot, ggarrange
# Nov 28
# Dallas Jordan

library(ggplot2)
library(dplyr)
library(sjPlot)
library(ggpubr)
library(broom)

# No OL -------------------------------------------------------------------

poly_path <- file.choose()
noOL_polys <- readRDS(poly_path)

SppA_95 <- subset(noOL_polys, name=="SppA_95")
SppA_50 <- subset(noOL_polys, name=="SppA_50")
SppA_95_df <- tidy(SppA_95)
SppA_50_df <- tidy(SppA_50)

SppB_95 <- subset(noOL_polys, name=="SppB_95")
SppB_50 <- subset(noOL_polys, name=="SppB_50")
SppB_95_df <- tidy(SppB_95)
SppB_50_df <- tidy(SppB_50)

col_list = c("#440154FF","#21908CFF")

noOL_95 <- ggplot()+
  geom_polygon(data=SppA_95_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  geom_polygon(data=SppB_95_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  scale_color_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  scale_fill_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  labs(color = "Species", fill = "Species", title="")+
  coord_cartesian(xlim=c(-200,200),ylim=c(-100,100))+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))
noOL_95

noOL_50 <- ggplot()+
  geom_polygon(data=SppA_50_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  geom_polygon(data=SppB_50_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.4)+
  scale_color_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  scale_fill_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  labs(color = "Species", fill = "Species",title="")+
  coord_cartesian(xlim=c(-200,200),ylim=c(-100,100))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
noOL_50


# Some OL -----------------------------------------------------------------

poly_path <- file.choose()
someOL_polys <- readRDS(poly_path)

SppA_95 <- subset(someOL_polys, name=="SppA_95")
SppA_50 <- subset(someOL_polys, name=="SppA_50")
SppA_95_df <- tidy(SppA_95)
SppA_50_df <- tidy(SppA_50)

SppB_95 <- subset(someOL_polys, name=="SppB_95")
SppB_50 <- subset(someOL_polys, name=="SppB_50")
SppB_95_df <- tidy(SppB_95)
SppB_50_df <- tidy(SppB_50)

col_list = c("#440154FF","#21908CFF")

someOL_95 <- ggplot()+
  geom_polygon(data=SppA_95_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  geom_polygon(data=SppB_95_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  scale_color_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  scale_fill_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  labs(color = "Species", fill = "Species")+
  coord_cartesian(xlim=c(-200,200),ylim=c(-100,100))
someOL_95

someOL_50 <- ggplot()+
  geom_polygon(data=SppA_50_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  geom_polygon(data=SppB_50_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.4)+
  scale_color_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  scale_fill_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  labs(color = "Species", fill = "Species")+
  coord_cartesian(xlim=c(-200,200),ylim=c(-100,100))
someOL_50

  

# Total OL -------------------------------------------------------------------
poly_path <- file.choose()
totOL_polys <- readRDS(poly_path)

SppA_95 <- subset(totOL_polys, name=="SppA_95")
SppA_50 <- subset(totOL_polys, name=="SppA_50")
SppA_95_df <- tidy(SppA_95)
SppA_50_df <- tidy(SppA_50)

SppB_95 <- subset(totOL_polys, name=="SppB_95")
SppB_50 <- subset(totOL_polys, name=="SppB_50")
SppB_95_df <- tidy(SppB_95)
SppB_50_df <- tidy(SppB_50)

col_list = c("#440154FF","#21908CFF")
col_list2 = c("#440154FF","#440154FF","#21908CFF")

totOL_95 <- ggplot()+
  geom_polygon(data=SppA_95_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  geom_polygon(data=SppB_95_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  scale_color_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  scale_fill_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  labs(color = "Species", fill = "Species")+
  coord_cartesian(xlim=c(-200,200),ylim=c(-100,100))
totOL_95

# to fix the fact that there are 2 polygons for Spp.A
SppA_50_df2 <- subset(SppA_50_df, group=="homerange2.1")
SppA_50_df3 <- subset(SppA_50_df, group=="homerange2.2")

totOL_50 <- ggplot()+
  # you're not showing the sliver because fuck the legends, is it important anyway?
  #geom_polygon(data=SppA_50_df3,aes(x=long,y=lat, col=factor(group),fill=factor(group)),show.legend = F, alpha=0.1)+
  geom_polygon(data=SppA_50_df2,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.1)+
  geom_polygon(data=SppB_50_df,aes(x=long,y=lat, col=factor(group),fill=factor(group)), alpha=0.4)+
  scale_color_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  scale_fill_manual(labels=c("Spp. A", "Spp. B"),values=col_list)+
  labs(color = "Species", fill = "Species")+
  coord_cartesian(xlim=c(-200,200),ylim=c(-100,100))
totOL_50


# plot together -----------------------------------------------------------

library(ggpubr)

combined_figure <- ggarrange(noOL_95 +rremove("xlab"), noOL_50+rremove("xlab")+rremove("ylab"),
                             someOL_95+rremove("xlab"),someOL_50+rremove("xlab")+rremove("ylab"),
                             totOL_95,totOL_50+rremove("ylab"), 
                             common.legend=T, legend="bottom",
                             labels = c("95th UD", "50th UD"),
                             ncol=2, nrow=3)
combined_figure

