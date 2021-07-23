# normalize UDs
# provided by Melinda, need to adapt to make work for my work
# April 2021

s <- sum(values(grid_ct_sppi), na.rm=TRUE)
values(grid_ct_sppi) <- values(grid_ct_sppi)/s
df <- data.frame(v=values(grid_ct_sppi), id=seq(1,ncell(grid_ct_sppi)))
df.2 <- df[rev(order(df$v)),]
df.2$cs <- cumsum(df.2[,‘v’])

# (grid_ct_sppi is density raster)

#CREATE BREAKS
df.2$c<-cut(df.2$cs,breaks=seq(0,1,by=0.05),labels=seq(100,5,by=-5))
df.3 <- df.2[order(df.2$id),]
values(grid_ct_sppi) <- as.numeric(df.3$c)
plot(grid_ct_sppi)

# basically breaks density by 5, so then you can call for 95 or 50 or etc
# (well density goes from 0 to 1, by 0.05, then grid is labeled 0-100 for density)
# 
# So, first code snippet normalizes density from 0-1 and then second snippet
# breaks it into bins labeled 0-100 by 5, so you can call different density grids 
# (i hope im explaining that in a way thats understandable!)