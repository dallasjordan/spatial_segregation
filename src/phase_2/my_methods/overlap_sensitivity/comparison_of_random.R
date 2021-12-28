# in order to plot polygons, first fortify the data
all_data@data$id <- rownames(all_data@data)
# create a data.frame from our spatial object
all_data_df <- fortify(all_data, region = "id")
# merge the "fortified" data with the data from our spatial object
alldatadf <- merge(all_data_df, all_data@data,
                   by = "id")

ggplot(data = alldatadf, aes(x = long, y = lat, group = group, fill="red")) +
  geom_polygon()

plot(all_data)

ggplot(all_data, aes(x = long, y = lat, group = "name")) + 
  geom_polygon(colour='black', fill='white')


ggplot(all_data, aes(x = x, y = y)) +
  geom_polygon(aes(fill = "name", group = id))

SppA_95 <- subset(all_data, name=="SppA_95")
SppA_50 <- subset(all_data, name=="SppA_50")

plot(all_data)
plot(SppA_95)
plot(SppA_50)

]]
