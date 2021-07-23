##############################################################
################# Home Ranges and UD estimators ##############
##############################################################


# Ecology in R workshop - Russell J Gray


# Home ranges are often used with movement data to answer different
# questions about an animals behavioral and spatial ecology
# There are several common home range and Utilization Distribution (UD)
# estimators that are used by ecologists, here we will go over a 
# few of them and how to determine what is right for you

# We will start again by importing and organizing the example buffalo
# data from the last lesson to use as our data set for this excersise. 

setwd("C:/Users/Russe/Desktop/Ecology in R/Lesson 4 Movement Data")


install.packages("move")
install.packages("caTools")
install.packages("spatialEco")
install.packages("reshape2")
install.packages("tibble")

library(move)
library(adehabitatHR)
library(caTools)
library(spatialEco)
library(reshape2)
library(tibble)
library(sp)
library(ggplot2)
library(dplyr)
library(lubridate)
library(mapview)
library(cowplot)
library(ggspatial)

#--------------------------------------------

# We will use the example dataset used by Benhamou (2011) as provided
# by the adehabitatHR vignette. The move package allows us to query movement
# data from the MoveBank database, however, not all of the data is available
# as open-source, and many of them have strict useage rights that restrict
# use to exploratory analysis only. 
data(buffalo)

# --------------------------------------------------------------------------


# to make the data resemble raw data we would usually have from a telemetry
# study, we will first extract and select the data we want
dat <- buffalo[["traj"]][[1]]
# if we don't call the dplyer:: package first, we will run an ambiguous error
# this is because one of the other packages we're currently using 
dat <- dplyr::select(dat, x, y, date)

# lets change our date to a universal format, readable by R
dat$datetime <- as.POSIXct(dat$date, format = "%Y-%m-%d %hh:%mm:%ss")

# Now we can isolate the date only and add that as the date column
dat$date <- as.Date(dat$datetime, format="%Y-%m-%d")

# add the CRS projection to the data and create a spatial object
sp.UTM <- SpatialPoints(cbind(dat$x, dat$y),
                        proj4string = CRS("+init=epsg:32632"))

# reproject the spatial object to latlong
sp.latlong <- spTransform(sp.UTM,
                          CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))


#~~~~~~~~~~~~~~~~~~~~ MCP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Usually when representing the general home range of a species, 
# weather it be for land conservation or species monitoring purposes, 
# we use Minimum Convex Polygons (MCP) for our movement data

# using the adehabitatHR package, we will make an MCP for our animal
# and make the total homerange area output in km squared for the units
MCP <- mcp(sp.UTM, percent = 100, unin = "m", unout = "km2")

# Check our MCP and points
mapview(MCP, col.regions = "grey", map.type = "OpenStreetMap") + 
  mapview(sp.latlong, alpha = 0.1, cex = 0.5)

# now we will get the total area of this home range estimate and
# save it in a value object called MCP.area
MCParea <- mcp.area(sp.UTM, percent = 100, 
                    unin = "m", unout = "km2")


# ~~~~~~~~~~~~~~~~~~~~~~~~ KDE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# next we will create Kernel Density Estimator (KDE) polygons representing 
# the uitization density of the homerange. There are several options 
# available for KDE called "smoothing factors", however the most 
# commonly used are href, Least Squared Cross Validation, and manual
# plugin for the smoothing value. the polygons created by KDE
# are referred to as "isopleths" a 99% isopleth contour represents the
# total range, a 95% isopleth contour represents the most likely range
# and the 50% isopleth contour represents the "Core area" or area with
# the most abundance of points

# First we will create KDE with smoothing factor href at 99%, 95%, 
# and 50% isopleth contours using kernelUD() and getverticeshr()
Khref <- kernelUD(sp.UTM, h = 'href', grid = 800, extent = 2.2)
Khref50 <- getverticeshr(Khref, 50, unin = "m", unout = "km2")
Khref95 <- getverticeshr(Khref, 95, unin = "m", unout = "km2")
Khref99 <- getverticeshr(Khref, 99, unin = "m", unout = "km2")

# then we will get the areas of each
kernelareashref <- kernel.area(Khref, percent = c(50,95,99), 
                               unin = "m", unout = "km2")


mapview(Khref99) + mapview(Khref95) + mapview(Khref50)
# Now we will create KDE with smoothing factor LSCV at 99%, 95%, 
# and 50% isopleth contours using kernelUD() and getverticeshr()
# You will notice a warning of "non-converging algorithm". This 
# is common with LSCV when there are many points close together
# this estimator is typically meant for VHF telemetry data and not
# GPS telemetery data
KLSCV <- kernelUD(sp.UTM, h = 'LSCV', grid = 800, extent = 2.2)
KLSCV50 <- getverticeshr(KLSCV, 50, unin = "m", unout = "km2")
KLSCV95 <- getverticeshr(KLSCV, 95, unin = "m", unout = "km2")
KLSCV99 <- getverticeshr(KLSCV, 99, unin = "m", unout = "km2")

# and again, we will get the areas of each
kernel.areasLSCV <- kernel.area(KLSCV, percent = c(50,95,99), 
                               unin = "m", unout = "km2")


# Finally, we will create KDE with a manual smoothing factor  
# of 100 at 99%, 95%, and 50% isopleth contours using kernelUD() 
# and getverticeshr()
# Note: Occasionally there will be an error saying that the 
# grid/extent is too small, so we have to do a bit of tweaking
# of the parameters for it to generate all contours. This happens kind of often
# when making kernels, and there are a few common workarounds like expanding the
# grid, but none of them work 100% of the time. So it's sometimes left to
# trial and error, increasing both the grid and extent. If you use lat/lon
# as the input for your spatial objects, sometimes using UTM instead 
# mitigates the issue
Kh100 <- kernelUD(sp.UTM, h = 100, grid = 800, extent =2.2)
Kh10050 <- getverticeshr(Kh100, 50, unin = "m", unout = "km2")
Kh10095 <- getverticeshr(Kh100, 95, unin = "m", unout = "km2")
Kh10099 <- getverticeshr(Kh100, 99, unin = "m", unout = "km2")

# and get the areas of each
kernel.areash100 <- kernel.area(Kh100, percent = c(50,95,99), 
                               unin = "m", unout = "km2")


# lets map each of our isopleth sets now, and see how much of a difference
# there is between their Utilization Distribution estimates.

# check Khref - wide expansion of point range probability
mapview(Khref99, col.regions = "light grey", map.type = "OpenStreetMap") + 
  mapview(Khref95, col.regions = "grey") + 
  mapview(Khref50, col.regions = "dark grey") +
  mapview(sp.latlong, alpha = 0.1, cex = 0.5)

# check KLSCV - tight constrained range of point probability
mapview(KLSCV99, col.regions = "light grey", map.type = "OpenStreetMap") + 
  mapview(KLSCV95, col.regions = "grey") + 
  mapview(KLSCV50, col.regions = "dark grey") +
  mapview(sp.latlong, alpha = 0.1, cex = 0.5)

# check Khref - A bit less patchy, not compeltely overfit like Khref
mapview(Kh10099, col.regions = "light grey", map.type = "OpenStreetMap") + 
  mapview(Kh10095, col.regions = "grey") + 
  mapview(Kh10050, col.regions = "dark grey") +
  mapview(sp.latlong, alpha = 0.1, cex = 0.5)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~ dBBMM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# dynamic Brownian Bridge Movement Models are some of the most accurate 
# models that can be used for animal movement data The models incorporate 
# both spatial and temporal data and create movement paths and UD 
# estimates based on those paths in time while KDE only expands 
# probabiliy of usage based on independant points without temporal 
# autocorrelation being accounted for. 

# we'll start with by assigning our animal an ID, lets call it "buffalo"
dat$ID <- "buffalo"

# dBBMM also use GPS accuracy to increase the confidence output of GPS
# data. Since we don't have the accuracy data for these animals, we will
# set it to a standardized error of 5m
set_loc.error <- 5

# Remove duplicate dates, in case there are any
dat$datetime<-unique(dat$datetime, .keep_all= FALSE)


# convert data into a move object ready for dBBMM calculation
# make sure you add in the correct CRS for the data, we will be using
# our UTM data for dBBMM not lat/long data
# Calculating trajectory:
move.obj <- move(x = dat$x, y = dat$y, proj = CRS("+init=epsg:32632"), 
                 time = as.POSIXct(dat$datetime))

# review the move object
move.obj

# set the window and margin size
ws <- 13
mrg <- 5

# Set the grid that the Utilisation Distribution will be calculated accross 
# similar to what we had to do in the kernelUD function
set_grid.ext <- 2.2
set_dimsize <- 800

# Calculating the dBBMM 
dbbmm <- brownian.bridge.dyn(object = move.obj, # move object
                             location.error = set_loc.error, # locations error
                             margin = mrg, window.size = ws,  # window and margin
                             ext = set_grid.ext, dimSize = set_dimsize,  # grid
                             verbose = F) 

# The dBBMM object is created as a raster, we will have to run the next few lines 
# of code to turn it into a UD object
dbbmm.sp <- as(dbbmm, "SpatialPixelsDataFrame")
dbbmm.sp.ud <- new("estUD", dbbmm.sp)
dbbmm.sp.ud@vol = FALSE
dbbmm.sp.ud@h$meth = "dBBMM"
dbbmm.ud <- getvolumeUD(dbbmm.sp.ud, standardize = TRUE) 


# we can now pull out isopleth contours
dBBMM.050 <- getverticeshr(dbbmm.ud, percent = 50, unin = "m", unout = "km2")
dBBMM.095 <- getverticeshr(dbbmm.ud, percent = 95, unin = "m", unout = "km2")
dBBMM.099 <- getverticeshr(dbbmm.ud, percent = 99, unin = "m", unout = "km2")


# and Finally we will create a value object of the home range area
dBBMM.areas <- kernel.area(dbbmm.ud, percent = c(50,95,99),
                           unin = "m", unout = "km2")

# Now we can map the ispopleths to see what the dBBMM created
mapview(dBBMM.099, col.regions = "light grey", map.type = "OpenStreetMap") + 
  mapview(dBBMM.095, col.regions = "grey") + 
  mapview(dBBMM.050, col.regions = "dark grey") +
  mapview(sp.UTM, alpha = 0.1, cex = 0.5)

# As we can see, the dBBMM shows much more connectivity and less patchy areas 
# than the KDE isopleths.


# ~~~~~~~~~~~~~~~~~~~~~ MODEL VALIDATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# one of the most commonly used metrics for relative goodness-of-fit of
# Utilization Distribution estimators is the Area Under Curve (AUC), just
# like the validation used for Species Distribution Modelling. 
# Another important technique used to examine UD isopleth contours is 
# complexity of the polygons. Typically the best model will not only have
# the *HIGHEST* AUC value, but will also have a tradeoff between type I errors 
# (oversmoothing - UD representing areas that aren't occupied) and Type II 
# errors (undersmoothing - UD lacking areas that may be occupied), represented
# by the complexity of the isopleth contours

# First we will make our AUC table. We do this by making a raster of 
# our locations over our Utilization Distribution object. Then, we pull 
# the actual vs. null values from the dataframe, and with that we 
# calculate the AUC using wilcoxon and ROC methods. However we will
# remove the KLSCV value from this analysis since it failed to converge
# and is therefore a NULL value. Also, since MCP are just a generalized 
# home range representative and not a UD estimator, it will also be ommited
# from these analyses

#########this is for dBBMM
nlocrast <- count.points(sp.UTM,dbbmm.ud) 
kerneldata <- dbbmm.ud@data$n # vector containing volume contour (= predicted) values
pointdata <- nlocrast@data$x 
pointdata <- ifelse(pointdata>=1,1,0) # vector containing location (= actual) values 
AUCdBBMM<-colAUC(kerneldata, pointdata, plotROC=FALSE, alg=c("Wilcoxon","ROC"))



#########this is for Khref
Khrefud <- getvolumeUD(Khref, standardize=FALSE)
nlocrast<-count.points(sp.UTM,Khrefud) 
kerneldata <- Khrefud@data$n # vector containing volume contour (= predicted) values
pointdata <- nlocrast@data$x 
pointdata <- ifelse(pointdata>=1,1,0) # vector containing location (= actual) values 
AUCKhref <- colAUC(kerneldata, pointdata, plotROC=FALSE, alg=c("Wilcoxon","ROC"))

#########this is for KLSCV
AUCKLSCV <- NA

#########this is for Kh100 (remember we used UTM for this estimator)
Kh100ud <- getvolumeUD(Kh100, standardize=FALSE)
nlocrast<-count.points(sp.UTM,Kh100ud) 
kerneldata <- Kh100ud@data$n # vector containing volume contour (= predicted) values
pointdata <- nlocrast@data$x 
pointdata <- ifelse(pointdata>=1,1,0) # vector containing location (= actual) values 
AUCKh100 <- colAUC(kerneldata, pointdata, plotROC=FALSE, alg=c("Wilcoxon","ROC"))


# Now we will create a data frame with our data
AUCdf <- cbind(AUCKhref[1], AUCKLSCV[1], AUCKh100[1],
               AUCdBBMM[1])
# convert it to a dataframe so it can have columns and row names
AUCdf <- as.data.frame(AUCdf)
# name the columns with the corresponding estimation method
names(AUCdf) <- c("Khref","KLSCV","Kh100","dBBMM")
# print and view the AUC table
print(AUCdf)

# As we can see from our AUC evaluation the Kh100 appears to be 
# the best fit model, with dBBMM in second and the Khref last (which is
# usually the case... basically, don't use khref). However, to make sure
# we choose the best-of-fit model we should also take into consideration
# the complexity and liklihood of type I and type II errors at each
# isopleth contour. To do this, we need to divide the perimeter of each 
# polygon by its total area (complexity = perimeter/total area)

# first, we will get the perimeter of each isopleth contour using the 
# polyPerimeter() function from the spatialEco package


#dbbmm perimeter
pdbbmm99<-as(dBBMM.099, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000
pdbbmm95<-as(dBBMM.095, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000
pdbbmm50<-as(dBBMM.050, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000

#KDEhref perimeter
pKhref99<-as(Khref99, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000
pKhref95<-as(Khref95, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000
pKhref50<-as(Khref50, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000

#KDELSCV perimeter
pKLSCV99<-as(KLSCV99, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000
pKLSCV95<-as(KLSCV95, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000
pKLSCV50<-as(KLSCV50, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000

#Kh100 perimeter
pKh10099<-as(Kh10099, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000
pKh10095<-as(Kh10095, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000
pKh10050<-as(Kh10050, "sf") %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length()/1000



#-------------------------------------

# Now we will get the information for the total area of each polygon 

#dbbmm area 
adbbmm99<-dBBMM.areas[3]
adbbmm95<-dBBMM.areas[2]
adbbmm50<-dBBMM.areas[1]

#KDEhref area
aKhref99<-kernelareashref[3]
aKhref95<-kernelareashref[2]
aKhref50<-kernelareashref[1]

#KDELSCV area
aKLSCV99<-kernel.areasLSCV[3]
aKLSCV95<-kernel.areasLSCV[2]
aKLSCV50<-kernel.areasLSCV[1]

#Kh100 area
aKh10099<-kernel.areash100[3]
aKh10095<-kernel.areash100[2]
aKh10050<-kernel.areash100[1]

#---------------------------------------------
# Now we calculate the complexity
# Complexity = perimeter(m) / total area (km)

#dbbmm complexity
cdbbmm99<-pdbbmm99 / adbbmm99
cdbbmm95<-pdbbmm95 / adbbmm95
cdbbmm50<-pdbbmm50 / adbbmm50

#KDEhref complexity
cKhref99<-pKhref99 / aKhref99
cKhref95<-pKhref95 / aKhref95
cKhref50<-pKhref50 / aKhref50

#KDELSCV complexity
cKLSCV99<-pKLSCV99 / aKLSCV99
cKLSCV95<-pKLSCV95 / aKLSCV95
cKLSCV50<-pKLSCV50 / aKLSCV50

#Kh100 complexity
cKh10099<-pKh10099 / aKh10099
cKh10095<-pKh10095 / aKh10095
cKh10050<-pKh10050 / aKh10050


# Save the complexity values as a data frame
Complexity <- cbind(cdbbmm99,cdbbmm95,cdbbmm50,
                    cKhref99,cKhref95,cKhref50,
                    cKLSCV99,cKLSCV95,cKLSCV50,
                    cKh10099,cKh10095,cKh10050)
# convert it to a dataframe so it can have columns and row names
Complexity <- as.data.frame(Complexity)
# name the columns with the corresponding estimation method
names(Complexity) <- c("dbbmm99","dbbmm95","dbbmm50",
                       "Khref99","Khref95","Khref50",
                       "KLSCV99","KLSCV95","KLSCV50",
                       "Kh10099","Kh10095","Kh10050")
# rename rows for the sake of neatness
row.names(Complexity) <- "Complexity"

Complexity

# now lets create a more meaningful data frame to review our data
# from both analyses
model.analysis <- data.frame(1:4)
model.analysis$Method <- c("Khref", "KLSCV", "Kh100", "dBBMM")
model.analysis$AUC <- c(AUCdf[1,])
model.analysis$comp_99 <- c(Complexity$Khref99, Complexity$KLSCV99,
                            Complexity$Kh10099, Complexity$dbbmm99) 
model.analysis$comp_95 <- c(Complexity$Khref95, Complexity$KLSCV95,
                            Complexity$Kh10095, Complexity$dbbmm95) 
model.analysis$comp_50 <- c(Complexity$Khref50, Complexity$KLSCV50,
                            Complexity$Kh10050, Complexity$dbbmm50) 


# To properly examine our analysis, we have to understand that large 
# complexity values indicate oversmoothing, while smaller values may
# indicate undersmoothing. First we will run code to disable scientific 
# notation so we can see the values more clearly
options(scipen=999)
print(model.analysis)

# Now we can save our model analysis as a .csv file
# first we need to switch AUC to a numeric value
model.analysis$AUC <- as.numeric(model.analysis$AUC)
model.analysis
write.csv(model.analysis, "model_analysis.csv")

# We can see from our values that Khref has large complexity 
# and is likely to be oversmoothing all isopleth contours, we can visualize
# this oversmoothing by plotting our MCP and 99% khref polygon
# and see that it massively breaches the maximum 100% extent of our points
mapview(MCP, col.regions = "red") + mapview(Khref99)

# we can better visualize the trends by making a ggplot graph
# we first need to "melt()" the dataframe to condense the data
melt.analysis <- melt(model.analysis, id = c("Method", "AUC"))
#remove the first 4 rows
melt.analysis <- melt.analysis[-c(1,2,3,4),]

# now plot the figure
analysis.plot <- ggplot(melt.analysis, aes(x = variable, y = value)) +
  geom_line(aes(color = Method, group = Method), size=0.8)+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold.italic"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"))+
  xlab("Countour (%)") + ylab("Complexity")+
  labs(title = "Model Complexity at 50%, 95%, & 99%", size = 10)+
  scale_x_discrete(labels=c("99%", "95%", "50%"))

analysis.plot

# We can also see that the dBBMM model has a slightly higher complexity 
# for the 99% isopleth than kLSCV and Kh100, this is likely because the 
# model creates a temporal path with the polygons and should not necessarily
# be considered as overfitting the data.

# dBBMM appears to be the tradeoff between the models other than the heightened
# value at 99%, as its values fall between both KLSCV and Kh100 for the lower
# isopleth contours.


# ~~~~~~~~~~~~~~~~~~~~~ Figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ggplot is weird with UTM spatial objects sometimes, so we will transform
# all of our spatialpolygons to latlong first
MCP <- spTransform(MCP,
                   CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
#-------------------------------
Khref99 <- spTransform(Khref99,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
Khref95 <- spTransform(Khref95,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
Khref50 <- spTransform(Khref50,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
#-----------------------------
KLSCV99 <- spTransform(KLSCV99,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
KLSCV95 <- spTransform(KLSCV95,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
KLSCV50 <- spTransform(KLSCV50,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
#------------------------------
Kh10099 <- spTransform(Kh10099,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
Kh10095 <- spTransform(Kh10095,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
Kh10050 <- spTransform(Kh10050,
                       CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
#------------------------------
dBBMM.099 <- spTransform(dBBMM.099,
                        CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
dBBMM.095 <- spTransform(dBBMM.095,
                         CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
dBBMM.050 <- spTransform(dBBMM.050,
                         CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))



#dBBMM
dBBMM.plot <- ggplot() +
  geom_spatial_polygon(data = MCP, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 2, color = "black", fill = NA) +
  geom_spatial_polygon(data = dBBMM.050, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 1, color = "black", fill="black") +
  geom_spatial_polygon(data = dBBMM.095, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 2, color = "black", fill="black") +
  geom_spatial_polygon(data = dBBMM.099, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 3, color = "black", fill="black") +
  geom_point(data = as.data.frame(sp.latlong),
             aes(x = coords.x1, y = coords.x2),
             colour = "black", pch = 16, size = 0.5, alpha=0.2) + 
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", title = "dBBMM")+
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text.y = element_text(face = 2, size = 9, colour = "black",
                                   margin = margin(t = 0, r = 3, b = 0, l = 10)),
        axis.text.x = element_text(face = 2, size = 9, colour = "black",
                                   margin = margin(t = 3, r = 0, b = 5, l = 0)),
        axis.title.y = element_text(face = 2),
        axis.title.x = element_text(face = 2),
        axis.ticks = element_line(colour = "black"))
dBBMM.plot

#-----------------------------------------------------

#Khref
Khref.plot <- ggplot() +
  geom_spatial_polygon(data = MCP, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 2, color = "black", fill = NA) +
  geom_spatial_polygon(data = Khref99, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 1, color = "black") +
  geom_spatial_polygon(data = Khref95, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 2, color = "black") +
  geom_spatial_polygon(data = Khref50, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 3, color = "black") +
  geom_point(data = as.data.frame(sp.latlong),
             aes(x = coords.x1, y = coords.x2),
             colour = "black", pch = 16, size = 0.5, alpha=0.2) + 
  theme_bw() +
  labs(x = "Longitude", y = "Latitude",title = "Khref")+
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text.y = element_text(face = 2, size = 9, colour = "black",
                                   margin = margin(t = 0, r = 3, b = 0, l = 10)),
        axis.text.x = element_text(face = 2, size = 9, colour = "black",
                                   margin = margin(t = 3, r = 0, b = 5, l = 0)),
        axis.title.y = element_text(face = 2),
        axis.title.x = element_text(face = 2),
        axis.ticks = element_line(colour = "black"))
Khref.plot 


#--------------------------------------------------------------

#KLSCV
KLSCV.plot <- ggplot() +
  geom_spatial_polygon(data = MCP, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 2, color = "black", fill = NA) +
  geom_spatial_polygon(data = KLSCV99, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 1, color = "black") +
  geom_spatial_polygon(data = KLSCV95, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 2, color = "black") +
  geom_spatial_polygon(data = KLSCV50, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 3, color = "black") +
  geom_point(data = as.data.frame(sp.latlong),
             aes(x = coords.x1, y = coords.x2),
             colour = "black", pch = 16, size = 0.5, alpha=0.2) + 
  theme_bw() +
  labs(x = "Longitude", y = "Latitude",title = "KLSCV")+
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text.y = element_text(face = 2, size = 9, colour = "black",
                                   margin = margin(t = 0, r = 3, b = 0, l = 10)),
        axis.text.x = element_text(face = 2, size = 9, colour = "black",
                                   margin = margin(t = 3, r = 0, b = 5, l = 0)),
        axis.title.y = element_text(face = 2),
        axis.title.x = element_text(face = 2),
        axis.ticks = element_line(colour = "black"))
KLSCV.plot 

#--------------------------------------------------------------
#KLSCV
Kh100.plot <- ggplot() +
  geom_spatial_polygon(data = MCP, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 2, color = "black", fill = NA) +
  geom_spatial_polygon(data = Kh10099, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 1, color = "black") +
  geom_spatial_polygon(data = Kh10095, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 2, color = "black") +
  geom_spatial_polygon(data = Kh10050, aes(x = long, y = lat, group = group), size = 0.5,
                       alpha = 0.2, linetype = 3, color = "black") +
  geom_point(data = as.data.frame(sp.latlong),
             aes(x = coords.x1, y = coords.x2),
             colour = "black", pch = 16, size = 0.5, alpha=0.2) + 
  theme_bw() +
  labs(x = "Longitude", y = "Latitude",title = "Kh100")+
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text.y = element_text(face = 2, size = 9, colour = "black",
                                   margin = margin(t = 0, r = 3, b = 0, l = 10)),
        axis.text.x = element_text(face = 2, size = 9, colour = "black",
                                   margin = margin(t = 3, r = 0, b = 5, l = 0)),
        axis.title.y = element_text(face = 2),
        axis.title.x = element_text(face = 2),
        axis.ticks = element_line(colour = "black"))
Kh100.plot 

#-------------------------------------------
#group all plots together in a 2x2 grid
pg <- plot_grid(Khref.plot, KLSCV.plot, Kh100.plot, dBBMM.plot,
              nrow = 2, ncol = 2)
pg


# Now plot all all figures with the complexity figure
pg.1 <- plot_grid(pg, analysis.plot, nrow = 1, ncol = 2)
pg.1

model.analysis
#~~~~~~~~~~~~~~~~ HOME RANGE VALUES ~~~~~~~~~~~~~~~~~~~~~~~

# In order to visualize the variations in home range values,
# we will plot them using a barchart in ggplot

# first we will need to create and melt a data frame with the 
homeranges <- as.data.frame(cbind(dBBMM.areas, kernelareashref,
                    kernel.areasLSCV, kernel.areash100))
homeranges <- tibble::rownames_to_column(homeranges, "VALUE")

# melt for plotting
homeranges.melt <- melt(homeranges, id = "VALUE")

#plot the data
p3 <- ggplot(homeranges.melt, aes(x=variable, y=value, fill=VALUE)) +
  geom_bar(stat="identity", position = "dodge", col = "black") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold.italic"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"))+
  xlab("Model") + ylab("Home Range")+
  labs(title = "Home Range size at 50%, 95%, & 99% (km^2)", size = 10,
       subtitle = "")+
  scale_x_discrete(labels=c("dBBMM", "Khref", "KLSCV", "Kh100"))

p3


# combine all plots
pg2 <- plot_grid(analysis.plot, p3, nrow = 2, ncol = 1) #combine analysis and home ranges
pg3 <- plot_grid(pg, pg2, nrow = 1, ncol = 2) # combine all plots
pg3 # print plots (use zoom function to expand the window)




# END




