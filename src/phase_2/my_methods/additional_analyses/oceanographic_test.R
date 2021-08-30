
# Stack all wind datasets.
# Each stack will be associated with a date and hourly timestamp (24 per day)
nc_dir<-'C:/Users/dalla/Downloads/'
setwd(nc_dir)

wind_files <- list.files(pattern='*.nc') 

wind_t1 <- read_ncdf(wind_files[3])
wind_t1 # Make sure you got the right stuff!
times_t1<-st_get_dimension_values(wind_t1, "time")  # stores times from each file 
wind_t1_u<-as(wind_t1[1,,,], "Raster") # refer to U component
wind_t1_v<-as(wind_t1[2,,,], "Raster") # v component 


wind_t2 <- read_ncdf(wind_files[2])
wind_t2
times_t2<-st_get_dimension_values(wind_t2, "time")  
wind_t2_u<-as(wind_t2[1,,,], "Raster")
wind_t2_v<-as(wind_t2[2,,,], "Raster")


# Stack rasters  ....................
u_stack<-stack(wind_t1_u, wind_t2_u)
v_stack<-stack(wind_t1_v, wind_t2_v)

all_times<-rbind(as.tibble(times_t1), as.tibble(times_t2))

