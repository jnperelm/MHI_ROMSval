library(tidyverse)
library(raster)
library(pracma)
library(reshape2)

# Import matched ROMS-STR data

setwd("Matched_STR_ROMS")
df <- list.files(pattern = ".rds") %>%
  map(readRDS) %>% 
  bind_rows()

df <- df %>% drop_na(ROMS_temp)
# Create month and year columns
df<- df %>%
  dplyr::mutate(year = lubridate::year(TIMESTAMP_3h), 
                month = lubridate::month(TIMESTAMP_3h))



# Import bathymetry netCDF
setwd("C:/Users/jessica.perelman/Documents/EDS_ROMS_TempProj/ROMS_subtemp_validation/ROMSval_Rproject")
bathy = "Env_data/bathy_Hawaii.nc"

bathy_raster <- brick(bathy, varname ="altitude")
bathy_raster; class(bathy_raster)
rug_raster = aggregate(bathy_raster, fact = 2, fun = 'sd') # fact = 2 means calculate sd from immediately adjoining cells
rug_raster; class(rug_raster)

rug_df = as.data.frame(rug_raster, xy=TRUE)
names(rug_df) <- c("LONGITUDE","LATITUDE","rugosity")

# Match values via 2D interpolation using a for loop for dates

rugosity <- acast(rug_df,LATITUDE~LONGITUDE,value.var='rugosity')
lon <- as.numeric(colnames(rugosity))
lat <- as.numeric(rownames(rugosity))
df$rugosity <- interp2(lon, lat, rugosity, xp=df$LONGITUDE, yp=df$LATITUDE, method = "nearest")

# Regrid bathymetry from 0.01 degreee to 0.1 degree

# resample grid
higher_res_grid = raster(ncol = 70, nrow = 50)  # create new raster grid with 0.1x0.1 grid size
bb = extent(-161, -154, 18, 23) # set lat lon extent
extent(higher_res_grid) <- bb # set lat lon extent
# resample
bathy_resampled <- resample(bathy_raster, higher_res_grid, method = "bilinear")
plot(bathy_resampled$layer)
rug_resampled = aggregate(bathy_resampled, fact = 2, fun = 'sd') # fact = 2 means calculate sd from immediately adjoining cells
rug_resampled; class(rug_resampled)

#convert back to dataframe
rug_resampled_df = as.data.frame(rug_resampled, xy=TRUE)
names(rug_resampled_df) <- c("LONGITUDE","LATITUDE","rugosity_0.1")

# Match values via 2D interpolation using a for loop for dates

rugosity_0.1 <- acast(rug_resampled_df,LATITUDE~LONGITUDE,value.var='rugosity_0.1')
lon <- as.numeric(colnames(rugosity_0.1))
lat <- as.numeric(rownames(rugosity_0.1))
df$rugosity_0.1 <- interp2(lon, lat, rugosity_0.1, xp=df$LONGITUDE, yp=df$LATITUDE, method = "nearest")

##########################
# Match Ocean Nino Index #
##########################

library(rsoi)
df$month2 = month.abb[df$month]
oni <- download_oni()
names(oni)[names(oni) == 'Month'] <- 'month2'
names(oni)[names(oni) == 'Year'] <- 'year'
df <- merge(df, oni[,c(1,2,5)], by = c("month2", "year"))

###############################
# Calculate distance to shore #
###############################
library(sf)
library(geosphere)
library(osmdata)

# convert data to sf object
df_sf <- df %>% st_as_sf(coords = c('LONGITUDE','LATITUDE')) %>% 
  st_set_crs(4326)
df_sf = unique(df_sf[,c(4,6,10,17)])

# get bounding box for osm data download and download coastline data for this area
osm_box <- getbb(place_name = "Hawaii") %>%
  opq () %>% 
  add_osm_feature("natural", "coastline") %>% 
  osmdata_sf() 

# use dist2Line from geosphere - only works for WGS84 data
dist <- geosphere::dist2Line(p = st_coordinates(df_sf), 
                             line = st_coordinates(osm_box$osm_lines)[,1:2])
dist = data.frame(dist)
names(dist)[c(2,3)] <- c("LONGITUDE", "LATITUDE")

#combine initial data with distance to coastline
df_new <- cbind(df_sf,dist)
df_new = left_join(df, unique(df_new[-25,c(1,3,4)]), by = c("DEPTH", "OCC_SITEID")) # remove duplicate Oahu site

# Save new file as metadata
saveRDS(df_new[,-c(1,20)], file = paste("Matched_STR_ROMS/Combined/STR_ROMS_matched.rds", sep = ""))


