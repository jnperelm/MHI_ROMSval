# Script to extract ROMS temperature data and match with in situ subsurface data
# Uses rerddapXtracto package for 3D variables

rm(list = ls())

library(tidyverse)
library(lubridate)
library(rerddap)
library(rerddapXtracto)
library(reshape2)

# Create list of STR files and file basenames
data = list.files(path = "STR_MHI_2013_2019/Recovered2019", pattern = "*.rds")
basenames = gsub("\\.rds$","", list.files(path = "STR_MHI_2013_2019/Recovered2019", pattern="\\.rds$"))

for(i in 1:length(data)) {
  
  # i = 5
  
  # Import dataset
  path = "STR_MHI_2013_2019/Recovered2019/"
  df <- readRDS(paste0(path,data[i]))
  df$TIMESTAMP = as.POSIXct(df$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  head(df)
  df = subset(df, TIMESTAMP > "2007-07-11 21:00:00")
  df$TIMESTAMP = round_date(df$TIMESTAMP,unit="1 minute")
  
  # Average observations from +/- 10 minutes (or +/- 1 hr) around every 3-hr time point
  numeric <- df %>% 
    mutate(TIMESTAMP_3h = as.POSIXct(ifelse(TIMESTAMP <= (floor_date(TIMESTAMP, unit = "3hours") + 600), # +/- 10 mins on either side of every 3 hours 
    #mutate(TIMESTAMP_3h = as.POSIXct(ifelse(TIMESTAMP <= (floor_date(TIMESTAMP, unit = "3hours") + 3600), # +/- 1 hr on either side of every 3 hours
                                            floor_date(TIMESTAMP, unit = "3hours"),
                                            ifelse(TIMESTAMP >= (ceiling_date(TIMESTAMP, unit = "3hours") - 3600),
                                                   ceiling_date(TIMESTAMP, unit = "3hours"),NA)), 
                                     origin = "1970-01-01", tz = "UTC"))  %>% 
    subset(!is.na(TIMESTAMP_3h)) %>%
    group_by(TIMESTAMP_3h) %>% 
    summarise_at(vars(INSTRUMENTSN, MISSIONINSTRUMENTID, LATITUDE, LONGITUDE, DEPTH,
                      TEMP_C), mean)
  
  character <- df %>% 
    mutate(TIMESTAMP_3h = as.POSIXct(ifelse(TIMESTAMP <= (floor_date(TIMESTAMP, unit = "3hours") + 600), # +/- 2.5 mins on either side of every 3 hours
                                            floor_date(TIMESTAMP, unit = "3hours"),
                                            ifelse(TIMESTAMP >= (ceiling_date(TIMESTAMP, unit = "3hours") - 600),
                                                   ceiling_date(TIMESTAMP, unit = "3hours"),NA)),
                                     origin = "1970-01-01", tz = "UTC"))  %>% 
    subset(!is.na(TIMESTAMP_3h)) %>%
    group_by(TIMESTAMP_3h) %>% 
    summarise_at(vars(OCC_SITEID, REGION, LOCATION, INSTRUMENTTYPE), ~ paste(.))
  
  df_3h <- inner_join(unique(character), numeric, by = c("TIMESTAMP_3h"))
  df_3h <- df_3h[df_3h$TIMESTAMP_3h <= as.POSIXct("2017-05-04 21:00:00", tz = "UTC"),]
  
  # Pull metadata from ERDDAP server and set lon, lat, depth, and time ranges to extract
  dataInfo <- info('roms_hiig_reanalysis')
  parameter <- 'temp'
  xcoord <- c(range(df_3h$LONGITUDE))
  ycoord <- c(range(df_3h$LATITUDE))
  zcoord <- c(range(df_3h$DEPTH))
  tcoord <- c(range(df_3h$TIMESTAMP_3h))
  
  # Extract variables from ERDDAP that match survey data range
  
  library(rerddapXtracto)
  extract = rxtracto_3D(dataInfo, parameter, xcoord = xcoord, ycoord = ycoord,
                        zcoord = zcoord, tcoord = tcoord, zName = "depth")
  new = extract$temp
  dimnames(new) <- list(extract$longitude, extract$latitude, extract$depth, as.character(extract$time))
  
  new2 = melt(new, varnames = c("lon", "lat", "depth", "TIMESTAMP_3h"), value.name = "ROMS_temp")
  new2$TIMESTAMP_3h = as.POSIXct(new2$TIMESTAMP_3h, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  head(new2)
  
  # Merge Extracted ROMS data with STR data
  
  df_new = merge(df_3h, new2[,c(4,5)], by = "TIMESTAMP_3h")
  
  # Save matched dataset as new .csv file
  
  write.csv(df_new, paste("Matched_STR_ROMS/",
                          basenames[i], "-ROMS.csv", sep = ""), row.names = FALSE)
}

