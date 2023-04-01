rm(list = ls())

library(tidyverse)
library(raster)
library(pracma)
library(reshape2)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(RColorBrewer)
library(viridis)
library(scales)

# Import rainfall netCDF

precip_path = "C:/Users/jessica.perelman/Documents/GitHub/env_data_summary/Data_Download/Precipitation_CHIRPS_Daily/Island_Level_Data/"
files = list.files(precip_path,pattern='*.nc', full.names = TRUE)

list1 = NULL
list2 = NULL
for (i in 1:length(files)) {
  # i = 4
  precip_raster <- brick(files[i], varname ="precip")
  precip_raster; class(precip_raster)
  #m_indices <- format(as.Date(names(precip_raster), format = "X%Y.%m.%d"), format = "%m")
  y_indices <- format(as.Date(names(precip_raster), format = "X%Y.%m.%d"), format = "%Y")
  m_indices <- format(as.Date(names(precip_raster), format = "X%Y.%m.%d"), format = "%m")
  #precip_raster = stackApply(precip_raster, y_indices, fun = 'mean')
  precip_raster2 = stackApply(precip_raster, m_indices, fun = 'mean')
  df = as.data.frame(precip_raster, xy=TRUE)
  df2 = as.data.frame(precip_raster2, xy=TRUE)
  #list1[[i]] <- df
  list2[[i]] <- df2
}
#list[[2]] <- list[[2]][,-3]
precip_df <- do.call(rbind, list2)
precip_df <- precip_df %>% mutate_all(~ifelse(is.nan(.), NA, .))
precip_df = precip_df[,c(1:8,10:14,9)] # reorder

# plot by year/month
setwd("C:/Users/jessica.perelman/Documents/EDS_ROMS_TempProj/ROMS_subtemp_validation/Env_data/rainfall_annual")
for (i in 3:length(names(precip_df))) {
  
  # i = 3
  
  tiff(file=paste(names(precip_df)[i],".tif"), res=300, height=80, width=170, units="mm")
  x <- ggplot() + geom_tile(data=precip_df, aes(x=x,y=y,fill=precip_df[,i]), height=0.05,width=0.05) +
    scale_fill_viridis(option = "D", na.value="white", lim = c(0.01, 18.59), name = "Precipitation (mm/day)") +
    ggtitle(names(precip_df)[i]) +
    ylab("Latitude") + xlab("Longitude") + 
    theme_bw() + scale_y_continuous(breaks = c(19,20,21,22)) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 14),
          axis.title = element_text(size= 14), legend.title = element_text(size=14),
          legend.text = element_text(size=12))
  
  print(x)
  dev.off()
}

library(plyr)
mean_mon <- data.frame(colMeans(precip_df[,c(3:14)], na.rm=T))
mean_mon$month <- c(1:12)
colnames(mean_mon)[1] <- "precip"
ggplot(data=mean_ann) + geom_line(aes(x=month, y=precip, group=1)) + theme_bw() + theme(panel.grid = element_blank())
