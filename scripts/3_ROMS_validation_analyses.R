###################################################
# ROMS subsurface temperature validation analyses #
# JN Perelman, KR Tanaka, J Smith, H Barkley      #
###################################################
rm(list = ls())

library(marmap)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(RColorBrewer)
library(viridis)
library(scales)
library(zoo)
library(lubridate)
library(openair)
library(mgcv)
library(raster)
library(reshape2)
library(pracma)
library(tidync)
library(metR)
library(ggmap)

# Import matched temperature data

df <- read.csv("Matched_STR_ROMS/Combined/STR_ROMS_matched.csv")
df <- df %>% drop_na(ROMS_temp)
df$TIMESTAMP_3h <- as.POSIXct(df$TIMESTAMP_3h, tz="UTC")
df = df[order(df$TIMESTAMP_3h),]

# remove Kaneohe Bay sites-- these are masked in the MHI ROMS model
df <- df %>% subset(OCC_SITEID != "OCC-OAH-004" &
                      OCC_SITEID != "OCC-OAH-005" &
                      OCC_SITEID != "OCC-OAH-006")

# Create month and year columns
df<- df %>%
  dplyr::mutate(year = lubridate::year(TIMESTAMP_3h), 
                month = lubridate::month(TIMESTAMP_3h),
                day = lubridate::day(TIMESTAMP_3h))

######################
# Plot STR locations #
######################

df_map <- df %>% group_by(OCC_SITEID) %>% summarise_at(vars(LONGITUDE, LATITUDE, DEPTH), list(~mean(., na.rm=T)))
df_map = merge(df_map, unique(df[,c(3,5)]), by = "OCC_SITEID")
df_map$depth_discrete = 
  ifelse(df_map$DEPTH < 10, 5,
         ifelse(df_map$DEPTH > 10 & df_map$DEPTH < 20, 15,
                ifelse(df_map$DEPTH > 20, 25,NA)))


#reading in a shapefile with multiple polygons
library(rgdal)
hawaii <-readOGR("Env_data/coastline_hawaii/Coastline.shp")
hawaii <- spTransform(hawaii, "+init=epsg:4326") # need to transform lat and long coordinates
set.seed(13)
ggplot() + geom_polygon(data = hawaii, aes(x=long, y=lat, fill = group), color = "black") + 
  scale_fill_manual(values = rep("lightgray",14)) +
  geom_jitter(data = df_map, aes(x=LONGITUDE, y=LATITUDE, color = factor(depth_discrete)), size= 3.4, width = 0.038, height = 0.038) +
  scale_color_manual(values = c("#66CCFF","#3366CC","#003399")) +
  ylab("Latitude") + xlab("Longitude") + 
  coord_sf(xlim = c(-160.1, -155), ylim = c(18.5, 22.3)) +
  theme_bw() + scale_y_continuous(breaks = c(19,20,21,22)) +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 14),
                     axis.title = element_text(size= 14), legend.title = element_text(size=14),
                     legend.text = element_text(size=12))

unique(df$OCC_SITEID)

# map by island to see individual locations
Island_Extents <- subset(read.csv("~/EDS_ROMS_TempProj/ROMS_subtemp_validation/Island_Extents.csv"),
                                  region == "MHI")
names(Island_Extents)[names(Island_Extents) == 'island.code'] <- 'LOCATION'
df_map = merge(df_map, Island_Extents[,c(2:6)], by = "LOCATION")

df_map %>%
  group_by(LOCATION) %>% 
  ggplot() +
  annotation_map(map_data("world")) +
  geom_point(data = df_map, aes(LONGITUDE, LATITUDE, fill = LOCATION), shape = 21, alpha = 0.5, size=3) +
  geom_rect(mapping = aes(
    xmin = x_min,
    xmax = x_max,
    ymin = y_min,
    ymax = y_max,
    fill = LOCATION,
    color = LOCATION), alpha = 0.01) + theme_bw() +
  facet_wrap(~LOCATION, scales = "free") +
  ylab("Latitude") + xlab("Longitude") +
  theme(legend.position = "none") 

###################################
# Plot raw time series (mean)     #
# Deming Regression               #
###################################

# Plot ROMS and STR time series
# plot time series
df %>%
  mutate(yrmoday = as.Date(paste0(df$year, "-", df$month, "-", df$day), "%Y-%m-%d")) %>%
  group_by(yrmoday) %>%
  summarise(ROMS = mean(ROMS_temp, na.rm=T),
            STR = mean(TEMP_C, na.rm=T)) %>%
  mutate(month = month(yrmoday)) %>%
  ggplot() +
    geom_line(aes(x=yrmoday, y=ROMS), color="#00CCCC", lwd=0.6) + ylim(c(23,29)) +
    #geom_point(aes(x=yrmoday, y=STR), size = 0.9, shape = 21) +
    #geom_errorbar(aes(x=yrmo, ymin=STR-STR_err, ymax=STR+STR_err)) +
    geom_line(aes(x=yrmoday, y=STR)) + 
    #geom_point(aes(x=yrmoday, y=ROMS), size = 0.9) +
    #geom_errorbar(aes(x=yrmo, ymin=ROMS-ROMS_err, ymax=ROMS+ROMS_err)) +
    #scale_x_continuous(breaks = c(2011,2012,2013,2014,2015,2016,2017), labels = c(2011,2012,2013,2014,2015,2016,2017)) +
     scale_y_continuous(breaks = c(23,24,25,26,27,28,29), labels = c(23,24,25,26,27,28,29)) +
     theme_bw() + theme(panel.grid = element_blank(),
                       axis.text = element_text(size = 12),
                       axis.title = element_text(size = 14)) + xlab("") + ylab("Temperature")

# Deming Regression & plot
library(SimplyAgree)
CVroms = (sd(df$ROMS_temp, na.rm=T)/mean(df$ROMS_temp, na.rm=T))*100
CVstr = (sd(df$TEMP_C, na.rm=T)/mean(df$TEMP_C, na.rm=T))*100
err_ratio = ((CVroms*mean(df$ROMS_temp))^2)/((CVstr*mean(df$TEMP_C))^2)
d = df[,c('TEMP_C','ROMS_temp')]

dem = dem_reg(x = 'TEMP_C',
              y = 'ROMS_temp',
              data = d,
              error.ratio = err_ratio,
              weighted = FALSE)
dem
#save(dem, "deming.RData")
#plot(dem, pch = 19, color = alpha("#666666", 0.06))

ggplot(df, aes(x=TEMP_C, y=ROMS_temp)) + 
  geom_abline(slope = 1, intercept = 0, col = "black", lty = "dotted", lwd=0.9) +
  geom_point(color = alpha("#666666", 0.06)) +
  xlim(c(21.5,30.5)) + ylim(c(21.5,30.5)) +
  geom_abline(slope = dem$model$coef[2], intercept = dem$mod$coef[1], col = "red", lwd=0.6) +
  geom_density_2d(aes(color = ..level..)) +
  theme_bw() + theme(panel.grid = element_blank())

# Linear regression & plot
#mod = lm(df$ROMS_temp~df$TEMP_C)
#summary(mod)

# correlation plot
#png(paste0("Figures/MS_v3/Fig2_corrplot.png"), height = 5, width = 5, res = 500, units = "in")
#plot(df$TEMP_C, df$ROMS_temp, pch = 19, col = alpha("#666666", 0.06), xlim = c(21.5,30.5),
#     ylim = c(21.5,30.5), xlab = expression(paste("MHI STR Temp (",degree*C,")")),
#     ylab = expression(paste("MHI ROMS Temp (",degree*C,")")))
#abline(0,1, lwd=1.3, lty = "dashed")
#par(new = T)
#abline(a = dem$model$coef[1], b = dem$mod$coef[2], col = "red", lwd=1.3)
#dev.off()


# correlations by week for 3-hrly, daily data
r2 = df %>%
  mutate(week = cut.Date(as.Date(TIMESTAMP_3h), breaks = "1 month")) %>%
  group_by(week) %>%
  summarise(r2_3hr = summary(lm(ROMS_temp~TEMP_C))$r.squared)

df %>%
  mutate(week = cut.Date(as.Date(TIMESTAMP_3h), breaks = "1 week")) %>%
  subset(week %in% "2014-03-03") %>%
  ggplot() + geom_point(aes(TEMP_C, ROMS_temp)) + theme_bw() + ggtitle("2014-03-03")

##############################################
# Correlations between biweekly means and cv #
# (relevant to coral reef research)          #
##############################################

biweekly = df %>%
  mutate(day = day(TIMESTAMP_3h), biweek = cut.Date(as.Date(TIMESTAMP_3h), breaks = "2 weeks")) %>%
  group_by(biweek) %>%
  summarise(ROMS_biweek_mean = mean(ROMS_temp, na.rm=T),
            STR_biweek_mean = mean(TEMP_C, na.rm=T),
            ROMS_biweek_cv = (sd(ROMS_temp, na.rm=T)/mean(ROMS_temp, na.rm=T))*100,
            STR_biweek_cv = (sd(TEMP_C, na.rm=T)/mean(TEMP_C, na.rm=T))*100,
            STR_biweek_sd = sd(TEMP_C, na.rm=T)*100)

biweekly$biweek = as.Date(biweekly$biweek)
biweekly$month = month(biweekly$biweek)
plot(biweekly$biweek, biweekly$STR_biweek_cv, pch = 19)
biweekly_summer = subset(biweekly, month >6 & month <11)

# plot timeseries of biweekly CV
plot(biweekly_summer$biweek, biweekly_summer$ROMS_biweek_cv, pch = 19, col="#00CCCC", ylim = c(0.54,3.35))
par(new=T)
plot(biweekly_summer$biweek, biweekly_summer$STR_biweek_cv, pch = 19, ylim = c(0.54,3.35))


biweek_cv_daily <- df %>%
  mutate(day = as.Date(TIMESTAMP_3h)) %>%
  group_by(day) %>%
  summarise(ROMS = mean(ROMS_temp, na.rm=T), STR = mean(TEMP_C, na.rm=T)) %>%
  mutate(biweek = cut.Date(day, breaks = "2 weeks")) %>%
  group_by(biweek) %>%
  summarise(ROMS_biweek_cv = (sd(ROMS, na.rm=T)/mean(ROMS, na.rm=T))*100,
            STR_biweek_cv = (sd(STR, na.rm=T)/mean(STR, na.rm=T))*100) 

biweek_cv_daily$biweek = as.Date(biweek_cv_daily$biweek)
plot(biweek_cv_daily$biweek, biweek_cv_daily$STR_biweek_cv)

  
mod1 = lm(biweekly$ROMS_biweek_mean~biweekly$STR_biweek_mean)
mod2 = lm(biweekly$ROMS_biweek_cv~biweekly$STR_biweek_cv)

plot(biweekly$STR_biweek_mean, biweekly$ROMS_biweek_mean, pch = 19, col = "black",
     ylim = c(23, 29), xlim = c(23,29), 
     xlab = expression(paste("Biweekly Mean STR Temp (",degree*C,")")),
     ylab = expression(paste("Biweekly Mean ROMS Temp (",degree*C,")")))
abline(0,1, lty = "dashed", lwd=1.3, col = "darkgray")
par(new = T)
abline(a = mod1$coefficients[1], b = mod1$coefficients[2], col = "red", lwd=1.3)

plot(biweekly$STR_biweek_cv, biweekly$ROMS_biweek_cv, pch = 19, col = "black",
     ylim = c(0.52, 3.5), xlim = c(0.52,3.5),
     xlab = expression(paste("STR Biweekly CV")),
     ylab = expression(paste("ROMS Biweekly CV")))
abline(0,1, lty = "dashed", lwd=1.3, col = "darkgray")
par(new = T)
abline(a = mod2$coefficients[1], b = mod2$coefficients[2], col = "red", lwd=1.3)

summary(mod1)
summary(mod2)


############################################
# Overall skill metrics:                   #
# Correlation, bias, absolute error, RMSE  #
############################################

df_skill <- df %>% 
  summarise(corr = cor(TEMP_C, ROMS_temp),
            err = mean(ROMS_temp-TEMP_C),
            abs_err = mean(abs(ROMS_temp-TEMP_C)),
            rmse = sqrt(mean((ROMS_temp-TEMP_C)^2)))

# calcualte & plot RMSE  by site

df_site = df %>% group_by(OCC_SITEID) %>% summarise_at(vars(LONGITUDE, LATITUDE), funs(mean(., na.rm=T)))
site = df %>% group_by(OCC_SITEID) %>% summarise(corr = cor(TEMP_C, ROMS_temp),
                                                 rmse = sqrt(mean((ROMS_temp-TEMP_C)^2)),
                                                 err = mean(ROMS_temp-TEMP_C),
                                                 abs_err = mean(abs(ROMS_temp-TEMP_C)),
                                                 dep = mean(DEPTH))
site = merge(site, merge(df_site, unique(df[,c(3,5)]), by = "OCC_SITEID"),
                         by = "OCC_SITEID")
island = df %>% group_by(LOCATION) %>% summarise(corr = cor(TEMP_C, ROMS_temp))
month = df %>% group_by(month) %>% summarise(corr = cor(TEMP_C, ROMS_temp))
year = df %>% group_by(year) %>% summarise(corr = cor(TEMP_C, ROMS_temp))

# site
set.seed(13)
buffer = 0.08
site %>%
  group_by(LOCATION) %>% 
  summarise(x_min = min(LONGITUDE) - buffer,
            x_max = max(LONGITUDE) + buffer,
            y_min = min(LATITUDE) - buffer,
            y_max = max(LATITUDE) + buffer) %>%
  ggplot() + geom_polygon(data = hawaii, aes(x=long, y=lat, fill = group), color = "black") + 
  scale_fill_manual(values = rep("lightgray",14), guide = "none") +
  geom_jitter(data = site, aes(x=LONGITUDE, y=LATITUDE, colour = rmse), size= 3.4, width = 0.04, height = 0.04) +
  scale_colour_viridis(option = "D") +
  geom_rect(mapping = aes(
    xmin = x_min,
    xmax = x_max,
    ymin = y_min,
    ymax = y_max), alpha = 0.001) + theme_bw() + theme(panel.background = element_blank(),
                                                       panel.grid = element_blank(),
                                                       strip.background = element_rect(fill = "white"),
                                                       strip.text = element_text(size = 14),
                                                       axis.title = element_text(size = 12)) +
  #facet_wrap(~LOCATION, scales = "free") + 
  scale_y_continuous(labels = label_number(accuracy = 0.2)) +
  scale_x_continuous(labels = label_number(accuracy = 0.2)) +
  ylab("Latitude") + xlab("Longitude") + labs(fill = "RMSE")

set.seed(13)
site %>%
  group_by(LOCATION) %>% 
  summarise(x_min = min(LONGITUDE) - buffer,
            x_max = max(LONGITUDE) + buffer,
            y_min = min(LATITUDE) - buffer,
            y_max = max(LATITUDE) + buffer) %>%
  ggplot() + geom_polygon(data = hawaii, aes(x=long, y=lat, fill = group), color = "black") + 
  scale_fill_manual(values = rep("lightgray",14), guide = "none") +
  geom_jitter(data = site, aes(x=LONGITUDE, y=LATITUDE, colour = err), size= 3.4, width = 0.04, height = 0.04) +
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-0.16,0.378)) +
  geom_rect(mapping = aes(
    xmin = x_min,
    xmax = x_max,
    ymin = y_min,
    ymax = y_max), alpha = 0.001) + theme_bw() + theme(panel.background = element_blank(),
                                                       panel.grid = element_blank(),
                                                       strip.background = element_rect(fill = "white"),
                                                       strip.text = element_text(size = 14),
                                                       axis.title = element_text(size = 12)) +
  #facet_wrap(~LOCATION, scales = "free") + 
  scale_y_continuous(labels = label_number(accuracy = 0.2)) +
  scale_x_continuous(labels = label_number(accuracy = 0.2)) +
  ylab("Latitude") + xlab("Longitude") + labs(fill = "Error")

site %>%
  group_by(LOCATION) %>% 
  summarise(x_min = min(LONGITUDE) - buffer,
            x_max = max(LONGITUDE) + buffer,
            y_min = min(LATITUDE) - buffer,
            y_max = max(LATITUDE) + buffer) %>%
  ggplot() +
  annotation_map(map_data("world")) +
  geom_point(data = site, aes(LONGITUDE, LATITUDE, fill = abs_err), shape = 21, size=4) +
  scale_fill_distiller(palette="Reds", direction=1) +
  geom_rect(mapping = aes(
    xmin = x_min,
    xmax = x_max,
    ymin = y_min,
    ymax = y_max), alpha = 0.001) + theme_bw() + theme(panel.background = element_blank(),
                                                       panel.grid = element_blank(),
                                                       strip.background = element_rect(fill = "white"),
                                                       strip.text = element_text(size = 14),
                                                       axis.title = element_text(size = 12)) +
  facet_wrap(~LOCATION, scales = "free") + 
  scale_y_continuous(labels = label_number(accuracy = 0.2)) +
  scale_x_continuous(labels = label_number(accuracy = 0.2)) +
  ylab("Latitude") + xlab("Longitude") + labs(fill = "Abs. Error")

###################
# Taylor Diagrams #
###################

# depth
TaylorDiagram(df, 
              obs = "TEMP_C", 
              mod = "ROMS_temp", 
              group = "DEPTH", 
              normalise = TRUE,
              cols = "viridis",
              cex=2)

# rugosity; rugosity
TaylorDiagram(df, 
              obs = "TEMP_C", 
              mod = "ROMS_temp", 
              group = "rugosity", 
              normalise = TRUE,
              cols = "viridis",
              cex=2)

# island
TaylorDiagram(df, 
              obs = "TEMP_C", 
              mod = "ROMS_temp",
              group = "LOCATION", 
              normalise = TRUE,
              cols = "hue",
              cex =2)

# season
df$m = factor(month.abb[df$month], levels = month.abb)
df = df %>% 
  mutate(s = factor(ifelse(m == "Jan" | m == "Feb" | m == "Mar", "Jan-Mar", 
                         ifelse(m == "Apr" | m == "May" | m == "Jun", "Apr-Jun",
                                ifelse(m == "Jul" | m == "Aug" | m == "Sep", "Jul-Sep", "Oct-Dec"))),levels = c("Jan-Mar","Apr-Jun","Jul-Sep","Oct-Dec")))
d2 = subset(df, year != 2010)
TaylorDiagram(d2, 
              obs = "TEMP_C", 
              mod = "ROMS_temp",
              group = "s", 
              normalise = TRUE,
              cols = "brewer1",
              cex =2,
              key.title = "Season")

# year
d2$y = factor(d2$year)
TaylorDiagram(d2, 
              obs = "TEMP_C", 
              mod = "ROMS_temp",
              group = "y", 
              normalise = TRUE,
              cols = "increment",
              cex =2,
              key.title = "YEAR")

# distance from shore
TaylorDiagram(df, 
              obs = "TEMP_C", 
              mod = "ROMS_temp",
              group = "distance", 
              normalise = TRUE,
              cols = "viridis",
              cex =2,
              key.title = "Distance Offshore (m)")

#######################################################
# GAMs: modelling spatiotemporal variability of error #
#######################################################

df$abs_err = abs(df$ROMS_temp-df$TEMP_C)
df$err = (df$ROMS_temp-df$TEMP_C)
df$yr_day = lubridate::yday(df$TIMESTAMP_3h)
hist(df$err) # normal distribution

# Generalized additive model (GAM) -- used mgcv() package
d2 = subset(df, year != 2010)

mod = gam(err ~ s(DEPTH, k=3) + s(distance, k=3) + s(yr_day, bs="cc", k=4) +
            factor(year) + te(LONGITUDE, LATITUDE, k=4), data = d2, 
          family = gaussian(link = "identity"), na.action = 'na.fail')

sp_mod = gam(err ~ te(LONGITUDE, LATITUDE, k=3), data = df, family = gaussian(link = "identity"), na.action = 'na.fail')

summary(mod)
plot(mod, rug=T)
plot(mod, select = 1, rug=T, shade=TRUE, ylim = c(-0.33,0.48), ylab = "Effect on MHIA Bias", xlab = "Depth (m)")
plot(mod, select = 2, rug=T, shade=TRUE, ylim = c(-0.33,0.48), ylab = "Effect on MHIA Bias", xlab = "Distance (m)")
plot(mod, select = 3, rug=T, shade=TRUE, ylim = c(-0.33,0.48), ylab = "Effect on MHIA Bias", xlab = "Day of Year")
plot(mod, select = 4, rug=T, cex.axis=3, cex.labs=3, cex.main = 3)
termplot(mod, rug=T, col.term = 1, ylim = c(-0.33,0.48))
gam.check(mod)
library(pammtools)
gg_tensor(sp_mod)
vis.gam(mod, view = c("LONGITUDE","LATITUDE"))

# Look at dredge
library(MuMIn)
mod.dredge = dredge(mod)

# plot residuals to look for spatial autocorrelation
df$resid = mod$residuals
df_mean <- df %>% group_by(OCC_SITEID, LOCATION) %>% summarise_at(vars(LONGITUDE, LATITUDE, resid), list(~mean(., na.rm=T)))
  
ggplot() + geom_point(df_mean, mapping = aes(x=LONGITUDE, y=LATITUDE, color = resid), size=3) +
  scale_color_gradientn(colors = c(low="blue",mid="white",high="red"))

df_mean %>%
  group_by(LOCATION) %>% 
  ggplot() +
  annotation_map(map_data("world")) +
  geom_point(data = df_mean, aes(LONGITUDE, LATITUDE, fill = resid), shape = 21, alpha = 0.9, size=3) +
  scale_fill_gradientn(colors = c(low="blue",mid="white",high="red")) + theme_bw() +
  facet_wrap(~LOCATION, scales = "free") +
  ylab("Latitude") + xlab("Longitude") +
  theme(legend.position = "none")   

###########################################
# Evaluate Rainfall patterns--            #
# Did not include this in the final paper #
# as precipitation data is not correlated #
# with runoff/river outflow               #
###########################################

load("~/EDS_ROMS_TempProj/ROMS_subtemp_validation/ROMSval_Rproject/Env_data/rainfall/outputs/EDS_Timeseries_2022-11-10.Rdata")
df$err = df$ROMS_temp - df$TEMP_C
df$abs_err = abs(df$ROMS_temp - df$TEMP_C)
library(zoo)
library(dplyr)
df_month <- df %>% group_by(month, year) %>% summarize_all(list(~mean(., na.rm=T)))
df_month$month <- sprintf('%02d', df_month$month)
df_month$myr = as.Date(as.yearmon(paste0(df_month$month, "-", df_month$year), "%m-%Y"))
sd_rain = sd(df_month$mean_Precipitation_CHIRPS_Daily_DY01)
mean_rain = mean(df_month$mean_Precipitation_CHIRPS_Daily_DY01)
sd_err = sd(df_month$err)
mean_err = mean(df_month$err)
df_month$rain_zscore = (df_month$mean_Precipitation_CHIRPS_Daily_DY01 - mean_rain)/sd_rain
df_month$err_zscore = (df_month$err - mean_err)/sd_err
ggplot(data = df_month) +
  geom_line(aes(x=myr, y=rain_zscore)) + 
  geom_line(aes(x=myr, y=err_zscore), color = "red") +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) + xlab("") + ylab("Z Score")

# subset for shallowest sites
shallow_month <- df %>%
  subset(DEPTH < 6) %>%
  group_by(month, year) %>% 
  summarize_all(list(~mean(., na.rm=T)))

shallow_month$month <- sprintf('%02d', shallow_month$month)
shallow_month$myr = as.Date(as.yearmon(paste0(shallow_month$month, "-", shallow_month$year), "%m-%Y"))
sd_rain = sd((shallow_month$mean_Precipitation_CHIRPS_Daily_DY01))
mean_rain = mean(shallow_month$mean_Precipitation_CHIRPS_Daily_DY01)
sd_err = sd((shallow_month$err))
mean_err = mean(shallow_month$err)
shallow_month$rain_zscore = (shallow_month$mean_Precipitation_CHIRPS_Daily_DY01 - mean_rain)/sd_rain
shallow_month$err_zscore = (shallow_month$err - mean_err)/sd_err
ggplot(data = shallow_month) +
  geom_line(aes(x=myr, y=rain_zscore)) + 
  geom_line(aes(x=myr, y=err_zscore), color = "red") +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.line.y.right = element_line(color = "red"), 
                     axis.ticks.y.right = element_line(color = "red"),
                     axis.text.y.right = element_text(color = "red"), 
                     axis.title.y.right = element_text(color = "red"),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) + xlab("") + ylab("Z Score")

# look for correlation/lags
library(gtools)
plot(running(shallow_month$rain_zscore, shallow_month$err_zscore, fun=cor, width=3), type = "b", ylab = "Corr")
ccf(shallow_month$rain_zscore, shallow_month$err_zscore, main="")

#############################################
# Evaluate temperature anomalies,           #
# number of days above bleaching threshold, #
# power spectral density plots--            #
# Management implications                   #
#############################################

# Temperature anomalies

climatologies <- df %>%
  dplyr::mutate(month = lubridate::month(TIMESTAMP_3h)) %>%
  group_by(month) %>%
  summarise(ROMS_clim = mean(ROMS_temp, na.rm=T), STR_clim = mean(TEMP_C, na.rm=T))
anomalies <- df %>%
  mutate(yrmo = as.Date(as.yearmon(paste0(df$month, "-", df$year), "%m-%Y"))) %>%
  group_by(yrmo) %>%
  summarise(ROMS = mean(ROMS_temp, na.rm=T), STR = mean(TEMP_C, na.rm=T)) %>%
  mutate(month = month(yrmo))

anomalies = merge(anomalies, climatologies, by = "month")
anomalies$ROMS_anom = anomalies$ROMS - anomalies$ROMS_clim
anomalies$STR_anom = anomalies$STR - anomalies$STR_clim

# plot time series
ggplot(data = anomalies) +
  geom_line(aes(x=yrmo, y=STR_anom), lty = "dashed") + 
  geom_line(aes(x=yrmo, y=ROMS_anom)) +
  geom_hline(yintercept = 0, lwd=0.2) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) + xlab("") + ylab("Temperature Anomaly")


# Calculate & plot number of days above bleaching threshold (MMM + 1 degree C)
# This is done by calculating MMM sequentially each year
# Remove years that bleached!

MMM_2011 <- df %>%
  subset(year<= 2011) %>%
  mutate(month = month(as.Date(TIMESTAMP_3h)), year = year(as.Date(TIMESTAMP_3h))) %>%
  group_by(month, year) %>%
  summarize(roms = mean(ROMS_temp, na.rm=T), str = mean(TEMP_C, na.rm=T))
ROMS12_max = max(MMM_2011$roms, na.rm=T)
STR12_max = max(MMM_2011$str, na.rm=T)

MMM_2012 <- df %>%
  subset(year<= 2012) %>%
  mutate(month = month(as.Date(TIMESTAMP_3h)), year = year(as.Date(TIMESTAMP_3h))) %>%
  group_by(month, year) %>%
  summarize(roms = mean(ROMS_temp, na.rm=T), str = mean(TEMP_C, na.rm=T))
ROMS13_max = max(MMM_2012$roms, na.rm=T)
STR13_max = max(MMM_2012$str, na.rm=T)

MMM_2013 <- df %>%
  subset(year<= 2013) %>%
  mutate(month = month(as.Date(TIMESTAMP_3h)), year = year(as.Date(TIMESTAMP_3h))) %>%
  group_by(month, year) %>%
  summarize(roms = mean(ROMS_temp, na.rm=T), str = mean(TEMP_C, na.rm=T))
ROMS14_max = max(MMM_2013$roms, na.rm=T)
STR14_max = max(MMM_2013$str, na.rm=T)

MMM_2014 <- df %>%
  subset(year<= 2013) %>%
  mutate(month = month(as.Date(TIMESTAMP_3h)), year = year(as.Date(TIMESTAMP_3h))) %>%
  group_by(month, year) %>%
  summarize(roms = mean(ROMS_temp, na.rm=T), str = mean(TEMP_C, na.rm=T))
ROMS15_max = max(MMM_2014$roms, na.rm=T)
STR15_max = max(MMM_2014$str, na.rm=T)

MMM_2015 <- df %>%
  subset(year<= 2013) %>%
  mutate(month = month(as.Date(TIMESTAMP_3h)), year = year(as.Date(TIMESTAMP_3h))) %>%
  group_by(month, year) %>%
  summarize(roms = mean(ROMS_temp, na.rm=T), str = mean(TEMP_C, na.rm=T))
ROMS16_max = max(MMM_2015$roms, na.rm=T)
STR16_max = max(MMM_2015$str, na.rm=T)

MMM_2016 <- df %>%
  subset(year!= 2014 & year!=2015 & year!=2017) %>%
  mutate(month = month(as.Date(TIMESTAMP_3h)), year = year(as.Date(TIMESTAMP_3h))) %>%
  group_by(month, year) %>%
  summarize(roms = mean(ROMS_temp, na.rm=T), str = mean(TEMP_C, na.rm=T))
ROMS17_max = max(MMM_2016$roms, na.rm=T)
STR17_max = max(MMM_2016$str, na.rm=T)

bleaching1 <- df %>%
  group_by(date = as.Date(TIMESTAMP_3h, format = "%Y-%m-%d")) %>%
  summarize(ROMS_daily = mean(ROMS_temp, na.rm=T), STR_daily = mean(TEMP_C, na.rm=T)) %>%
  mutate(year = year(date), month = month(date)) %>%
  mutate(ROMS_above_thr = ifelse(year==2012,
      ifelse(ROMS_daily > (ROMS12_max + 1), "1", "0"),
      ifelse(year==2013, 
             ifelse(ROMS_daily > (ROMS13_max + 1), "1", "0"),
             ifelse(year==2014,
                    ifelse(ROMS_daily > (ROMS14_max + 1), "1", "0"),
                    ifelse(year==2015,
                           ifelse(ROMS_daily > (ROMS15_max + 1), "1", "0"),
                           ifelse(year==2016,
                                  ifelse(ROMS_daily > (ROMS16_max + 1), "1", "0"),
                                  ifelse(year==2017,
                                         ifelse(ROMS_daily > (ROMS17_max + 1), "1", "0"),NA)))))),
      STR_above_thr = ifelse(year==2012,
                              ifelse(STR_daily > (STR12_max + 1), "1", "0"),
                              ifelse(year==2013, 
                                     ifelse(STR_daily > (STR13_max + 1), "1", "0"),
                                     ifelse(year==2014,
                                            ifelse(STR_daily > (STR14_max + 1), "1", "0"),
                                            ifelse(year==2015,
                                                   ifelse(STR_daily > (STR15_max + 1), "1", "0"),
                                                   ifelse(year==2016,
                                                          ifelse(STR_daily > (STR16_max + 1), "1", "0"),
                                                          ifelse(year==2017,
                                                                 ifelse(STR_daily > (STR17_max + 1), "1", "0"),NA)))))),
      yrmo = as.Date(as.yearmon(paste0(month, "-", year), "%m-%Y")))

bleaching1 = merge(bleaching1, anomalies[,c(2,3,4)], by = "yrmo")

bleaching2 = bleaching1 %>%
  mutate(ROMS_diff = ifelse(year == 2012 & ROMS_above_thr == 1,
                            (ROMS-(ROMS12_max+1)), 
                            ifelse(year == 2013 & ROMS_above_thr == 1,
                                   (ROMS-(ROMS13_max+1)),
                                   ifelse(year == 2014 & ROMS_above_thr == 1,
                                          (ROMS-(ROMS14_max+1)),
                                          ifelse(year == 2015 & ROMS_above_thr == 1,
                                                 (ROMS-(ROMS15_max+1)),
                                                 ifelse(year == 2016 & ROMS_above_thr == 1,
                                                        (ROMS-(ROMS16_max+1)),
                                                        ifelse(year == 2017 & ROMS_above_thr == 1,
                                                               (ROMS-(ROMS17_max+1)), NA)))))),
         STR_diff = ifelse(year == 2012 & STR_above_thr == 1,
                            (STR-(STR12_max+1)), 
                            ifelse(year == 2013 & STR_above_thr == 1,
                                   (STR-(STR13_max+1)),
                                   ifelse(year == 2014 & STR_above_thr == 1,
                                          (STR-(STR14_max+1)),
                                          ifelse(year == 2015 & STR_above_thr == 1,
                                                 (STR-(STR15_max+1)),
                                                 ifelse(year == 2016 & STR_above_thr == 1,
                                                        (STR-(STR16_max+1)),
                                                        ifelse(year == 2017 & STR_above_thr == 1,
                                                               (STR-(STR17_max+1)), NA))))))) %>%
  group_by(yrmo) %>%
  summarize(ROMS_days_above = length(ROMS_above_thr[ROMS_above_thr == 1]),
            STR_days_above = length(STR_above_thr[STR_above_thr ==1]),
            ROMS_severity = ROMS_days_above*mean(ROMS_diff, na.rm=T),
            STR_severity = STR_days_above*mean(STR_diff, na.rm=T))

roms = cbind(bleaching2[,c(1,2,4)], "ROMS")
names(roms) <- c('yrmo', 'days','severity','type')
str = cbind(bleaching2[,c(1,3,5)], "STR")
names(str) <- c('yrmo', 'days', 'severity', 'type')
bleaching3 = rbind(roms[roms$yrmo >= "2012-01-01",],
                   str[str$yrmo >= "2012-01-01",])

bleaching3 = subset(bleaching3, yrmo >= "2014-01-01")
ggplot(data = bleaching3, aes(x=yrmo, y=days, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#00CCCC", "black")) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) + xlab("") + ylab("Days above bleaching threshold")

######################################################
# Evaluate patterns in power spectral density plots  #   
# information on estimating the spectral density:    #
# https://online.stat.psu.edu/stat510/lesson/12/12.1 #
######################################################

library(afc)
library(ggplot2)
library(stats)
library(astsa)

# average temperature by 3-hour

# Here, we plot periodograms, which displays information about the periodic components
# of our temperature time series.The periodograms provide information on relative
# strengths that explain the variation in the time series.
# The spectral density is a frequency domain that characterizes a stationary time series.
# We shall "smooth" the spectral density by using centered moving averages
# (called the Daniell Kernal method, aka a non parametric estimation of spectral density)
# average between two time periods of time t in either direction
# This command supplies co-efficients that are applied as weights in averaging the original data

### Plot Spectral Density for “Cycles Per Unit Time”
# Same as above (where frequency axis is cycles per sampling interval).
# However, now convert frequency axis to cycles per unit time.
# Do this by extracting frequency that R returns and dividing by the length of the sample interval.
# Also, multiple the spectral density by 2 so that the area under the periodogram actually equals
# the variance of the time series.
# Note: This is the preferred unit for the spectral density plot (so go with these periodograms)

k = kernel("daniell",4)

# Define sampling interval (3 hours)

si = 3

# calculate “cycles per unit time”

# all
STR.psd = mvspec(df$TEMP_C,k,log="yes",plot=TRUE) # PSD calculation
STR.spx = log(STR.psd$freq)/si    # extract frequency
STR.spy = 2*(log(STR.psd$spec))   # divide by length of sampling interval

ROMS.psd = mvspec(df$ROMS_temp,k,log="yes",plot=TRUE) # PSD calculation
ROMS.spx = log(ROMS.psd$freq)/si    # extract frequency
ROMS.spy = 2*(log(ROMS.psd$spec))   # divide by length of sampling interval


# Plot PSD for “cycles per unit time”

# plot
plot(STR.spy~STR.spx, type="l", xlim = c(-3.6, -0.25), ylim = c(-15, 16), ylab = "", xlab = "")
par(new=T)
plot(ROMS.spy~ROMS.spx, type="l", col="#00CCCC", xlim = c(-3.6, -0.25), ylim = c(-15, 16))

# plot small segments to identify locations of peaks on the x axis
x11()
plot(STR.spy~STR.spx, type="l", xlim = c(-1.5, -1.4), ylim = c(-15, 16), ylab = "", xlab = "")
par(new=T)
plot(ROMS.spy~ROMS.spx, type="l", col="#00CCCC", xlim = c(-1.5, -1.4), ylim = c(-15, 16))
locator(n=1, type ="l")
# [1] -1.11005 # 9 hr
# [1] -1.341742 # 11.5 hr
# [1] -1.425477 # 12.5 hr
# [1] -1.617014  # 15.1 hr
# [1] -1.656209  # 15.7 hr
# [1] -1.79  # 18 hr
# [1] -2.06  # 24 hr

# Calculate period
1/(exp(-1.425325)) # *3 gives you time period

# Therefore, vertical lines should be drawn at these points:
line_9hr = log(1/3.03451)
line_11.5hr = log(1/3.825702)
line_12.5hr = log(1/4.159842)
line_15.1hr = log(1/5.038024)
line_15.7hr = log(1/5.239411)
line_18hr = log(1/5.870853)
line_24hr = log(1/7.84597)

# Plot PSD for “cycles per unit time”
plot(STR.spy~STR.spx, type="l", xlim = c(-3.55, -0.2310491), ylim = c(-13.43, 18.9), ylab = "", xlab = "")
par(new=T)
plot(ROMS.spy~ROMS.spx, type="l", col="#00CCCC", xlim = c(-3.55, -0.2310491), ylim = c(-13.43, 18.9))
abline(v= line_9hr, lty="dotted", lwd=0.3) 
abline(v= line_11.5hr, lty="dotted", lwd=0.3) 
abline(v= line_12.5hr, lty="dotted", lwd=0.3) 
abline(v= line_15.1hr, lty="dotted", lwd=0.3) 
abline(v= line_15.7hr, lty="dotted", lwd=0.3) 
abline(v= line_18hr, lty="dotted", lwd=0.3) 
abline(v= line_24hr, lty="dotted", lwd=0.3) 


#######################################################################################
# Make map of ROMS 10, 20, and 30m temperature anomaly rel. to surface for Sept. 2015 #
#######################################################################################

#load("Env_data/roms_hiig_reanalysis_1.RData")
#load("Env_data/roms_hiig_reanalysis_10.RData")
#load("Env_data/roms_hiig_reanalysis_20.RData")
#load("Env_data/roms_hiig_reanalysis_30.RData")

#df<- df %>%
#  mutate(year = lubridate::year(time), 
#         month = lubridate::month(time)) %>%
#  group_by(month, year, latitude, longitude) %>%
#  summarise(temperature = mean(temp, na.rm=T))

#df <- df %>% subset(year >= 2010)
#climatologies <- df %>%
#  group_by(month, latitude, longitude) %>%
#  summarise(temp_clim = mean(temperature, na.rm=T))

#anomalies = merge(df, climatologies, by = c("month", "latitude","longitude"))
#anomalies$temp_anom = anomalies$temperature - anomalies$temp_clim
#sept2015_20m = subset(anomalies, year == 2015 & month == 9)

# load pre-saved data for Sept 2015 anomalies
roms_1m = readRDS("Env_data/ROMS_sept2015_1m.rds")
roms_10m = readRDS("Env_data/ROMS_sept2015_10m.rds")
roms_20m = readRDS("Env_data/ROMS_sept2015_20m.rds")
roms_30m = readRDS("Env_data/ROMS_sept2015_30m.rds")
# normalize
roms_10m$temp_anom_norm = roms_10m$temp_anom - roms_1m$temp_anom
roms_20m$temp_anom_norm = roms_20m$temp_anom - roms_1m$temp_anom
roms_30m$temp_anom_norm = roms_30m$temp_anom - roms_1m$temp_anom

# plot anomaly differences maps (diff. in temp. anomaly relative to surface anomaly)
ggplot() +
  geom_tile(roms_10m, mapping=aes(x=longitude,y=latitude,fill=temp_anom_norm)) +
  #scale_fill_viridis(option = "B", limits = c(-0.563,0.416), breaks = c(-0.4,-0.2,0.0,0.2,0.4)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-0.563,0.416)) +
  ggtitle("Sept 2015: 10m") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

# density plots
ggplot() + 
  geom_density(roms_1m, mapping=aes(x=temp_anom, colour = "1m")) +
  geom_density(roms_10m, mapping=aes(x=temp_anom, colour = "10m")) +
  geom_density(roms_20m, mapping=aes(x=temp_anom, colour = "20m")) +
  geom_density(roms_30m, mapping=aes(x=temp_anom, colour = "30m")) +
  scale_color_manual(name = "", 
                     values = c("1m" = "black", "10m" = "maroon", "20m" = "green", "30m" = "purple"),
                     guide = guide_legend(override.aes = list(linetype = 1, shape = 16))) +
  theme_bw() + theme(panel.grid = element_blank())
  

#saveRDS(sept2015_1m, "ROMS_sept2015_1m.rds")
#saveRDS(sept2015_10m, "ROMS_sept2015_10m.rds")
#saveRDS(sept2015_20m, "ROMS_sept2015_20m.rds")
#saveRDS(sept2015_30m, "ROMS_sept2015_30m.rds")

#######
# End #
#######
