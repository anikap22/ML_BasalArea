

##############
##
## ML predictions of BAI
##
##
##############

# Response variable BAI from FIA

# Covariates: MAT, MAP, N deposition, SMAP, NDVI (NASA), EVI (NASA), FPAR/LAI (MODIS),
# ET (MODIS), GPP/NPP (MODIS), soil thickness (ORNL), soil C (ORNL), global soil products (ORNL), 
# soil P (ORNL), soil type (ORNL)

# Other covariates that could be obtained from FIA?? (e.g. diversity)

require(tidyr)
require(dplyr)
require(ggplot2)
require(plyr) #for count function
require(MODISTools)
require(lubridate)
library(raster)
library(sp)
library(rgdal)
library(tmap)
library(ncdf4)

options(scipen = 6)

setwd('C:/Users/Anika/Documents/GradSchool/BAI_ML/')


### Import FIA with BAI calculated ----------------------

FIA <- readRDS("../FIA_CompetitionModel/output/FIA_plot_withgroups.RDS")
FIA <- FIA %>% dplyr::select(BA1, BA2, n, stdage, elev, soilC1, soilC2, understoryC1, understoryC2, litterC1,
                             litterC2, mt1, mt2, BA_change, LAT, LON, dep, sm, texture, fortypsm, STATECD,
                             t, MAT, MAP)
FIA$BAI <- FIA$BA_change/FIA$t

#write.csv(FIA[,c("LAT","LON")], file="fia_lat_lon.csv") #save lat lon for extracting other data

### Add N deposition ------------------------------
# already added but could do a better job (e.g. at several time points, average over number of years)


### Add SMAP from NASA ------------------------
# I think this is already added as sm
# could do a better job (e.g. average over number of years)


### Access MODIS --------------------------------
# access modis through modisTools: https://cran.r-project.org/web/packages/MODISTools/vignettes/modistools-vignette.html
# do approx. 10 year averages for MODIS data

products <- mt_products()
head(products)
View(products)
# ndvi/evi: MOD13Q1 (16 day)
# fpar/lai: MOD15A2H (8 day)
# et: MOD16A2 (8 day)
# gpp: MOD17A2H (8 day)
# npp: MOD17A3H (1 year)

bands <- mt_bands(product = "MOD13Q1")
head(bands)
#MOD13A3

dates <- mt_dates(product = "MOD13Q1", lat = 42, lon = -110)
head(dates)



### Add NDVI and EVI from MODIS ----------------------

FIA$NDVI <- NULL
FIA$EVI <- NULL

for(i in 1:nrow(FIA)){
  print(100*i/nrow(FIA)) #loop counter
  # find closest date to start and end measurements in FIA
  dates <- mt_dates(product = "MOD13Q1", lat = FIA$LAT[i], lon = FIA$LON[i])
  dates$calendar_date <- date(dates$calendar_date)
  startdate <- lubridate::date(date_decimal(FIA$mt1[i], tz = "UTC"))
  enddate <- lubridate::date(date_decimal(FIA$mt2[i], tz = "UTC"))
  fia_start <- dates$calendar_date[which(abs(dates$calendar_date - startdate) == min(abs(dates$calendar_date - startdate)))]
  #fia_end <- dates$calendar_date[which(abs(dates$calendar_date - enddate) == min(abs(dates$calendar_date - enddate)))]
  fia_end <- fia_start + years(1) #was taking way too long getting 10+ years of data
  if(length(fia_end) > 1){
    fia_end <- fia_end[1]
  }
  if(length(fia_start) > 1){
    fia_start <- fia_start[1]
  }
  
  test_ndvi <- mt_subset(product = "MOD13Q1",
                         lat = FIA$LAT[i],
                         lon = FIA$LON[i],
                         band = c("250m_16_days_NDVI","250m_16_days_EVI"),
                         start = fia_start,
                         end = fia_end,
                         internal = TRUE,
                         progress = TRUE)
  
  FIA$NDVI[i] <- mean(test_ndvi[test_ndvi$band == "250m_16_days_NDVI","value"], na.rm=T)
  FIA$EVI[i] <- mean(test_ndvi[test_ndvi$band == "250m_16_days_EVI","value"], na.rm=T)
}



### Add FPAR/LAI from MODIS -----------------------

FIA$FPAR <- NULL
FIA$LAI <- NULL

for(i in 1:nrow(FIA)){
  print(100*i/nrow(FIA)) #loop counter
  # find closest date to start and end measurements in FIA
  dates <- mt_dates(product = "MOD15A2H", lat = FIA$LAT[i], lon = FIA$LON[i])
  dates$calendar_date <- date(dates$calendar_date)
  startdate <- lubridate::date(date_decimal(FIA$mt1[i], tz = "UTC"))
  enddate <- lubridate::date(date_decimal(FIA$mt2[i], tz = "UTC"))
  fia_start <- dates$calendar_date[which(abs(dates$calendar_date - startdate) == min(abs(dates$calendar_date - startdate)))]
  #fia_end <- dates$calendar_date[which(abs(dates$calendar_date - enddate) == min(abs(dates$calendar_date - enddate)))]
  fia_end <- fia_start + years(2) #was taking way too long getting 10+ years of data
  if(length(fia_end) > 1){
    fia_end <- fia_end[1]
  }
  if(length(fia_start) > 1){
    fia_start <- fia_start[1]
  }
  
  fpar <- mt_subset(product = "MOD15A2H",
                         lat = FIA$LAT[i],
                         lon = FIA$LON[i],
                         band = c("Fpar_500m","Lai_500m"),
                         start = fia_start,
                         end = fia_end,
                         internal = TRUE,
                         progress = TRUE)
  
  FIA$FPAR[i] <- mean(fpar[fpar$band == "Fpar_500m","value"], na.rm=T)
  FIA$LAI[i] <- mean(fpar[fpar$band == "Lai_500m","value"], na.rm=T)
}



### Add ET from MODIS ------------------------

FIA$ET <- NULL

for(i in 1:nrow(FIA)){
  print(100*i/nrow(FIA)) #loop counter
  # find closest date to start and end measurements in FIA
  dates <- mt_dates(product = "MOD16A2", lat = FIA$LAT[i], lon = FIA$LON[i])
  dates$calendar_date <- date(dates$calendar_date)
  startdate <- lubridate::date(date_decimal(FIA$mt1[i], tz = "UTC"))
  enddate <- lubridate::date(date_decimal(FIA$mt2[i], tz = "UTC"))
  fia_start <- dates$calendar_date[which(abs(dates$calendar_date - startdate) == min(abs(dates$calendar_date - startdate)))]
  #fia_end <- dates$calendar_date[which(abs(dates$calendar_date - enddate) == min(abs(dates$calendar_date - enddate)))]
  fia_end <- fia_start + years(2) #was taking way too long getting 10+ years of data
  if(length(fia_end) > 1){
    fia_end <- fia_end[1]
  }
  if(length(fia_start) > 1){
    fia_start <- fia_start[1]
  }
  
  et <- mt_subset(product = "MOD16A2",
                    lat = FIA$LAT[i],
                    lon = FIA$LON[i],
                    band = "ET_500m",
                    start = fia_start,
                    end = fia_end,
                    internal = TRUE,
                    progress = TRUE)
  
  FIA$FPAR[i] <- mean(et[et$band == "ET_500m","value"], na.rm=T)
}


### Add GPP from MODIS ---------------------

FIA$GPP <- NULL

for(i in 1:nrow(FIA)){
  print(100*i/nrow(FIA)) #loop counter
  # find closest date to start and end measurements in FIA
  dates <- mt_dates(product = "MOD17A2H", lat = FIA$LAT[i], lon = FIA$LON[i])
  dates$calendar_date <- date(dates$calendar_date)
  startdate <- lubridate::date(date_decimal(FIA$mt1[i], tz = "UTC"))
  enddate <- lubridate::date(date_decimal(FIA$mt2[i], tz = "UTC"))
  fia_start <- dates$calendar_date[which(abs(dates$calendar_date - startdate) == min(abs(dates$calendar_date - startdate)))]
  #fia_end <- dates$calendar_date[which(abs(dates$calendar_date - enddate) == min(abs(dates$calendar_date - enddate)))]
  fia_end <- fia_start + years(2) #was taking way too long getting 10+ years of data
  if(length(fia_end) > 1){
    fia_end <- fia_end[1]
  }
  if(length(fia_start) > 1){
    fia_start <- fia_start[1]
  }
  
  gpp <- mt_subset(product = "MOD17A2H",
                  lat = FIA$LAT[i],
                  lon = FIA$LON[i],
                  band = "Gpp_500m",
                  start = fia_start,
                  end = fia_end,
                  internal = TRUE,
                  progress = TRUE)
  
  FIA$GPP[i] <- mean(gpp[gpp$band == "Gpp_500m","value"], na.rm=T)
}


### Add NPP from MODIS ---------------------

FIA$NPP <- NULL

for(i in 1:nrow(FIA)){
  print(100*i/nrow(FIA)) #loop counter
  # find closest date to start and end measurements in FIA
  dates <- mt_dates(product = "MOD17A3H", lat = FIA$LAT[i], lon = FIA$LON[i])
  dates$calendar_date <- date(dates$calendar_date)
  startdate <- lubridate::date(date_decimal(FIA$mt1[i], tz = "UTC"))
  enddate <- lubridate::date(date_decimal(FIA$mt2[i], tz = "UTC"))
  fia_start <- dates$calendar_date[which(abs(dates$calendar_date - startdate) == min(abs(dates$calendar_date - startdate)))]
  #fia_end <- dates$calendar_date[which(abs(dates$calendar_date - enddate) == min(abs(dates$calendar_date - enddate)))]
  fia_end <- fia_start + years(5) #was taking way too long getting 10+ years of data
  if(length(fia_end) > 1){
    fia_end <- fia_end[1]
  }
  if(length(fia_start) > 1){
    fia_start <- fia_start[1]
  }
  
  npp <- mt_subset(product = "MOD17A3H",
                   lat = FIA$LAT[i],
                   lon = FIA$LON[i],
                   band = "Npp_500m",
                   start = fia_start,
                   end = fia_end,
                   internal = TRUE,
                   progress = TRUE)
  
  FIA$NPP[i] <- mean(npp[npp$band == "Npp_500m","value"], na.rm=T)
}



### Add soil thickness from ORNL -----------------------

soilthick <- raster("data_raw/Global_Soil_Regolith_Sediment_1304/data/average_soil_and_sedimentary-deposit_thickness.tif")

#extract
FIAxy <- FIA[,c("LON","LAT")]
colnames(FIAxy) <- c("x","y") #lon=y, lat=x (reversed...?)
soilthickextract <- raster::extract(soilthick, FIAxy)

#visualize
tm_shape(soilthick) + tm_raster(palette = "Greys") + 
  tm_legend(outside = TRUE, text.size = .8) 

#put soil thickness into FIA df
FIA$soilthick <- soilthickextract



### Add soil C from ORNL (could bring in other c metrics from this file) ------------------------

soilC <- readOGR(dsn=path.expand("data_raw/WEST_SOIL_CARBON_1238/data/statsgo1poly"),
                            layer="statsgo1poly")

FIAxy <- FIA[,c("LON","LAT")]
pts <- SpatialPoints(FIAxy,
                     proj4string = CRS(proj4string(soilthick)))
pts_t <- spTransform(pts,
                     crs(soilC))

soilCextract <- over(pts_t, soilC) #want OC0_20cm which is organic C in top 20cm soil (I think)

#put soil C into FIA df
FIA$soilC_ORNL <- soilCextract$OC0_20cm

#visualize -- not working; read about tm_shape()
tm_shape(soilC) + tm_shape(palette = "Greys") +
  tm_legend(outside = TRUE, text.size = 0.8)



### Add soil P from ORNL (could bring in other p metrics from this file) -----------------------------

ncfname <- "data_raw/GLOBAL_PHOSPHORUS_DIST_MAP_1223/data/pforms_den.nc"
p_raster <- brick(ncfname, varname="tot") #getting total P (also have data for labile, occluded, secondary minearl P, apatite p)

#visualize
tm_shape(p_raster) + tm_raster(palette = "Greys") + 
  tm_legend(outside = TRUE, text.size = .8) 
plot(p_raster)

#extract (can only extract 1000 values at a time for some reason...)
FIAxy <- FIA[,c("LON","LAT")]
soilpextract <- NULL
k <- 0
while(k < nrow(FIAxy)){
  if(k+1000 < nrow(FIAxy)){
    temp <- raster::extract(p_raster, FIAxy[(1+k):(1000+k),])
  }else{
    temp <- raster::extract(p_raster, FIAxy[(1+k):nrow(FIAxy),])
  }
  soilpextract <- rbind(soilpextract, temp)
  k <- k + 1000
}

#put soil P into FIA df
FIA$soilP_ORNL <- as.vector(soilpextract)



### Add soil type from ORNL --------------------------

soiltype <- raster("data_raw/ZOBLERSOILDERIVED_540/data/z_soiltype.dat")
#extract
FIAxy <- FIA[,c("LON","LAT")]
colnames(FIAxy) <- c("x","y") #lon=y, lat=x (reversed...?)
soiltypeextract <- raster::extract(soiltype, FIAxy)

#put soil type into FIA df
FIA$soiltype <- soiltypeextract


### Export compiled data --------------------------




