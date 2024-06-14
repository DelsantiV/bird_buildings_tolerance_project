setwd("D:/Projet woizos")
library(sf)
library(tidyverse)
library(fs)
library(terra)

STOC_positions_vect <- vect(paste("STOC/Points_location_full.shp", sep= ""))
STOC_positions <- st_read(paste("STOC/Points_location_full.shp", sep= ""))

Temp_raster = project(rast("Temperature and precipitations/rangetemp_1km_eu.tif"), "EPSG:4326")
Precip_raster = project(rast("Temperature and precipitations/precip_1km_eu.tif"), "EPSG:4326")


Temp_values <- extract(Temp_raster, STOC_positions_vect)
Precip_values <- extract(Precip_raster, STOC_positions_vect)

Temp_values <- Temp_values %>% mutate(ID = unlist(lapply(Temp_values$ID, function(ID) STOC_positions[ID,]$ID_point))) %>% filter(!is.na(rangetemp_1km_eu)) 
Temp_values <- rename(Temp_values,c(ID_point = ID, Temp_range = rangetemp_1km_eu))

Precip_values <- Precip_values %>% mutate(ID = unlist(lapply(Precip_values$ID, function(ID) STOC_positions[ID,]$ID_point))) %>% filter(!is.na(precip_1km_eu)) 
Precip_values <- rename(Precip_values,c(ID_point = ID, Precip = precip_1km_eu))

STOC_temp_and_precip <- st_drop_geometry(STOC_positions) %>% dplyr::filter(ID_point %in% Precip_values$ID_point) %>% dplyr::select(c(ID_point)) %>%
  left_join(Temp_values, by = "ID_point") %>% left_join(Precip_values, by = "ID_point")

write.csv(STOC_temp_and_precip, "STOC/Temp_and_precip_table_.csv")
