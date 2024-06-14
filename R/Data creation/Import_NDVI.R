setwd("D:/Projet woizos")
library(sf)
library(tidyverse)
library(fs)
library(terra)
library(ggplot2)

year = "2008"
STOC_year = "2008"
NDVI_raster = rast(paste("NDVI/NDVI_composite_",year,".tif", sep=""))
NDVI_raster$NDVI <- NDVI_raster[[1]]/10000
NDVI_raster$NDVI_1 <- NDVI_raster[[2]]/10000
NDVI_raster$NDVI_2 <- NDVI_raster[[3]]/10000
NDVI_raster$NDVI_variance <- NDVI_raster[[4]]/1e8
STOC_positions <- st_read(paste("STOC/Points_location_full_",STOC_year,".shp", sep= ""))
STOC_positions_vect <- vect(STOC_positions)

buffer_250 <- terra::buffer(STOC_positions_vect, 250)
buffer_500 <- terra::buffer(STOC_positions_vect, 500)

NDVI_values_250 <- extract(NDVI_raster, buffer_250, fun = mean)
NDVI_values_500 <- extract(NDVI_raster, buffer_500, fun = mean)


NDVI_values_250 <- NDVI_values_250 %>% mutate(ID = unlist(lapply(NDVI_values_250$ID, function(ID) STOC_positions[ID,]$ID_point))) %>% filter(!is.na(NDVI)) 
NDVI_values_250 <- rename(NDVI_values_250,c(ID_point = ID, NDVI_mean_250 = NDVI, NDVI_min_250 = NDVI_1, NDVI_max_250 = NDVI_2, NDVI_var_250 = NDVI_variance))

NDVI_values_500 <- NDVI_values_500 %>% mutate(ID = unlist(lapply(NDVI_values_500$ID, function(ID) STOC_positions[ID,]$ID_point))) %>% filter(!is.na(NDVI))
NDVI_values_500 <- rename(NDVI_values_500,c(ID_point = ID, NDVI_mean_500 = NDVI, NDVI_min_500 = NDVI_1, NDVI_max_500 = NDVI_2, NDVI_var_500 = NDVI_variance))

STOC_NDVI <- st_drop_geometry(STOC_positions) %>% dplyr::filter(ID_point %in% NDVI_values_250$ID_point) %>% dplyr::select(c(ID_point)) %>%
  left_join(NDVI_values_250, by = "ID_point") %>% left_join(NDVI_values_500, by = "ID_point")

write.csv(STOC_NDVI, paste("STOC/NDVI_table_",year,".csv", sep = ""))
