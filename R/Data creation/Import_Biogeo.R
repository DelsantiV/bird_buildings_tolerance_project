setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)
library(terra)
library(ggplot2)

year = "2018"
STOC_positions = st_read(paste("STOC/Points_location_full_",year,".shp", sep= ""))

Biogeo_data <- st_read("Biogeo/2016/Biogeo_FR_WGS84.shp")
biogeo_list <- apply(st_intersects(STOC_positions, Biogeo_data, sparse = F),1,
                     function(row){Biogeo_data[which(row),]$code})
biogeo_list[which(lengths(biogeo_list) == 0)] <- "Undefined"
STOC_Biogeo <- st_drop_geometry(STOC_positions) %>% dplyr::select(c(ID_point)) %>%
  add_column(biogeo_zone = unlist(biogeo_list))

write.csv(STOC_Biogeo, paste("STOC/Biogeo_table_",year,".csv", sep = ""))