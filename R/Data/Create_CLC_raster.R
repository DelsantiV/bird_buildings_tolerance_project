setwd("C:/Users/vincent/Documents/Stage Vincent")
library(raster)
library(sf)
library(tidyverse)
library(fs)
library(terra)

CLC_data <- st_read("Donnees Land cover/2012/CLC12_FR_RGF_SHP/CLC06R_FR_RGF.shp")
CLC_data <- st_transform(CLC_data, crs = 4326)
CLC_raster <- raster(crs = crs(CLC_data), vals = 0, resolution = c(0.0005, 0.0005), ext = extent(st_bbox(CLC_data))) %>%
  raster::rasterize(CLC_data, ., field = as.double(CLC_data$CODE_06))
conversion_m <- cbind(from = c(0,200,300,400,500, 900), to = c(199,299,399,499,599, 1000), becomes = c(1, 2, 3, 4, 5, 0))
simplified_raster <- reclassify(CLC_raster, conversion_m)
#levels(simplified_raster)[[1]]$desc <- c("No data", "Urbain", "Agriculture", "Natural", "Wet", "Water")
#plot(simplified_raster, col = c("white", "red", "yellow", "forestgreen", "cyan", "navy"))

writeRaster(simplified_raster, "Donnees Land cover/CLC_simplified_raster_2006.tif", overwrite=TRUE)


