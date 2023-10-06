setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)
library(terra)

CLC_data <- st_read("Donnees Land cover/2012/CLC12_FR_RGF_SHP/CLC12_FR_RGF.shp")
CLC_data <- st_transform(CLC_data, crs = 4326)
CLC_vect <- terra::vect(CLC_data)
CLC_vect$CODE_12 <- as.double(CLC_vect$CODE_12)
empty_ras <- rast(CLC_vect, nrows = 21931*4, ncols = 24241*4)
CLC_raster <- terra::rasterize(CLC_vect, empty_ras, "CODE_12")
m <- c(0, 199, 1,
       200, 299, 2,
       300, 399, 3,
       400, 499, 4,
       500, 599, 5,
       600, 1000, 0)
conversion_m <- matrix(m, ncol=3, byrow=TRUE)
simplified_raster <- classify(CLC_raster, conversion_m, include.lowest=TRUE, brackets=TRUE)

terra::writeRaster(simplified_raster, "Donnees Land cover/CLC_simplified_raster_2012.tif", overwrite=TRUE)

#project(simplified_raster,"EPSG:32631")
