setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(terra)


write_simplified_CLC <- function(year){
  CLC_data <- st_read(paste("Donnees Land cover/",year,"/Data_CLC.shp",sep=""))
  CLC_data <- st_transform(CLC_data, crs = 4326)
  CLC_vect <- terra::vect(CLC_data)
  CLC_vect$CLC_type <- as.double(CLC_vect$CLC_type)
  empty_ras <- rast(CLC_vect, nrows = 21931*4, ncols = 24241*4)
  CLC_raster <- terra::rasterize(CLC_vect, empty_ras, "CLC_type")
  m <- c(0, 199, 1,
         200, 299, 2,
         300, 399, 3,
         400, 499, 4,
         500, 599, 5,
         600, 1000, 0)
  conversion_m <- matrix(m, ncol=3, byrow=TRUE)
  simplified_raster <- classify(CLC_raster, conversion_m, include.lowest=TRUE, brackets=TRUE)
  
  terra::writeRaster(simplified_raster, paste("Donnees Land cover/CLC_simplified_raster_",year,".tif",sep=""), overwrite=TRUE)
}

for (year in c("2006", "2012", "2018")){
  write_simplified_CLC(year)
}
