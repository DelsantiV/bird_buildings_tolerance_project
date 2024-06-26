setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)
library(terra)
library(ggplot2)
library(raster)

NA_to_zero <- function(number){
  if (is.na(number)){to_return <- 0} else{to_return <- as.double(number)}
  to_return
}

year = "2018"
STOC_year = "2018"
STOC_positions <- st_read(paste("STOC/Points_location_full_",STOC_year,".shp", sep= ""))
STOC_positions_vect <- vect(STOC_positions)

buffer_250 <- terra::buffer(STOC_positions_vect, 250)
buffer_500 <- terra::buffer(STOC_positions_vect, 500)

simplified_raster <- raster(paste("Donnees Land cover/CLC_simplified_raster_",year,".tif", sep = "")) # Vient de Create_CLC_raster
plot(simplified_raster)

CLC_STOC_250_raw <- raster::extract(rast(simplified_raster), buffer_250)
names(CLC_STOC_250_raw) <- c("Num","CLC_250")

CLC_STOC_250 <- CLC_STOC_250_raw %>%
  mutate(CLC_250 = case_when(
    CLC_250 == 1 ~ "Urban",
    CLC_250 == 2 ~ "Agriculture",
    CLC_250 == 3 ~ "Forest and natural",
    CLC_250 == 4 ~ "Water",
    CLC_250 == 5 ~ "Water",
    TRUE ~ "No data")) %>% 
  group_by(Num) %>% 
  summarise(Contains_Urban_250 = NA_to_zero(table(CLC_250)["Urban"]) > 10,
            Contains_Agri_250 = NA_to_zero(table(CLC_250)["Agriculture"]) > 10,
            Contains_Forest_250 = NA_to_zero(table(CLC_250)["Forest and natural"]) > 10,
            Contains_Water_250 = NA_to_zero(table(CLC_250)["Water"]) > 10,
            CLC_250 = names(which.max(table(CLC_250))), ) %>%
  add_column(ID_point = STOC_positions_vect$ID_point)
CLC_STOC_250$Num <- NULL

CLC_STOC_500_raw <- raster::extract(rast(simplified_raster), buffer_500) 
names(CLC_STOC_500_raw) <- c("Num","CLC_500")
                             
CLC_STOC_500 <- CLC_STOC_500_raw %>%
  mutate(CLC_500 = case_when(
    CLC_500 == 1 ~ "Urban",
    CLC_500 == 2 ~ "Agriculture",
    CLC_500 == 3 ~ "Forest and natural",
    CLC_500 == 4 ~ "Water",
    CLC_500 == 5 ~ "Water",
    TRUE ~ "No data")) %>% 
  group_by(Num) %>% 
  summarise(Contains_Urban_500 = NA_to_zero(table(CLC_500)["Urban"]) > 10,
            Contains_Agri_500 = NA_to_zero(table(CLC_500)["Agriculture"]) > 10,
            Contains_Forest_500 = NA_to_zero(table(CLC_500)["Forest and natural"]) > 10,
            Contains_Water_500 = NA_to_zero(table(CLC_500)["Water"]) > 10,
            CLC_500 = names(which.max(table(CLC_500))), ) %>%
  add_column(ID_point = STOC_positions_vect$ID_point)
CLC_STOC_500$Num <- NULL

STOC_CLC <-  st_drop_geometry(STOC_positions) %>% dplyr::select(c(ID_point)) %>%
  left_join(CLC_STOC_250, by = "ID_point") %>% left_join(CLC_STOC_500, by = "ID_point") 

write.csv(STOC_CLC,paste("STOC/CLC_table_",STOC_year,".csv",sep=""))