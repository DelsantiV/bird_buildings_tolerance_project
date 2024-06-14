setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)
library(terra)

gc()

Dep_data = st_read("Donnees geo/Departements/DEPARTEMENT.shp")

year= "2008"
data_files <- dir_ls(paste("Test Download/",year,"_Batiment_Only", sep = ""))
st_batiment <- st_read(data_files[4], promote_to_multi = FALSE, type = 3)
st_batiment <- st_batiment[c('ID','HAUTEUR','geometry')]
st_batiment <- st_batiment %>% filter(as.numeric(HAUTEUR) > 3)
st_batiment$geometry <- st_zm(st_batiment$geometry, drop = T, what = "ZM")
st_batiment$surface <- st_area(st_batiment$geometry)
st_batiment$geometry <- st_centroid(st_batiment$geometry)
Dep <- Dep_data[which(st_intersects(st_batiment[3,], Dep_data, sparse = FALSE)),]$INSEE_DEP
st_batiment$INSEE_DEP <- rep(Dep)


for (file in data_files[-1]){
  if (grepl('.SHP',file)){
    data_batiment_dep <- st_read(file, promote_to_multi = FALSE, type = 3)
    data_batiment_dep <- data_batiment_dep[c('ID','HAUTEUR','geometry')]
    data_batiment_dep <- data_batiment_dep %>% filter(as.numeric(HAUTEUR) > 3)
    data_batiment_dep$geometry <- st_zm(data_batiment_dep$geometry, drop = T, what = "ZM")
    data_batiment_dep$surface <- st_area(data_batiment_dep$geometry)
    data_batiment_dep$geometry <- st_centroid(data_batiment_dep$geometry)
    test_index = 1
    Dep = character(0)
    while (length(Dep) ==0){
      Dep <- Dep_data[which(st_intersects(data_batiment_dep[test_index,], Dep_data, sparse = FALSE)),]$INSEE_DEP
      test_index <- test_index+1
    }
    data_batiment_dep$INSEE_DEP <- rep(Dep)
    st_batiment <- rbind(st_batiment, data_batiment_dep)
  }
}
#st_batiment$position = st_batiment["position"][[1]]
#st_batiment <- st_as_sf((st_batiment),remove = FALSE, wkt = "position", crs = st_crs(2154))
st_batiment_wgs84 <- st_transform(x = st_batiment, crs = 4326)
data_batiment_dep <- NULL
st_batiment_wgs84$volume <- st_batiment_wgs84$surface * st_batiment_wgs84$HAUTEUR

st_write(st_batiment_wgs84,paste("Donnees bati/Extracted_bat_centroids_",year,".shp", sep=""), append = FALSE)
saveRDS(st_batiment_wgs84,file=paste("Donnees bati/Extracted_bat_centroids_",year,".RDS", sep = ""))
#st_read("Données bâti/Extracted_bat_France.shp")