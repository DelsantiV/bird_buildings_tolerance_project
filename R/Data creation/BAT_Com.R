setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(purrr)

year = "2008"
st_batiment <- st_read(paste("Donnees bati/Extracted_bat_centroids_",year,"_COM.shp", sep = ""))
Com_data <- st_read("Donnees geo/Communes/COM_FRONT.shp")
#Com_data <- st_transform(x = Com_data, crs = 4326)
#Com_list <- apply(st_contains(Com_data, st_batiment, sparse = FALSE), 2, function(col) {(Com_data[which(col), ]$INSEE_COM)[1]})
#no_com = is.na(Com_list)
#st_batiment$INSEE_COM <- Com_list

Dep_list = unique(Com_data$INSEE_DEP)
for (Dep in Dep_list){
  st_batiment_dep <- st_batiment %>% filter(INSEE_DEP == Dep)
  print(paste(which(Dep_list == Dep)," /",length(Dep_list)))
  #Com_data_dep <- Com_data[which(Com_data$INSEE_DEP == Dep),]
  #st_batiment[st_batiment_dep_ind,]$INSEE_COM <- unlist(apply(st_contains(Com_data_dep, st_batiment[st_batiment_dep_ind,], sparse = FALSE), 2, function(col) {(Com_data_dep[which(col),]$INSEE_COM)[1]}))
  st_write(st_batiment_dep,paste("Donnees bati/DEP_",year,"/Extracted_bat_centroids_",year,"_DEP",Dep,".shp", sep=""), append = FALSE)
}

