setwd("~/Vincent")
library(sf)
library(tidyverse)
library(doParallel)
library(foreach)
library(pbapply)
library(fs)
library(doSNOW)
library(progress)

# Read data
gc()
Com_data <- st_read("Donnees geo/Communes/COM_FRONT.shp")
Dep_list <- unique(Com_data$INSEE_DEP)


get_bat_df <- function(Com, Com_front, year){
  data_files <- dir_ls(paste("Donnees bati/DEP_",year, sep = ""))
  com_list <- unlist(append(strsplit(Com_front,","), Com))
  dep_list <- unique(unlist(lapply(com_list,function(text) substring(text,1,2))))
  bat_df <- st_read(data_files[which(grepl(paste(dep_list[1],".shp",sep=""),data_files))], quiet = T)
  for (dep in dep_list[-1]){
    bat_df_sup <- st_read(data_files[which(grepl(paste(dep,".shp",sep=""),data_files))], quiet = T) 
    bat_df <- rbind(bat_df,bat_df_sup)
  }
  bat_df
}

calculate_bat_area_in_range <- function(buffer_dist, bat_row){
  area_bat <- bat_row$surface
  dist <- bat_row$dist_to_point
  bat_eq_radius <- sqrt(area_bat/pi)
  if (buffer_dist > dist + bat_eq_radius){
    final_area <- area_bat
    final_volume <- bat_row$volume
  }
  else {
    d1 <- (buffer_dist^2 - bat_eq_radius^2 + dist^2)/(2*dist)
    d2 <- dist-d1
    final_area <- (buffer_dist^2)*acos(d1/buffer_dist) - d1*sqrt(buffer_dist^2-d1^2) + (bat_eq_radius^2)*acos(d2/bat_eq_radius) - d2*sqrt(bat_eq_radius^2-d2^2)
    final_volume <- bat_row$volume * final_area/area_bat
  }
  c(final_area,final_volume)
}

find_vol_and_dist <- function(df_row, bat_df){
  gc(verbose=F)
  Com <- df_row$INSEE_COM
  Com_front <- df_row$COM_FRONT
  bats_in_com <- bat_df %>% filter(INSEE_COM %in% unlist(append(strsplit(Com_front,","), Com)))
  
  if (length(st_is_within_distance(df_row,bats_in_com,250)[1][[1]])){
    bat_in_500 <- st_is_within_distance(df_row,bats_in_com,500)
    bats_in_range <- bats_in_com[unlist(bat_in_500),]
  } else if (length(st_is_within_distance(df_row,bats_in_com,500)[1][[1]])){
    bat_in_1000 <- st_is_within_distance(df_row,bats_in_com,1000)
    bats_in_range <- bats_in_com[unlist(bat_in_1000),]
  } else{
    bats_in_range <- bats_in_com
  }
  all_dists <- st_distance(st_sfc(df_row$geometry, crs = 4326), bats_in_range)
  units(all_dists) <- NULL
  sorted_dist <- sort(all_dists - sqrt(bats_in_range$surface/pi), decreasing = F, index.return = T)
  dist_matrix <- sorted_dist$x
  dist_ind <- sorted_dist$ix
  bats_in_250 <- bats_in_com[dist_ind[dist_matrix < 250],] %>% add_column(dist_to_point = dist_matrix[dist_matrix < 250])
  bats_in_500 <- bats_in_com[dist_ind[dist_matrix < 500],] %>% add_column(dist_to_point = dist_matrix[dist_matrix < 500])
  if (dim(bats_in_250)[1]){
    bat_infos_250 <- apply(bats_in_250, 1, function(bat_row) calculate_bat_area_in_range(250, bat_row))
    bat_surf_250 <- sum(bat_infos_250[1,])
    bat_vol_250 <- sum(bat_infos_250[2,])
  } else {
    bat_surf_250 <- 0
    bat_vol_250 <- 0
  }
  
  if (dim(bats_in_500)[1]){
    bat_infos_500 <- apply(bats_in_500, 1, function(bat_row) calculate_bat_area_in_range(500, bat_row))
    bat_surf_500 <- sum(bat_infos_500[1,])
    bat_vol_500 <- sum(bat_infos_500[2,])
  } else{
    bat_surf_500 <- 0
    bat_vol_500 <- 0
  }  
  list(dist_matrix[1],bat_surf_250, bat_vol_250,bat_surf_500, bat_vol_500)
}


for (year in c("2008", "2013", "2018")){
  Sites_positions <- st_read(paste("STOC/Points_location_full_",year,".shp",sep=""), crs = 4326)
  cat("\n\n\n", year,"\n\n")
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1]-5)
  registerDoSNOW(cl)
  clusterExport(cl, list("find_vol_and_dist", "get_bat_df","Sites_positions", "data_files"))
  clusterEvalQ(cl, c(library('sf'), library('tidyverse'), library('fs')))
  length(Sites_positions$ID)
  
  start.time = Sys.time()
  pb <- txtProgressBar(min=1, max=length(Sites_positions$ID_point), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  bat_vol_and_dist <- foreach(i = 1:length(Sites_positions$ID_point), .options.snow = opts) %dopar% find_vol_and_dist(Sites_positions[i,], get_bat_df(Sites_positions[i,]$INSEE_COM, Sites_positions[i,]$COM_FRONT, year= year))

  Sites_positions <- Sites_positions %>% 
    add_column(D_1 = unlist(lapply(bat_vol_and_dist, function(el) el[[1]]))) %>% 
    add_column(S_250 = unlist(lapply(bat_vol_and_dist, function(el) el[[2]]))) %>% 
    add_column(V_250 = unlist(lapply(bat_vol_and_dist, function(el) el[[3]]))) %>% 
    add_column(S_500 = unlist(lapply(bat_vol_and_dist, function(el) el[[4]]))) %>% 
    add_column(V_500 = unlist(lapply(bat_vol_and_dist, function(el) el[[5]])))
  
  names(Sites_positions)[(length(names(Sites_positions))-4):length(names(Sites_positions))] <- c(paste("D_1_",year,sep=""),paste("S_250_",year,sep=""),paste("V_250_",year,sep=""),paste("S_500_",year,sep=""),paste("V_500_",year,sep=""))
  
  #stop cluster
  stopCluster(cl)
  
  st_write(Sites_positions,paste("STOC/Points_bat_infos",year,".shp",sep=""), append = F)
  total.time <- Sys.time() - start.time
  cat("\n\nTime : ",total.time,"\n")
}
