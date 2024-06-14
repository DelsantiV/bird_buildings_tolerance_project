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
Sites_positions <- st_read("Martin/Positions.shp", crs = 4326)
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

find_vol_and_dist <- function(df_row, buffer_dist, bat_df){
  gc(verbose=F)
  Com <- df_row$INSEE_COM
  Com_front <- df_row$COM_FRONT
  bats_in_com <- bat_df %>% filter(INSEE_COM %in% unlist(append(strsplit(Com_front,","), Com)))
  bat_in_dist <- st_is_within_distance(df_row,bats_in_com,buffer_dist)
  bat_vol <- sum(bats_in_com[unlist(bat_in_dist),]$volume)
  bat_surf <- sum(bats_in_com[unlist(bat_in_dist),]$surface)
  all_dists <- st_distance(st_sfc(df_row$geometry, crs = 4326),bats_in_com)
  units(all_dists) <- NULL
  sorted_dist <- sort(all_dists - sqrt(bats_in_com$surface/pi), decreasing = F, index.return = T)
  dist_matrix <- sorted_dist$x
  dist_ind <- sorted_dist$ix
  list(c(dist_matrix[1],dist_matrix[5],dist_matrix[10],dist_matrix[50]), bat_vol)
}

for (year in c("2008","2013", "2018")){
  cat("\n\n\n", year,"\n\n")
  gc(verbose=F)
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1]-2)
  registerDoSNOW(cl)
  clusterExport(cl, list("find_vol_and_dist", "get_bat_df","Sites_positions", "data_files"))
  clusterEvalQ(cl, c(library('sf'), library('tidyverse'), library('fs')))
  length(Sites_positions$ID)
  
  start.time = Sys.time()
  pb <- txtProgressBar(min=1, max=length(Sites_positions$ID), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  bat_vol_and_dist <- foreach(i = 1:length(Sites_positions$ID), .options.snow = opts) %dopar% find_vol_and_dist(Sites_positions[i,], 500,get_bat_df(Sites_positions[i,]$INSEE_COM, Sites_positions[i,]$COM_FRONT, year = year))
  
  dist_list <- lapply(bat_vol_and_dist, function(el) el[[1]])
  vol_list <- lapply(bat_vol_and_dist, function(el) el[[2]])
  Sites_positions <- Sites_positions %>% 
    add_column(V_500 = unlist(vol_list))%>% 
    add_column(D_1 = unlist(unlist(lapply(dist_list, function(el) el[1])))) %>% 
    add_column(D_5 = unlist(unlist(lapply(dist_list, function(el) el[2])))) %>% 
    add_column(D_10 = unlist(unlist(lapply(dist_list, function(el) el[3])))) %>% 
    add_column(D_50 = unlist(unlist(lapply(dist_list, function(el) el[4]))))
  
  names(Sites_positions)[(length(names(Sites_positions))-4):length(names(Sites_positions))] <- c(paste("V_500_",year,sep=""),paste("D_1_",year,sep=""),paste("D_5_",year,sep=""),paste("D_10_",year,sep=""),paste("D_50_",year,sep=""))
  
  #stop cluster
  stopCluster(cl)
  
  st_write(Sites_positions,paste("Martin/Positions_and_bat_infos_",year,"_corr.shp",sep=""), append = F)
  total.time <- Sys.time() - start.time
  cat("\n\nTime : ",total.time,"\n")
}
