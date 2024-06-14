setwd("C:/Users/vincent/Documents/Stage Vincent")
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
year = "2008"
data_files <- dir_ls(paste("Donnees bati/DEP_",year, sep = ""))
Sites_positions <- readRDS("Points_all_infos.RDS")
names(Sites_positions)[2:3] <- c("INSEE_COM", "COM_FRONT")
Com_data <- st_read("Donnees geo/Communes/COM_FRONT.shp")
Dep_list <- unique(Com_data$INSEE_DEP)


get_bat_df <- function(Com, Com_front, use_polygon_data = F){
  if (use_polygon_data){
    data_files_list <- data_files_full[which(grepl(".SHP", data_files_full))]
    dep_dec = 1
    filename <- "Departement_"
    }
  else {
    data_files_list <- data_files[which(grepl(".shp", data_files))]
    dep_dec = 0
    filename <- "DEP"
    }
  com_list <- unlist(append(strsplit(Com_front,","), Com))
  dep_list <- unique(unlist(lapply(com_list,function(text) substring(text,1,2))))
  selected_files <- data_files_list[which(grepl(paste(filename,as.integer(dep_list[1])+dep_dec,sep=""),data_files_list))]
  bat_df <- st_zm(st_read(selected_files[1], quiet = T), drop = T, what = "ZM")
  bat_df$NATURE <- NULL
  for (file in selected_files[-1]){
    bat_df_sup <- st_zm(st_read(file, quiet = T), drop = T, what = "ZM")
    bat_df_sup$NATURE <- NULL
    bat_df <- rbind(bat_df,bat_df_sup)    
  }
  if (st_crs(bat_df) != st_crs(4326)){bat_df <- st_transform(bat_df, crs = 4326)}
  bat_df
}

find_vol_and_dist <- function(df_row, buffer_dist, bat_df){
  Com <- df_row$INSEE_COM
  Com_front <- df_row$COM_FRONT
  bats_in_com <- bat_df %>% dplyr::filter(INSEE_COM %in% unlist(append(strsplit(Com_front,","), Com)))
  bat_in_dist <- st_is_within_distance(df_row,bats_in_com,buffer_dist)
  bat_vol <- sum(bats_in_com[unlist(bat_in_dist),]$volume)
  sorted_dist <- sort(st_distance(st_sfc(df_row$geometry, crs = 4326),bats_in_com) - sqrt(bats_in_com$surface/pi), decreasing = F, index.return = T)
  dist_matrix <- sorted_dist$x
  dist_ind <- sorted_dist$ix
  list(c(dist_matrix[1],dist_matrix[5],dist_matrix[10],dist_matrix[50]), bat_vol)
}

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-5)
registerDoSNOW(cl)
clusterExport(cl, list("find_vol_and_dist", "get_bat_df","Sites_positions", "data_files"))
clusterEvalQ(cl, c(library('sf'), library('tidyverse')))
length(Sites_positions$ID)

start.time = Sys.time()
pb <- txtProgressBar(min=1, max=length(Sites_positions$ID), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
bat_vol_and_dist <- foreach(i = 1:length(Sites_positions$ID), .options.snow = opts) %dopar% find_vol_and_dist(Sites_positions[i,], 500,get_bat_df(Sites_positions[i,]$INSEE_COM, Sites_positions[i,]$COM_FRONT))

dist_list <- lapply(bat_vol_and_dist, function(el) el[[1]])
vol_list <- lapply(bat_vol_and_dist, function(el) el[[2]])
Sites_positions <- Sites_positions %>% 
  add_column(V_500_2008 = unlist(vol_list))%>% 
  add_column(D_1_2008 = unlist(unlist(lapply(dist_list, function(el) el[1])))) %>% 
  add_column(D_5_2008 = unlist(unlist(lapply(dist_list, function(el) el[2])))) %>% 
  add_column(D_10_2008 = unlist(unlist(lapply(dist_list, function(el) el[3])))) %>% 
  add_column(D_50_2008 = unlist(unlist(lapply(dist_list, function(el) el[4]))))

#stop cluster
stopCluster(cl)

st_write(Sites_positions,"Martin/Positions_and_bat_infos_2008.shp", append = F)
total.time <- Sys.time() - start.time
cat("\n ",total.time)*

buffer_dist = 500
for (i in 2:20){ 
  df_row <- Sites_positions[i,]
  Com <- df_row$INSEE_COM
  Com_front <- df_row$COM_FRONT
  com_list <- unlist(append(strsplit(Com_front,","), Com))
  dep_list <- unique(unlist(lapply(com_list,function(text) substring(text,1,2))))
  cat("Point :",df_row$ID,"\n\n")
  for (year in c("2008","2013","2018")){
    data_files <- dir_ls(paste("Donnees bati/DEP_",year, sep = ""))
    data_files_full <- dir_ls(paste("Test download/",year,"_Batiment_Only/",sep=""))
    bat_df <- get_bat_df(Sites_positions[i,]$INSEE_COM, Sites_positions[i,]$COM_FRONT)
    bat_df_full  <- get_bat_df(Com, Com_front, use_polygon_data = T)
    bats_in_com <- bat_df %>% filter(INSEE_COM %in% unlist(append(strsplit(Com_front,","), Com)))
    bats_poly_point <- st_intersection(bat_df_full, st_buffer(df_row$geometry, buffer_dist))
    bat_in_dist <- st_is_within_distance(df_row,bats_in_com,buffer_dist)
    st_write(bats_in_com[unlist(bat_in_dist),], paste("Martin/Test_Bat/Bat_at_point_",df_row$ID,"_",year,".shp",sep= ""), quiet = T, append = F)
    st_write(bats_poly_point, paste("Martin/Test_Bat/Bat_at_point_",df_row$ID,"_",year,"_poly.shp",sep= ""), quiet = T, append = F)
    dist_matrix = sort(st_distance(st_sfc(df_row$geometry, crs = 4326),bats_in_com), decreasing = F)
    bat_vol <- sum(bats_in_com[unlist(bat_in_dist),]$volume)
    cat("Année",year,":\n")
    cat("Nb batiments commune :",dim(bats_in_com)[1],"   Nb batiments 500m :",dim(bats_in_com[unlist(bat_in_dist),])[1],"\nVolume total :",bat_vol,"\nSurface totale :",sum(bats_in_com[unlist(bat_in_dist),]$surface), "\nVolume moyen :",mean(bats_in_com[unlist(bat_in_dist),]$volume),"\nSurface moyenne :",mean(bats_in_com[unlist(bat_in_dist),]$surface),"\nHauteur moyenne :",mean(bats_in_com[unlist(bat_in_dist),]$HAUTEUR), "\n\n")
  }
  cat("--------------------------------------------------------------------\n\n")
}
