setwd("D:/Projet woizos")
library(raster)
library(sf)
library(tidyverse)
library(fs)
library(terra)

n_ini_sum = 0
n_sp_sum = 0
n_pos_sum = 0
n_NDVI_sum = 0
n_bat_sum = 0
n_final_sum = 0


sp_list <- readLines("STOC/sp_list_selected_V2.txt")
for (year in c("2008","2013","2018")){
  # Charger les données
  years_list <- as.character(seq(as.integer(year)-2, as.integer(year)+2))
  STOC_file = "STOC/data_FrenchBBS_point_France_100m_allSp_2001_2023fr.csv"
  STOC_positions = st_read(paste("STOC/Points_location_full_",year,".shp", sep= "")) # Vient de STOC_points_geom
  STOC_bat_info <- st_read(paste("STOC/Points_bat_infos",year,".shp", sep= ""), crs = 4326)  # Vient de Calculate_STOC_BAT_dist_parallel (tourné sur le cluster)
  names(STOC_bat_info)[(length(names(STOC_bat_info))-5):(length(names(STOC_bat_info))-1)] <- c("Dist_1_bat", "Bt_s_250", "Bt_v_250", "Bt_s_500", "Bt_v_500")
  STOC_bat_info <- STOC_bat_info %>% mutate(Dist_1_bat = if_else(Dist_1_bat < 1, 1, Dist_1_bat))
  NDVI_table <- read.csv(paste("STOC/NDVI_table_",year,".csv", sep = ""))[,-1]  # Vient de Import_NDVI
  CLC_table <- read.csv(paste("STOC/CLC_table_",year,".csv", sep = ""))   # Vient de Import_CLC
  Temp_and_precip_table <- read.csv("STOC/Temp_and_precip_table_.csv")    # Vient de Import_temp_and_prec
  
  
  STOC_data <- read.csv(STOC_file, sep = ';', dec = ',') %>% 
    dplyr::select(c(point, 
                    carre, 
                    annee, 
                    altitude, 
                    code_sp, 
                    nom_francais, 
                    abondance_filtre_tuckey))
  
  names(STOC_data) <- c("ID_point", 
                        "ID_carre", 
                        "annee", 
                        "altitude", 
                        "code_sp", 
                        "nom_sp", 
                        "Abundance")
  
  n_ini <-  dim(ungroup(STOC_data) %>% 
            filter(annee %in% years_list) %>% 
            group_by(ID_point) %>%
            summarise(ID_carre = unique(ID_carre), nb_sp = n_distinct(code_sp), altitude = unique(altitude)))[1]
  n_ini_sum <- n_ini_sum + n_ini
  STOC_data <- STOC_data  %>% 
    
    # Sélection des observations sur la bonne période d'années
    filter(annee %in% years_list) %>%
    
    # Sélection des espèces filtrées
    filter(code_sp %in% sp_list) %>%
    
    # Calcul des abondances moyennes sur la période d'années considérée
    group_by(across(-c(annee, Abundance))) %>% 
    summarise(Abundance = round(mean(Abundance)))
  
  
  n_sp <- n_ini - dim(ungroup(STOC_data) %>% 
            group_by(ID_point) %>%
            summarise(ID_carre = unique(ID_carre), nb_sp = n_distinct(code_sp), altitude = unique(altitude)))[1] 
  n_sp_sum <- n_sp_sum + n_sp
  STOC_data_point <- ungroup(STOC_data) %>% 
    group_by(ID_point) %>%
    summarise(ID_carre = unique(ID_carre), nb_sp = n_distinct(code_sp), altitude = unique(altitude)) %>%
    
    # Sélection des points filtrés par la position
    filter(ID_point %in% STOC_bat_info$ID_point) 
  
  n_pos <- n_ini - n_sp - dim(STOC_data_point)[1]
  n_pos_sum <- n_pos_sum + n_pos
  STOC_data_point <- STOC_data_point %>%
    
    # Ajout du NDVI
    filter(ID_point %in% NDVI_table$ID_point) %>%
    left_join(NDVI_table %>% 
                dplyr::select(c(ID_point,NDVI_mean_250, NDVI_var_250, NDVI_mean_500, NDVI_var_500)), 
              by = "ID_point") %>%
  
    # Ajout du CLC
    filter(ID_point %in% CLC_table$ID_point) %>%
    left_join(CLC_table %>% 
                dplyr::select(c(ID_point,CLC_250, CLC_500, Contains_Water_250, Contains_Forest_250, Contains_Urban_250, Contains_Agri_250)), 
              by = "ID_point") %>%
    
    # Ajout des données température range et précipitations
    filter(ID_point %in% Temp_and_precip_table$ID_point) %>%
    left_join(Temp_and_precip_table %>% 
                dplyr::select(c(ID_point,Temp_range, Precip)), 
              by = "ID_point") %>%
    
    # Ajout des données bâti
    filter(ID_point %in% STOC_bat_info$ID_point) %>%
    left_join(st_drop_geometry(STOC_bat_info) %>% 
                dplyr::select(c(ID_point, Dist_1_bat, Bt_s_250, Bt_v_250, Bt_s_500, Bt_v_500)), 
              by = "ID_point")
  
  n_NDVI <- n_ini - n_sp - n_pos - dim(STOC_data_point)[1]
  n_NDVI_sum <- n_NDVI_sum + n_NDVI
  STOC_data_point <- STOC_data_point %>% filter(Dist_1_bat < 2500)
  n_bat <- n_ini - n_sp -n_pos - n_NDVI - dim(STOC_data_point)[1]
  n_bat_sum <- n_bat_sum + n_bat
  n_final <- dim(STOC_data_point)[1]
  n_final_sum <- n_final_sum + n_final
  
  print(c(n_ini, n_sp, n_pos, n_NDVI, n_bat, n_final))
  
  STOC_data <- STOC_data %>% filter(ID_point %in% STOC_data_point$ID_point)
  
  
  STOC_data_ab <- STOC_data %>% group_by(ID_point,code_sp, nom_sp) %>% summarise(Abundance = round(mean(Abundance)))
  STOC_data_sp_ab <- STOC_data_ab %>% 
    group_by(across(c(code_sp,nom_sp))) %>%
    summarise(total_ab = sum(Abundance), Nb_sites = n(), max_ab = max(Abundance), mean_ab = mean(Abundance), std_ab = sd(Abundance)) %>%
    arrange(desc(total_ab))
  
  write.csv(STOC_data_sp_ab, file = paste("STOC/Stats_abondance_",year,".csv",sep=""))
  
  STOC_data_carre <- STOC_data %>% 
    group_by(ID_carre) %>%
    summarise(nb_sp = n_distinct(code_sp), all_sp = paste0(unique(nom_sp), collapse = ","), )
  
  STOC_data_point$nb_sp_carre <- unlist(lapply(STOC_data_point$ID_carre, function(carre){STOC_data_carre[which(STOC_data_carre$ID_carre == carre),]$nb_sp}))
  #STOC_data_point$all_sp_carre <- unlist(lapply(STOC_data_point$ID_carre, function(carre){STOC_data_carre[which(STOC_data_carre$ID_carre == carre),]$all_sp}))
  
  STOC_data_point <- STOC_data_point %>% 
    mutate(prop_sp = nb_sp/nb_sp_carre)
  
  write.csv(STOC_data_point, file = paste("STOC/STOC_points_data_ready_",year,".csv",sep=""))
  write.csv(STOC_data, file = paste("STOC/STOC_all_obs_",year,".csv",sep=""))
  
  STOC_data_point_uncounted <- STOC_data %>% uncount(Abundance)
  write.csv(STOC_data_point_uncounted, file = paste("STOC/STOC_points_data_uncounted_",year,".csv",sep=""))
}


cat("\n\n\n",c(n_ini_sum, n_sp_sum, n_pos_sum, n_NDVI_sum, n_bat_sum, n_final_sum))

# #STOC_positions06 <- STOC_positions[which(startsWith(STOC_positions$INSEE_COM,"06")),]
# STOC_bat_dist_06 <- read.csv(paste("STOC/STOC_bat_dist_06_",year,".csv", sep= ""))
# STOC_bat_dist_06 <- read.csv(paste("STOC/STOC_bat_vol_06_",year,".csv", sep= ""))
# STOC_bat_dist_06 <- STOC_bat_dist_06 %>% select(c(ID_point,Bat_dist))
# STOC_data_point_06 <- STOC_data_sp %>% filter(ID_point %in% STOC_bat_dist_06$ID_point)
# total_ab_06 <- STOC_data_point_06 %>%
#   group_by(code_sp, nom_francais) %>%
#   summarise(max_ab = max(abondance), mean_ab = mean(abondance), std_ab = sd(abondance), total_ab = sum(abondance)) %>%
#   arrange(desc(total_ab))
# high_ab_sp <- (total_ab_06 %>% filter(total_ab > 4))$code_sp
# STOC_data_point_06 <- STOC_data_point_06 %>%
#   filter(code_sp %in% high_ab_sp) %>%
#   left_join(STOC_bat_dist_06, by = "ID_point") %>%
#   left_join(STOC_bat_vol_06, by = "ID_point") %>%
#   arrange(code_sp) %>%
#   mutate(abondance = round(abondance)) %>%
#   mutate(Bat_dist = round(Bat_dist))
# STOC_data_point_06_ready <- STOC_data_point_06 %>% uncount(abondance)
# write.csv(STOC_data_point_06_ready, file = paste("STOC/STOC_points_data_with_dist_06_",year,".csv"))




annee_periods <- list()
annee_list = c(2008, 2013, 2018)
for (i in 1:length(annee_list)){
  annee = annee_list[i]
  annee_periods[[i]] <- as.character(seq(annee-2, annee+2))
}
#annee_periods[[length(annee_list) +1 ]] <- as.character(seq(as.integer(annee_periods[[length(annee_list)]][5])+1, 2023))
obs_in_period <- function(year_list, sp_row) {any(unlist(lapply(year_list, function(year) grepl(year, sp_row["annee"]))))}
obs_in_all <- function(sp_row) {all(unlist(lapply(annee_periods, function(year_list) obs_in_period(year_list, sp_row))))}

# Sélection des espèces observées sur toutes les périodes d'années
#STOC_positions <- STOC_positions[which(apply(STOC_positions,1,obs_in_all)),]
#STOC_data <- STOC_data[which(apply(STOC_data,1,obs_in_all)),]

# Sélection des espèces observées dans suffisamment de carrés différents
STOC_sp_carre <- STOC_data %>% dplyr::select(c(ID_carre,code_sp))
