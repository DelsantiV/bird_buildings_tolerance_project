setwd("C:/Users/vincent/Documents/Stage Vincent")
library(tidyverse)
library(fs)
library(ggplot2)
library(StatMatch)
library(ape)
library(PCAmixdata)
library(MASS)
library(gridExtra)
library(lme4)
library(DHARMa)


STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
code_sp_df <- STOC_sp_table %>% dplyr::select(code_sp, nom_fr)

year = "2013"
STOC_data_year <- read.csv(paste("STOC/STOC_points_data_ready_",year,".csv",sep=""))   # Vient de Observations_df_construction
STOC_data_year <- STOC_data_year %>% filter(Abundance > 0) %>% filter(code_sp != "COLLIV")
# Assemblages venant de Create_assemblage_df
STOC_assemblage_all <- read.csv(paste("STOC/Assemblage_",year,"_all.csv",sep = ""))
STOC_assemblage_aquatic <- read.csv(paste("STOC/Assemblage_",year,"_aquatic.csv",sep = ""))
aquatic_list <- names(STOC_assemblage_aquatic)[-1]
STOC_assemblage_raptor <- read.csv(paste("STOC/Assemblage_",year,"_raptor.csv",sep = ""))
raptor_list <- names(STOC_assemblage_raptor)[-1]
STOC_assemblage_other <- read.csv(paste("STOC/Assemblage_",year,"_other.csv",sep = ""))
other_sp_names <- names(STOC_assemblage_other)[-1]

sp_all <- unique(STOC_data_year$code_sp)
all_sp_names <- unlist(lapply(sp_all, function(code){code_sp_df[which(code_sp_df$code_sp == code),]$nom_fr}))


STOC_data_year <- STOC_data_year %>%
  dplyr::select(-c(X, NDVI_mean_500, NDVI_var_500, Bt_s_250, Bt_s_500)) %>%
  add_column(Dist_1_bat_log = log(STOC_data_year$Dist_1_bat))

# Moyenne des abondances sur les années d'observation
STOC_data_year <- ungroup(STOC_data_year) %>% 
  group_by(across(-c(annee, Abundance))) %>% 
  summarise(Abundance = round(mean(Abundance)))

Site_ab_sp_df_other <- ungroup(STOC_data_year) %>% dplyr::filter(code_sp %in% other_sp_names)
Site_ab_sp_df_aquatic <- ungroup(STOC_data_year) %>% dplyr::filter(code_sp %in% aquatic_list)
Site_ab_sp_df_raptor <- ungroup(STOC_data_year) %>% dplyr::filter(code_sp %in% raptor_list)


Sites_infos <- ungroup(STOC_data_year) %>% 
  dplyr::group_by(across(-c(code_sp, nom_sp, Abundance))) %>% 
  summarise()  
Sites_infos <- Sites_infos %>% 
  add_column(Bt_v_250_quart = Sites_infos$Bt_v_250, Bt_v_500_quart = Sites_infos$Bt_v_500)
Sites_infos$Bt_v_250_quart[Sites_infos$Bt_v_250_quart == 0] <- NA
Sites_infos$Bt_v_250_quart = ntile(Sites_infos$Bt_v_250_quart, 4)
Sites_infos$Bt_v_250_quart[is.na(Sites_infos$Bt_v_250_quart)] <- 0
Sites_infos$Bt_v_250_quart <- as.factor(Sites_infos$Bt_v_250_quart)
Sites_infos$Bt_v_500_quart[Sites_infos$Bt_v_500_quart == 0] <- NA
Sites_infos$Bt_v_500_quart = ntile(Sites_infos$Bt_v_500_quart, 4)
Sites_infos$Bt_v_500_quart[is.na(Sites_infos$Bt_v_500_quart)] <- 0
Sites_infos$Bt_v_500_quart <- as.factor(Sites_infos$Bt_v_500_quart)
Sites_infos$Contains_Water_250 <- as.factor(as.numeric(Sites_infos$Contains_Water_250))
Sites_infos$Contains_Forest_250 <- as.factor(as.numeric(Sites_infos$Contains_Forest_250))
Sites_infos$Contains_Urban_250 <- as.factor(as.numeric(Sites_infos$Contains_Urban_250))
Sites_infos$Contains_Agri_250 <- as.factor(as.numeric(Sites_infos$Contains_Agri_250))
Columns_to_scale <- c("NDVI_mean_250", "NDVI_var_250", "altitude","Dist_1_bat_log", "Temp_range", "Precip")
Sites_infos[,Columns_to_scale] <- scale(Sites_infos[,Columns_to_scale])
nb_points <- length(Sites_infos$ID_point)

Ab_df <- ungroup(Site_ab_sp_df_other) %>% group_by(code_sp) %>% summarise(Abundance = n_distinct(ID_point)) %>% arrange(desc(Abundance))
write.csv(Ab_df, "STOC/Retained_sp.csv")
write.csv(Sites_infos, paste("STOC/Infos_at_site_",year,".csv",sep=""))
