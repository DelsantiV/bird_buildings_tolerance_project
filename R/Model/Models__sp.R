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

year = "2018"
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


# Modèles par espèce
nb_iter = 100
Sites_infos <- ungroup(Site_ab_sp_df_other) %>% dplyr::filter(code_sp %in% other_sp_names) %>%
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
attribution_vect <- c(rep(T,round(nb_points*0.8)),rep(F,round(nb_points*0.2)))
attribution_matrix <- matrix(ncol = nb_points, nrow = nb_iter)
for (i in 1:nb_iter){
  attribution_matrix[i,] <- sample(attribution_vect)
}

Ab_df_other <- ungroup(Site_ab_sp_df_other) %>% group_by(code_sp) %>% summarise(Abundance = n_distinct(ID_point)) %>% arrange(desc(Abundance))
Ab_df_aquatic <- ungroup(Site_ab_sp_df_aquatic) %>% group_by(code_sp) %>% summarise(Abundance = n_distinct(ID_point)) %>% arrange(desc(Abundance))
Ab_df_raptor <- ungroup(Site_ab_sp_df_raptor) %>% group_by(code_sp) %>% summarise(Abundance = n_distinct(ID_point)) %>% arrange(desc(Abundance))
write.csv(Ab_df_other, "STOC/Retained_sp.csv")
write.csv(Sites_infos, paste("STOC/Infos_at_site_other_species_",year,".csv",sep=""))
write.csv(attribution_matrix, "STOC/Attribution_train_preds.csv")

# Run Poisson models per species
start_sp <- 20
end_sp <- 50
count <- 0
for (sp in Ab_df$code_sp[start_sp:end_sp]){
  count = count +1
  cat(sp," (",count,"/",end_sp-start_sp+1,")\n\n")
  sp_assemblage = STOC_assemblage_other[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage) %>% mutate(Bt_v_250 = Bt_v_250 / 1e6, Bt_v_500 = Bt_v_500/1e6)
  Results_base <- list()
  Results_volume <- list()
  Results_dist <- list()
  Results_dist_nl <- list()
  Estimates_base <- list()
  Estimates_volume <- list()
  Estimates_dist <- list()
  Estimates_dist_nl <- list()
  MAE_matrix <- matrix(ncol = 4, nrow = nb_iter) 
  for (i in 1:nb_iter){
    attribution <- attribution_matrix[i,]
    Train_points <- Sites_infos_sp[attribution,]
    Test_points <- Sites_infos_sp[!attribution,]
    
    #prefit.glm <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250, 
    #                          data = Train_points[1:4000,])
    
    fit.glm_sp_base <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + Contains_Forest_250 + Contains_Urban_250 + Contains_Agri_250, 
                              data = Train_points)
    Estimates_base[[i]] <- summary(fit.glm_sp_base)[["coefficients"]]
    Preds_base <- predict(fit.glm_sp_base, newdata = Test_points, type = "response")
    MAE_base <- mean(abs(Preds_base - Test_points$Abundance))
    base_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_base, Real_Abundances = Test_points$Abundance)
    Results_base[[i]] <- base_model_perf
    
    fit.glm_sp_volume <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + Bt_v_250_quart + Contains_Forest_250 + Contains_Urban_250 + Contains_Agri_250, 
                                data = Train_points)
    Estimates_volume[[i]] <- summary(fit.glm_sp_volume)[["coefficients"]]
    Preds_volume <- predict(fit.glm_sp_volume, newdata = Test_points, type = "response")
    MAE_volume <- mean(abs(Preds_volume - Test_points$Abundance))
    volume_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_volume, Real_Abundances = Test_points$Abundance)
    Results_volume[[i]] <- volume_model_perf
    
    fit.glm_sp_dist <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + Dist_1_bat_log + Contains_Forest_250 + Contains_Urban_250 + Contains_Agri_250, 
                              data = Train_points)
    Estimates_dist[[i]] <- summary(fit.glm_sp_dist)[["coefficients"]]
    Preds_dist <- predict(fit.glm_sp_dist, newdata = Test_points, type = "response")
    MAE_dist <- mean(abs(Preds_dist - Test_points$Abundance))
    dist_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_dist, Real_Abundances = Test_points$Abundance)
    Results_dist[[i]] <- dist_model_perf
    
    fit.glm_sp_dist_nl <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + poly(Dist_1_bat_log,2) + Contains_Forest_250 + Contains_Urban_250 + Contains_Agri_250, 
                                 data = Train_points)
    Estimates_dist_nl[[i]] <- summary(fit.glm_sp_dist_nl)[["coefficients"]]
    Preds_dist_nl <- predict(fit.glm_sp_dist_nl, newdata = Test_points, type = "response")
    MAE_dist_nl <- mean(abs(Preds_dist_nl - Test_points$Abundance))
    dist_nl_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_dist_nl, Real_Abundances = Test_points$Abundance)
    Results_dist_nl[[i]] <- dist_nl_model_perf
    
    MAE_matrix[i,]  <- c(MAE_base, MAE_volume, MAE_dist, MAE_dist_nl)
    
  }
  MAE_df <- as.data.frame(MAE_matrix)
  names(MAE_df) <- c("MAE Base", "MAE Volume", "MAE Dist", "MAE Dist non linear")
  write.csv(MAE_df,paste("Results/NB_models_MAE_",sp,"_more_param.csv", sep = ""))
  saveRDS(Estimates_base, paste("Results/NB_models_Estimates_base_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Estimates_volume, paste("Results/NB_models_Estimates_volume_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Estimates_dist, paste("Results/NB_models_Estimates_dist_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Estimates_dist_nl, paste("Results/NB_models_Estimates_dist_nl_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Results_base, paste("Results/NB_models_Results_base_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Results_volume, paste("Results/NB_models_Results_volume_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Results_dist, paste("Results/NB_models_Results_dist_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Results_dist_nl, paste("Results/NB_models_Results_dist_nl_",sp,"_more_param.RDS", sep = ""))
}


# Evaluate model trained on a year period on another one
evaluate_model <- function(sp, train_year, test_year, sp_cat = other){
  if (!dir_exists(paste("Results/Preds_over_year/",sp_cat,"/",sp,"/", sep=""))){
    dir.create(paste("Results/Preds_over_year/",sp_cat,"/",sp,"/", sep=""))
  }
  sp_assemblage_train = read.csv(paste("STOC/Assemblage_",train_year,"_",sp_cat,".csv",sep = ""))[,sp]
  sp_assemblage_test = read.csv(paste("STOC/Assemblage_",test_year,"_",sp_cat,".csv",sep = ""))[,sp]
  Sites_infos_train <- read.csv(paste("STOC/Infos_at_site_",train_year,".csv",sep=""))
  Sites_infos_test <- read.csv(paste("STOC/Infos_at_site_",test_year,".csv",sep=""))
  Sites_infos_train <- Sites_infos_train %>% add_column(Abundance = sp_assemblage_train)
  Sites_infos_test <- Sites_infos_test %>% add_column(Abundance = sp_assemblage_test)
  
  fit.glm_sp_base <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + Contains_Forest_250 + Contains_Urban_250 + Contains_Agri_250 + Contains_Water_250, 
                               data = Sites_infos_train)
  Estimates_base <- summary(fit.glm_sp_base)[["coefficients"]]
  Preds_base <- predict(fit.glm_sp_base, newdata = Sites_infos_test, type = "response")
  base_model_perf <- data.frame(Site = Sites_infos_test$ID_point, Predictions = Preds_base, Real_Abundances = Sites_infos_test$Abundance)
  write.csv(Estimates_base, paste("Results/Preds_over_year/",sp_cat,"/",sp,"/",sp,"_Estimates_base.csv", sep = ""))
  write.csv(base_model_perf, paste("Results/Preds_over_year/",sp_cat,"/",sp,"/",sp,"_Results_base.csv", sep = ""))
  
  fit.glm_sp_dist_nl <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + poly(Dist_1_bat_log,2) + Contains_Forest_250 + Contains_Urban_250 + Contains_Agri_250 + Contains_Water_250, 
                               data = Sites_infos_train)
  Estimates_dist_nl <- summary(fit.glm_sp_dist_nl)[["coefficients"]]
  Preds_dist_nl <- predict(fit.glm_sp_dist_nl, newdata = Sites_infos_test, type = "response")
  dist_nl_model_perf <- data.frame(Site = Sites_infos_test$ID_point, Predictions = Preds_dist_nl, Real_Abundances = Sites_infos_test$Abundance)
  write.csv(Estimates_dist_nl, paste("Results/Preds_over_year/",sp_cat,"/",sp,"/",sp,"_Estimates_dist_nl.csv", sep = ""))
  write.csv(dist_nl_model_perf, paste("Results/Preds_over_year/",sp_cat,"/",sp,"/",sp,"_Results_dist_nl.csv", sep = ""))
}

count <- 0
ind_start <- 1
ind_end <- 10
for (sp in Ab_df_raptor$code_sp[ind_start:ind_end]){
  count <- count +1 
  cat(sp," (",count+ind_start-1,"/",ind_end,")\n\n")
  evaluate_model(sp, "2013", "2018", sp_cat= "raptor")
}
