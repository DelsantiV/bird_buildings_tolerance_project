setwd("D:/Projet woizos")
library(tidyverse)
library(fs)
library(ggplot2)
library(MASS)
library(lme4)
library(splines)


STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
code_sp_df <- STOC_sp_table %>% dplyr::select(code_sp, nom_fr)
sp_list <- readLines("STOC/sp_list_selected.txt")

Sites_infos <- read.csv("STOC/Infos_at_site_all_years.csv")
Assemblage_df <- read.csv("STOC/Assemblage_all_years.csv")


# Modèles par espèce
nb_iter = 100
nb_points <- length(Sites_infos$ID_point)
attribution_vect <- c(rep(T,round(nb_points*0.8)),rep(F,round(nb_points*0.2)))
attribution_matrix <- matrix(ncol = nb_points, nrow = nb_iter)
for (i in 1:nb_iter){
  attribution_matrix[i,] <- sample(attribution_vect)
}


# Run Negative binomial models per species
calculate_model_sp_no_CLC <- function(sp, year = "All_years") {
  sp_assemblage <- Assemblage_df[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  Results_base <- list()
  Results_dist_agri <- list()
  Models_base <- list()
  Models_dist <- list()
  
  
  cat("\nSpecies modelling progress:\n")
  for (i in 1:nb_iter){
    setTxtProgressBar(pb,i)
    
    
    attribution <- attribution_matrix[i,]
    Train_points <- Sites_infos_sp[attribution,]
    Test_points <- Sites_infos_sp[!attribution,]
    
    fit.glm_sp_base <- glm.nb(Abundance~Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250, 
                              data = Train_points)
    Models_base[[i]] <- fit.glm_sp_base
    Preds_base <- predict(fit.glm_sp_base, newdata = Test_points, type = "response")
    base_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_base, Real_Abundances = Test_points$Abundance)
    Results_base[[i]] <- base_model_perf
    
    
  
    fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                   data = Train_points)
    Models_dist[[i]] <- fit.glm_sp_dist_agri
    Preds_dist_agri <- predict(fit.glm_sp_dist_agri, newdata = Test_points, type = "response")
    MAE_dist_agri <- mean(abs(Preds_dist_agri - Test_points$Abundance))
    dist_agri_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_dist_agri, Real_Abundances = Test_points$Abundance)
    Results_dist_agri[[i]] <- dist_agri_model_perf
    
    
  }
  cat("\n\n")
  close(pb)
  saveRDS(Results_base, paste("Models results/",year,"/NB_models_Results_base_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Results_dist_agri, paste("Models results/",year,"/NB_models_Results_dist_agri_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Models_base, paste("E:/Vincent/Woizos/Models results/",year,"/NB_models_base_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Models_dist, paste("E:/Vincent/Woizos/Models results/",year,"/NB_models_dist_agri_",sp,"_no_CLC.RDS", sep = ""))
}


main_pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                          max = length(sp_list), # Maximum value of the progress bar
                          style = 3,    # Progress bar style (also available style = 1 and style = 2)
                          width = 50,   # Progress bar width. Defaults to getOption("width")
                          char = "=")   # Character used to create the bar

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                   max = nb_iter, # Maximum value of the progress bar
                   style = 3,    # Progress bar style (also available style = 1 and style = 2)
                   width = 30,   # Progress bar width. Defaults to getOption("width")
                   char = "-")   # Character used to create the bar


for (j in 1:length(sp_list)){
  cat("\n\n",sp_list[j],":\n\nGeneral progress:\n")
  setTxtProgressBar(main_pb,j)
  cat("\n")
  calculate_model_sp_no_CLC(sp_list[j])
}
close(main_pb)