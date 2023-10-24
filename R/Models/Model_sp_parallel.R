setwd("~/Vincent/")
library(tidyverse)
library(fs)
library(ggplot2)
library(ggeffects)
library(MASS)
library(lme4)
library(pbapply)
library(doParallel)
library(doSNOW)
library(progress)
library(splines)

year = "2018"
STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
code_sp_df <- STOC_sp_table %>% dplyr::select(code_sp, nom_fr)
Sites_infos <- read.csv(paste("STOC/Infos_at_site_",year,".csv",sep=""))
bat_dist_sd <- sd(Sites_infos$Dist_1_bat)
Sites_infos <- Sites_infos %>% mutate(Dist_1_bat = Sites_infos$Dist_1_bat/bat_dist_sd) 
Site_ab_sp_df_all <- read.csv(paste("STOC/Assemblage_",year,"_all.csv",sep=""))
Site_ab_sp_df_other <- read.csv(paste("STOC/Assemblage_",year,"_other.csv",sep=""))
Site_ab_sp_df_aquatic <- read.csv(paste("STOC/Assemblage_",year,"_aquatic.csv",sep=""))
Site_ab_sp_df_raptor <- read.csv(paste("STOC/Assemblage_",year,"_raptor.csv",sep=""))

# Modèles par espèce
nb_iter = 100
nb_points <- length(Sites_infos$ID_point)
attribution_vect <- c(rep(T,round(nb_points*0.8)),rep(F,round(nb_points*0.2)))
attribution_matrix <- matrix(ncol = nb_points, nrow = nb_iter)
for (i in 1:nb_iter){
  attribution_matrix[i,] <- sample(attribution_vect)
}

# Run Negative binomial models per species
calculate_model_sp <- function(sp, sp_cat = "all") {
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  Results_base <- list()
  Results_volume <- list()
  Results_dist <- list()
  Results_dist_nl <- list()
  Estimates_base <- list()
  Estimates_volume <- list()
  Estimates_dist <- list()
  Estimates_dist_nl <- list()
  Confint_base <- list()
  Confint_volume <- list()
  Confint_dist <- list()
  Confint_dist_nl <- list()
  MAE_matrix <- matrix(ncol = 4, nrow = nb_iter) 
  for (i in 1:nb_iter){
    attribution <- attribution_matrix[i,]
    Train_points <- Sites_infos_sp[attribution,]
    Test_points <- Sites_infos_sp[!attribution,]
    
    #prefit.glm <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250, 
    #                          data = Train_points[1:4000,])
    
    fit.glm_sp_base <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250, 
                              data = Train_points)
    estims <- summary(fit.glm_sp_base)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_base[[i]] <- estims
    Preds_base <- predict(fit.glm_sp_base, newdata = Test_points, type = "response")
    MAE_base <- mean(abs(Preds_base - Test_points$Abundance))
    base_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_base, Real_Abundances = Test_points$Abundance)
    Results_base[[i]] <- base_model_perf
    Confint_base[[i]] <- confint(fit.glm_sp_base, predictors[which(estims[,"Std. Error"] < 10)], level = 0.95)
    
    fit.glm_sp_volume <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + Bt_v_250_quart, 
                                data = Train_points)
    estims <- summary(fit.glm_sp_volume)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_volume[[i]] <- estims
    Preds_volume <- predict(fit.glm_sp_volume, newdata = Test_points, type = "response")
    MAE_volume <- mean(abs(Preds_volume - Test_points$Abundance))
    volume_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_volume, Real_Abundances = Test_points$Abundance)
    Results_volume[[i]] <- volume_model_perf
    Confint_volume[[i]] <- confint(fit.glm_sp_volume, predictors[which(estims[,"Std. Error"] < 10)], level = 0.95)
    
    fit.glm_sp_dist <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + Dist_1_bat_log, 
                              data = Train_points)
    estims <- summary(fit.glm_sp_dist)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_dist[[i]] <- estims
    Preds_dist <- predict(fit.glm_sp_dist, newdata = Test_points, type = "response")
    MAE_dist <- mean(abs(Preds_dist - Test_points$Abundance))
    dist_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_dist, Real_Abundances = Test_points$Abundance)
    Results_dist[[i]] <- dist_model_perf
    Confint_dist[[i]] <- confint(fit.glm_sp_dist, predictors[which(estims[,"Std. Error"] < 10)], level = 0.95)
    
    fit.glm_sp_dist_nl <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + poly(Dist_1_bat_log,2), 
                                 data = Train_points)
    estims <- summary(fit.glm_sp_dist_nl)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_dist_nl[[i]] <- estims
    Preds_dist_nl <- predict(fit.glm_sp_dist_nl, newdata = Test_points, type = "response")
    MAE_dist_nl <- mean(abs(Preds_dist_nl - Test_points$Abundance))
    dist_nl_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_dist_nl, Real_Abundances = Test_points$Abundance)
    Results_dist_nl[[i]] <- dist_nl_model_perf
    Confint_dist_nl[[i]] <- confint(fit.glm_sp_dist_nl, predictors[which(estims[,"Std. Error"] < 10)], level = 0.95)
    
    
    MAE_matrix[i,]  <- c(MAE_base, MAE_volume, MAE_dist, MAE_dist_nl)
    
  }
  MAE_df <- as.data.frame(MAE_matrix)
  names(MAE_df) <- c("MAE Base", "MAE Volume", "MAE Dist", "MAE Dist non linear")
  write.csv(MAE_df,paste("Results/Results_by_sp/",sp_cat,"/NB_models_MAE_",sp,".csv", sep = ""))
  saveRDS(Estimates_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_base_",sp,".RDS", sep = ""))
  saveRDS(Estimates_volume, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_volume_",sp,".RDS", sep = ""))
  saveRDS(Estimates_dist, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_dist_",sp,".RDS", sep = ""))
  saveRDS(Estimates_dist_nl, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_dist_nl_",sp,".RDS", sep = ""))
  saveRDS(Results_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_base_",sp,".RDS", sep = ""))
  saveRDS(Results_volume, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_volume_",sp,".RDS", sep = ""))
  saveRDS(Results_dist, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_dist_",sp,".RDS", sep = ""))
  saveRDS(Results_dist_nl, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_dist_nl_",sp,".RDS", sep = ""))
  saveRDS(Confint_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_base_",sp,".RDS", sep = ""))
  saveRDS(Confint_volume, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_volume_",sp,".RDS", sep = ""))
  saveRDS(Confint_dist, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_dist_",sp,".RDS", sep = ""))
  saveRDS(Confint_dist_nl, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_dist_nl_",sp,".RDS", sep = ""))
}


# Run Negative binomial models per species
calculate_model_sp_more_param <- function(sp, sp_cat = "all") {
  
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  Results_base <- list()
  Results_volume <- list()
  Results_dist <- list()
  Results_dist_nl <- list()
  Estimates_base <- list()
  Estimates_volume <- list()
  Estimates_dist <- list()
  Estimates_dist_nl <- list()
  Confint_base <- list()
  Confint_volume <- list()
  Confint_dist <- list()
  Confint_dist_nl <- list()
  MAE_matrix <- matrix(ncol = 4, nrow = nb_iter) 
  for (i in 1:nb_iter){
    attribution <- attribution_matrix[i,]
    Train_points <- Sites_infos_sp[attribution,]
    Test_points <- Sites_infos_sp[!attribution,]
    
    #prefit.glm <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250, 
    #                          data = Train_points[1:4000,])
    
    fit.glm_sp_base <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Forest_250 + Contains_Water_250 + Contains_Urban_250 + Contains_Agri_250, 
                              data = Train_points)
    estims <- summary(fit.glm_sp_base)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_base[[i]] <- estims
    Preds_base <- predict(fit.glm_sp_base, newdata = Test_points, type = "response")
    MAE_base <- mean(abs(Preds_base - Test_points$Abundance))
    base_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_base, Real_Abundances = Test_points$Abundance)
    Results_base[[i]] <- base_model_perf
    Confint_base[[i]] <- confint(fit.glm_sp_base, predictors[which(estims[,"Std. Error"] < 10)], level = 0.95)
    
    fit.glm_sp_volume <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Bt_v_250_quart + Contains_Forest_250 + Contains_Water_250 + Contains_Urban_250 + Contains_Agri_250, 
                                data = Train_points)
    estims <- summary(fit.glm_sp_volume)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_volume[[i]] <- estims
    Preds_volume <- predict(fit.glm_sp_volume, newdata = Test_points, type = "response")
    MAE_volume <- mean(abs(Preds_volume - Test_points$Abundance))
    volume_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_volume, Real_Abundances = Test_points$Abundance)
    Results_volume[[i]] <- volume_model_perf
    Confint_volume[[i]] <- confint(fit.glm_sp_volume, predictors[which(estims[,"Std. Error"] < 10)], level = 0.95)
    
    fit.glm_sp_dist <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Dist_1_bat_log + Contains_Forest_250 + Contains_Water_250 + Contains_Urban_250 + Contains_Agri_250, 
                              data = Train_points)
    estims <- summary(fit.glm_sp_dist)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_dist[[i]] <- estims
    Preds_dist <- predict(fit.glm_sp_dist, newdata = Test_points, type = "response")
    MAE_dist <- mean(abs(Preds_dist - Test_points$Abundance))
    dist_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_dist, Real_Abundances = Test_points$Abundance)
    Results_dist[[i]] <- dist_model_perf
    Confint_dist[[i]] <- confint(fit.glm_sp_dist, predictors[which(estims[,"Std. Error"] < 10)], level = 0.95)
    
    fit.glm_sp_dist_nl <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + poly(Dist_1_bat_log,2, raw = T) + Contains_Forest_250 + Contains_Water_250 + Contains_Urban_250 + Contains_Agri_250, 
                                 data = Train_points)
    estims <- summary(fit.glm_sp_dist_nl)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_dist_nl[[i]] <- estims
    Preds_dist_nl <- predict(fit.glm_sp_dist_nl, newdata = Test_points, type = "response")
    MAE_dist_nl <- mean(abs(Preds_dist_nl - Test_points$Abundance))
    dist_nl_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_dist_nl, Real_Abundances = Test_points$Abundance)
    Results_dist_nl[[i]] <- dist_nl_model_perf
    Confint_dist_nl[[i]] <- confint(fit.glm_sp_dist_nl, predictors[which(estims[,"Std. Error"] < 10)], level = 0.95)
    
    MAE_matrix[i,]  <- c(MAE_base, MAE_volume, MAE_dist, MAE_dist_nl)
    
  }
  MAE_df <- as.data.frame(MAE_matrix)
  names(MAE_df) <- c("MAE Base", "MAE Volume", "MAE Dist", "MAE Dist non linear")
  write.csv(MAE_df,paste("Results/Results_by_sp/",sp_cat,"/NB_models_MAE_",sp,"_more_param.csv", sep = ""))
  saveRDS(Estimates_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_base_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Estimates_volume, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_volume_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Estimates_dist, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_dist_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Estimates_dist_nl, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_dist_nl_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Results_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_base_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Results_volume, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_volume_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Results_dist, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_dist_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Results_dist_nl, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_dist_nl_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Confint_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_base_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Confint_volume, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_volume_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Confint_dist, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_dist_",sp,"_more_param.RDS", sep = ""))
  saveRDS(Confint_dist_nl, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_dist_nl_",sp,"_more_param.RDS", sep = ""))
}


# Run Negative binomial models per species
calculate_model_sp_no_CLC <- function(sp, sp_cat = "all") {
  
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  Results_base <- list()
  Results_CLC <- list()
  Results_dist_agri <- list()
  Estimates_base <- list()
  Estimates_CLC <- list()
  Estimates_dist_agri <- list()
  Confint_base <- list()
  Confint_CLC <- list()
  Confint_dist_agri <- list()
  MAE_matrix <- matrix(ncol = 3, nrow = nb_iter) 
  for (i in 1:nb_iter){
    attribution <- attribution_matrix[i,]
    Train_points <- Sites_infos_sp[attribution,]
    Test_points <- Sites_infos_sp[!attribution,]
    
    fit.glm_sp_base <- glm.nb(Abundance~Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250, 
                              data = Train_points)
    estims <- summary(fit.glm_sp_base)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_base[[i]] <- estims
    Preds_base <- predict(fit.glm_sp_base, newdata = Test_points, type = "response")
    MAE_base <- mean(abs(Preds_base - Test_points$Abundance))
    base_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_base, Real_Abundances = Test_points$Abundance)
    Results_base[[i]] <- base_model_perf
    confi_int <- confint(fit.glm_sp_base, predictors[which(estims[,"Pr(>|z|)"] < 0.05)], level = 0.95)
    for (predic in predictors){
      if (!(predic %in% rownames(confi_int))){
        confi_int <- rbind(confi_int,c(NA, NA))
        rownames(confi_int)[length(rownames(confi_int))] <- predic
      }
    }
    Confint_base[[i]] <- confi_int
    
    
    fit.glm_sp_CLC <- glm.nb(Abundance~Temp_range + CLC_250 + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Forest_250 + Contains_Water_250 + Contains_Urban_250 + Contains_Agri_250, 
                              data = Train_points)
    estims <- summary(fit.glm_sp_CLC)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_CLC[[i]] <- estims
    Preds_CLC <- predict(fit.glm_sp_CLC, newdata = Test_points, type = "response")
    MAE_CLC <- mean(abs(Preds_CLC - Test_points$Abundance))
    CLC_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_CLC, Real_Abundances = Test_points$Abundance)
    Results_CLC[[i]] <- CLC_model_perf
    confi_int <- confint(fit.glm_sp_CLC, predictors[which(estims[,"Pr(>|z|)"] < 0.05)], level = 0.95)
    for (predic in predictors){
      if (!(predic %in% rownames(confi_int))){
        confi_int <- rbind(confi_int,c(NA, NA))
        rownames(confi_int)[length(rownames(confi_int))] <- predic
      }
    }
    Confint_CLC[[i]] <- confi_int
    
    
    fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                 data = Train_points)
    estims <- summary(fit.glm_sp_dist_agri)[["coefficients"]]
    predictors <- rownames(estims)
    Estimates_dist_agri[[i]] <- estims
    Preds_dist_agri <- predict(fit.glm_sp_dist_agri, newdata = Test_points, type = "response")
    MAE_dist_agri <- mean(abs(Preds_dist_agri - Test_points$Abundance))
    dist_agri_model_perf <- data.frame(Site = Test_points$ID_point, Predictions = Preds_dist_agri, Real_Abundances = Test_points$Abundance)
    Results_dist_agri[[i]] <- dist_agri_model_perf
    confi_int <- confint(fit.glm_sp_dist_agri, predictors[which(estims[,"Pr(>|z|)"] < 0.05)], level = 0.95)
    for (predic in predictors){
      if (!(predic %in% rownames(confi_int))){
        confi_int <- rbind(confi_int,c(NA, NA))
        rownames(confi_int)[length(rownames(confi_int))] <- predic
      }
    }
    Confint_dist_agri[[i]] <- confi_int
    
    MAE_matrix[i,]  <- c(MAE_base, MAE_CLC, MAE_dist_agri)
    
  }
  MAE_df <- as.data.frame(MAE_matrix)
  names(MAE_df) <- c("MAE Base", "MAE CLC", "MAE Dist Agri")
  write.csv(MAE_df,paste("Results/Results_by_sp/",sp_cat,"/NB_models_MAE_",sp,"_no_CLC.csv", sep = ""))
  saveRDS(Estimates_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_base_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Estimates_CLC, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_CLC_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Estimates_dist_agri, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Estimates_dist_agri_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Results_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_base_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Results_CLC, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_CLC_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Results_dist_agri, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Results_dist_agri_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Confint_base, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_base_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Confint_CLC, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_CLC_",sp,"_no_CLC.RDS", sep = ""))
  saveRDS(Confint_dist_agri, paste("Results/Results_by_sp/",sp_cat,"/NB_models_Confint_dist_agri_",sp,"_no_CLC.RDS", sep = ""))
}


plot_bat_eff <- function(sp){
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                   bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat", "Contains_Agri_250"))
  pred_bat_0 <- pred_bat %>% filter(group == 0)
  pred_bat_1 <- pred_bat %>% filter(group == 1)
  pred_bat_total <- as.tibble(pred_bat_0) %>% 
    rename(Dist_1_bat_0 = x, conf.high_0 = conf.high, conf.low_0 = conf.low, predicted_0 = predicted) %>% 
    add_column(Dist_1_bat_1 = pred_bat_1$x, conf.high_1 = pred_bat_1$conf.high, conf.low_1 = pred_bat_1$conf.low, predicted_1 = pred_bat_1$predicted)
  plot <- ggplot(pred_bat_total) + 
    geom_ribbon(aes(Dist_1_bat_0*bat_dist_sd, ymin = conf.low_0, ymax = conf.high_0), fill = "#EE3B3B", alpha = 0.3) +
    geom_line(aes(Dist_1_bat_0*bat_dist_sd, predicted_0, color = "No")) + 
    geom_ribbon(aes(Dist_1_bat_1*bat_dist_sd, ymin = conf.low_1, ymax = conf.high_1), fill = "#00B2EE", alpha = 0.3) +
    geom_line(aes(Dist_1_bat_1*bat_dist_sd, predicted_1, color = "Yes")) +
    scale_color_manual(name = "CLC Agriculture in 250m", values = c("Yes" = "#009ACD", "No" = "#CD3333")) +
    ggtitle(paste("Marginal effect of building distance -", (code_sp_df %>% filter(code_sp == sp))$nom_fr)) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance") +
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(c(max(pred_bat_total$predicted_0), max(pred_bat_total$predicted_1)))))
  print(plot)
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
  
  fit.glm_sp_dist_nl <- glm.nb(Abundance~CLC_250 + Temp_range + Precip + altitude + NDVI_mean_250 + poly(Dist_1_bat_log,2, raw = T) + Contains_Forest_250 + Contains_Urban_250 + Contains_Agri_250 + Contains_Water_250, 
                               data = Sites_infos_train)
  Estimates_dist_nl <- summary(fit.glm_sp_dist_nl)[["coefficients"]]
  Preds_dist_nl <- predict(fit.glm_sp_dist_nl, newdata = Sites_infos_test, type = "response")
  dist_nl_model_perf <- data.frame(Site = Sites_infos_test$ID_point, Predictions = Preds_dist_nl, Real_Abundances = Sites_infos_test$Abundance)
  write.csv(Estimates_dist_nl, paste("Results/Preds_over_year/",sp_cat,"/",sp,"/",sp,"_Estimates_dist_nl.csv", sep = ""))
  write.csv(dist_nl_model_perf, paste("Results/Preds_over_year/",sp_cat,"/",sp,"/",sp,"_Results_dist_nl.csv", sep = ""))
}


#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoSNOW(cl)
clusterExport(cl, "Sites_infos")
clusterEvalQ(cl, c(library('lme4'), library('MASS'), library('tidyverse'), library('fs'), library('splines')))


sp_list_other <- read_lines("STOC/sp_list_other.txt")
sp_list_all <- read_lines("STOC/sp_list_all.txt")
sp_list_agri <- read_lines("STOC/sp_list_agri.txt")
pb <- txtProgressBar(min=1, max=length(sp_list_all), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


#foreach(num = 1:49, .options.snow = opts) %dopar% calculate_model_sp(sp_list_all[num], sp_cat = "all")
#foreach(num = 1:length(sp_list_all), .options.snow = opts) %dopar% calculate_model_sp(sp_list_all[num], sp_cat = "all")
foreach(num = 1:length(sp_list_all), .options.snow = opts) %dopar% calculate_model_sp_no_CLC(sp_list_all[num], sp_cat = "all")

#foreach(num = 1:length(sp_list_all), .options.snow = opts) %dopar% evaluate_model(sp_list_all[num], "2013", "2018", sp_cat = "all")

stopCluster(cl)
