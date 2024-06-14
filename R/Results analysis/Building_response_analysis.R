setwd("D:/Projet woizos")
library(tidyverse)
library(fs)
library(MASS)
library(lme4)
library(ggeffects)
library(glmm.hp)
library(splines)

Sites_infos <- read.csv("STOC/Infos_at_site_2018.csv")
bat_dist_sd <- sd(Sites_infos$Dist_1_bat)
Sites_infos <- Sites_infos %>% add_column(Dist_1_bat_red = Sites_infos$Dist_1_bat/bat_dist_sd)
attribution_matrix <- read.csv("STOC/Attribution_train_preds.csv")[,-1]
STOC_assemblage <- column_to_rownames(read.csv(paste("STOC/Assemblage_2018_all.csv",sep = "")),"X")
STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv")

Type_table <- read.csv("STOC/group_species.csv") %>% 
  rename(Nom_fr = Name) %>% 
  dplyr::select(c(Nom_fr, Groupe))

load_model_estim <- function(sp,model_type, more_param = F, sp_cat = "other", CLC = T){
  if (more_param){file = paste("NB_models_Estimates_",model_type,"_",sp,"_more_param.RDS",sep="")
  } else if (!CLC){file = paste("NB_models_Estimates_",model_type,"_",sp,"_no_CLC.RDS",sep="")
  } else {file = paste("NB_models_Estimates_",model_type,"_",sp,".RDS",sep="")}
  if (fs::file_exists(paste("Results/Results_by_sp/",sp_cat,"/",file, sep =""))){
    readRDS(paste("Results/Results_by_sp/",sp_cat,"/",file, sep =""))
  }
  else{
    cat("No file found for ",model_type," model and species ",sp)
    NA
  }
}

load_model_results <- function(sp, model_type, more_param = F, sp_cat = "other", CLC = T){
  if (more_param){file = paste("NB_models_Results_",model_type,"_",sp,"_more_param.RDS",sep="")
  } else if (!CLC){file = paste("NB_models_Results_",model_type,"_",sp,"_no_CLC.RDS",sep="")
  } else {file = paste("NB_models_Results_",model_type,"_",sp,".RDS",sep="")}
  if (fs::file_exists(paste("Results/Results_by_sp/",sp_cat,"/",file, sep =""))){
    readRDS(paste("Results/Results_by_sp/",sp_cat,"/",file, sep =""))
  }
  else{
    #cat("No file found for ",model_type," model and species ",sp,"\n")
    NA
  }
}

get_model_all_mean_estims <- function(sp,model_type, more_param = F, sp_cat = "other", CLC = T){
  all_estims <- load_model_estim(sp,model_type, more_param = more_param, sp_cat = sp_cat, CLC = CLC)
  all_params <- row.names(all_estims[[1]])[-1]
  res <- list()
  i <- 0
  for (param in all_params){
    i <- i + 1
    param_mean <- mean(unlist(lapply(all_estims, function(el) el[param,]["Estimate"])))
    param_p_val <- mean(unlist(lapply(all_estims, function(el) el[param,]["Pr(>|z|)"])))
    res[[i]] <- param_mean
  }
  names(res) <- all_params
  res
}

load_model_confint <- function(sp, model_type, more_param = F, sp_cat = "other", CLC = T){
  if (more_param){file = paste("NB_models_Confint_",model_type,"_",sp,"_more_param.RDS",sep="")
  } else if (!CLC){file = paste("NB_models_Confint_",model_type,"_",sp,"_no_CLC.RDS",sep="")
  } else {file = paste("NB_models_Confint_",model_type,"_",sp,".RDS",sep="")}
  if (fs::file_exists(paste("Results/Results_by_sp/",sp_cat,"/",file, sep =""))){
    readRDS(paste("Results/Results_by_sp/",sp_cat,"/",file, sep =""))
  }
  else{
    #cat("No file found for ",model_type," model and species ",sp,"\n")
    NA
  }
}

get_pred_confint <- function(sp, model_type, more_param = F, sp_cat = "all", predictor = "poly(Dist_1_bat_log, 2)1", CLC = T){
  all_confint <- lapply(load_model_confint(sp, model_type, more_param = more_param, sp_cat = sp_cat, CLC = CLC), function(el) el[predictor,])
  low_bound <- min(unlist(lapply(all_confint, function(el) el[1])))
  names(low_bound) <- "Min"
  up_bound <- max(unlist(lapply(all_confint, function(el) el[2])))
  names(up_bound) <- "Max"
  c(low_bound, up_bound)
}


get_pred_info <- function(sp, Sites_infos_sp, pred = "NDVI_mean_250", model_type = "base"){
  fit.glm_sp_base <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250, 
                         data = Sites_infos_sp)
  fit.glm_sp <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                           bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                         data = Sites_infos_sp)
  coef_contrib1 <- glmm.hp(fit.glm_sp_base)$delta[pred,]["I.perc(%)"]
  coef_contrib2 <- glmm.hp(fit.glm_sp)$delta[pred,]["I.perc(%)"]
  coef_estim <- get_model_all_mean_estims(sp, model_type, sp_cat = "all_test", CLC = F)[pred]
  coef_confint <- get_pred_confint(sp, model_type, sp_cat = "all_test", CLC = F, predictor = pred)
  c(coef_estim, coef_confint, coef_contrib1, coef_contrib2)
}

Models_perf_df <- read.csv("Results/Models_perf_with_response_and_dist.csv")[,-1] %>% 
  filter(Model.quality == "Good") %>%
  rename(code_sp = Code_sp) %>%
  mutate(Response.to.building.corrected = Response.to.building.distance) %>% 
  mutate_at(vars(c('Response.to.building.corrected')), ~if_else(Better.model.with.building == "No", "Indifferent",.)) %>%
  left_join(Type_table, by = "Nom_fr")
sp_list <- Models_perf_df$code_sp
sp_list_not_better <- (Models_perf_df %>% filter(Better.model.with.building == "No"))

Col_names <- c("Code_sp",
               "Nom_fr",
               "Response.to.building.corrected",
               "NDVI_mean_estim",
               "NDVI_mean_min",
               "NDVI_mean_max",
               "NDVI_mean_contrib_base",
               "NDVI_mean_contrib_dist")
NDVI_impact <- as.data.frame(matrix(ncol = length(Col_names)))
names(NDVI_impact) <- Col_names

compt = 0
for (sp in Models_perf_df$code_sp[30:43]){
  compt <- compt+1
  cat(sp," (",compt,"/",length(Models_perf_df$code_sp[30:43]),")\n\n")
  sp_assemblage <- STOC_assemblage[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  NDVI_info <- get_pred_info(sp, Sites_infos_sp)
  NDVI_impact <- NDVI_impact %>% 
    add_row(Code_sp = sp,
            Nom_fr = (STOC_sp_table %>% filter(code_sp == sp))$nom_fr,
            Response.to.building.corrected = (Models_perf_df %>% filter(code_sp == sp))$Response.to.building.corrected,
            NDVI_mean_estim = NDVI_info[1],
            NDVI_mean_min = NDVI_info[2],
            NDVI_mean_max = NDVI_info[3],
            NDVI_mean_contrib_base = NDVI_info[4],
            NDVI_mean_contrib_dist = NDVI_info[5])
}
NDVI_impact <- NDVI_impact[-1,]
write.csv(NDVI_impact, "Results/NDVI_impact2.csv")



NDVI_mean_conf_int <- list()
NDVI_var_conf_int <- list()
Contains_Agri_conf_int <- list()
i <- 0
for (sp in sp_list){
  i <- i + 1
  NDVI_mean_conf_int[[i]] <- get_pred_confint(sp, "base", sp_cat = "all_test", CLC = F, predictor = "NDVI_mean_250")
  NDVI_var_conf_int[[i]] <- get_pred_confint(sp, "base", sp_cat = "all_test", CLC = F, predictor = "NDVI_var_250")
  Contains_Agri_conf_int[[i]] <- get_pred_confint(sp, "base", sp_cat = "all_test", CLC = F, predictor = "Contains_Agri_250")
}

Models_perf_detailed_df <- Models_perf_df %>% add_column(NDVI_mean_min = unlist(lapply(NDVI_mean_conf_int, function(el) el[[1]])), 
                                                NDVI_mean_max = unlist(lapply(NDVI_mean_conf_int, function(el) el[[2]])), 
                                                NDVI_var_min = unlist(lapply(NDVI_var_conf_int, function(el) el[[1]])), 
                                                NDVI_var_max = unlist(lapply(NDVI_var_conf_int, function(el) el[[2]])), 
                                                Contains_Agri_max = unlist(lapply(Contains_Agri_conf_int, function(el) el[[1]])), 
                                                Contains_Agri_min = unlist(lapply(Contains_Agri_conf_int, function(el) el[[2]])))
write.csv(Models_perf_detailed_df,"Results/Good_sp_detailed_estim.csv")

for (sp in (Models_perf_detailed_df %>% filter(Response.to.building.corrected == "Indifferent"))$code_sp){
  mean_preds <- lapply(get_model_all_mean_estims(sp, "base", sp_cat = "all_test", CLC = F),abs)
  sp_assemblage <- STOC_assemblage[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                   bs(Dist_1_bat_red, degree = 3, knots = quantile(Sites_infos$Dist_1_bat_red, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  cat(sp,":",names(which.max(mean_preds)),"\n",anova(fit.glm_sp_dist_agri),"\n\n")
}