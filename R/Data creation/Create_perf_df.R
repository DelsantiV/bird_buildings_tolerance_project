setwd("C:/Users/vincent/Documents/Stage Vincent")
library(tidyverse)
library(fs)
library(MASS)
library(gridExtra)
library(lme4)
library(sjPlot)
library(ggeffects)
library(broom.helpers)
library(GGally)
library(ggplot2)
library(splines)

STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv")
code_sp_df <- STOC_sp_table %>% dplyr::select(code_sp, nom_fr)
Sites_infos <- read.csv("STOC/Infos_at_site_2018.csv")
bat_dist_sd <- sd(Sites_infos$Dist_1_bat)
Sites_infos <- Sites_infos %>% add_column(Dist_1_bat_red = Sites_infos$Dist_1_bat/bat_dist_sd)
attribution_matrix <- read.csv("STOC/Attribution_train_preds.csv")[,-1]
STOC_assemblage <- column_to_rownames(read.csv(paste("STOC/Assemblage_2018_all.csv",sep = "")),"X")
STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
sp_list_other <- read_lines("STOC/sp_list_other.txt")
sp_list<- read_lines("STOC/sp_list_all.txt")

for (sp in sp_list){
  sp_assemblage <- STOC_assemblage[,sp]
  
  Sites_infos_sp <- Sites_infos %>% add_column(Presence = as.factor(unlist(lapply(sp_assemblage, function(ab) as.integer(ab > 0.9)))))
  plot <- ggplot(data = Sites_infos_sp, aes(x = Dist_1_bat)) + 
    geom_histogram(aes(fill = Presence), position = "fill", color = "white", binwidth = 25) + 
    scale_fill_manual(values = c("grey80", "grey30"), labels = c("No", "Yes")) +
    theme_bw() +
    labs(x = "Distance to closest building", y = "Proportion of sites where species is observed") +
    ggtitle(paste(sp, "presence distribution regarding distance to closest building"))
  
  print(plot)
}

plot_bat_eff_CLC_agri <- function(sp){
  sp_assemblage <- STOC_assemblage[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                   bs(Dist_1_bat_red, degree = 3, knots = quantile(Sites_infos$Dist_1_bat_red, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat_red", "Contains_Agri_250"))
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
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(c(max(pred_bat_total$predicted_0), max(pred_bat_total$predicted_1))))) + 
    theme_bw()
  print(plot)
  ggsave(paste("Results/Figures/Marginal_effects/",sp,"_Agri.png",sep=""))
}

plot_bat_eff_NDVI <- function(sp){
  sp_assemblage <- STOC_assemblage[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                   bs(Dist_1_bat_red, degree = 3, knots = quantile(Sites_infos$Dist_1_bat_red, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat_red", "NDVI_mean_250"))
  pred_bat$x <- pred_bat$x*bat_dist_sd
  plot <- plot(pred_bat) + 
    ggtitle(paste("Marginal effect of building distance -", (code_sp_df %>% filter(code_sp == sp))$nom_fr)) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance",) +
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(pred_bat$predicted))) + 
    theme_bw()
  print(plot)
  ggsave(paste("Results/Figures/Marginal_effects/",sp,"_NDVI.png",sep=""))
}

plot_bat_eff <- function(sp){
  sp_assemblage <- STOC_assemblage[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                   bs(Dist_1_bat_red, degree = 3, knots = quantile(Sites_infos$Dist_1_bat_red, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat_red"))
  pred_bat$x <- pred_bat$x*bat_dist_sd
  plot <- plot(pred_bat) + 
    ggtitle(paste("Marginal effect of building distance -", (code_sp_df %>% filter(code_sp == sp))$nom_fr)) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance",) +
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(pred_bat$predicted))) + 
    theme_bw()
  print(plot)
  ggsave(paste("Results/Figures/Marginal_effects/",sp,".png",sep=""))
}

for (sp in sp_list){
  plot_bat_eff(sp)
  plot_bat_eff_NDVI(sp)
  plot_bat_eff_CLC_agri(sp)
}

sp_cat = "all_test"

NA_to_zero <- function(number){
  if (is.na(number)){to_return <- 0} else{to_return <- as.double(number)}
  to_return
}

nbchar <- function(string){
  if (is.na(string)){len <- 2} else{len <- nchar(string)}
  len
}

estim_sign <- function(conf_inter){
  if (any(is.na(conf_inter))){
    to_return <- "Non significant"
  }else if (conf_inter[1] < 0 & conf_inter[2] < 0){
    to_return <- "Negative"
  } else if (conf_inter[1] > 0 & conf_inter[2] >0){
    to_return <- "Positive"
  }
  else{
    to_return <- "Non significant"
  }
  to_return
}


get_model_infos <- function(Results_sp_df){
  if (length(Results_sp_df)>1){
    Classif_eval_list <- lapply(Results_sp_df, calculate_PPV_and_NPV)
    all_PPV <- unlist(lapply(Classif_eval_list, function(el) el[1]))
    all_NPV <- unlist(lapply(Classif_eval_list, function(el) el[2]))
    all_MCC <- unlist(lapply(Classif_eval_list, function(el) el[3]))
    nb_no_PP <- as.double(table(all_MCC)["Undefined (PP = 0)"])
    nb_no_PN <- as.double(table(all_MCC)["Undefined (PN = 0)"])
    all_MCC <- as.double(all_MCC[!grepl("Undefined", all_MCC)])
    PPV_mean <- mean(all_PPV)
    PPV_sd <- sd(all_PPV)
    NPV_mean <- mean(all_NPV)
    NPV_sd <- sd(all_NPV)
    MCC_mean <- mean(all_MCC)
    MCC_sd <- sd(all_MCC)
    c(PPV_mean, PPV_sd, NPV_mean, NPV_sd, MCC_mean, MCC_sd, NA_to_zero(nb_no_PP), NA_to_zero(nb_no_PN))
  } else {NA}
}


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


calculate_PPV_and_NPV <- function(results_df){
  Presence_sites <- results_df %>% dplyr::filter(Real_Abundances > 0)
  P <- as.double((dim(Presence_sites)[1]))
  PP <- as.double(dim(results_df %>% dplyr::filter(Predictions>0.5))[1])
  TP <- as.double(dim(Presence_sites %>% dplyr::filter(Predictions>0.5))[1])
  FP <- PP - TP
  TPR <- TP/P
  if (PP == 0){PPV = 0} else {PPV <- TP/PP}
  No_ab_sites <- results_df %>% dplyr::filter(Real_Abundances == 0)
  N <- as.double((dim(No_ab_sites)[1]))
  PN <- as.double(dim(results_df %>% dplyr::filter(Predictions<0.5))[1])
  TN <- as.double(dim(No_ab_sites %>% dplyr::filter(Predictions<0.5))[1])
  FN <- PN - TN
  TNR <- TN/N
  if (PN == 0){NPV = 0} else {NPV <- TN/PN}
  if (PP == 0){MCC <- "Undefined (PP = 0)"
  } else if (PN == 0){MCC <- "Undefined (PN = 0)"
  } else {MCC <- (TP*TN + FP*FN)/sqrt(PP * P * PN * N)}
  list(PPV,NPV,MCC)
}

find_CLC_estim_sp <- function(sp, model_type, sp_cat = "other", more_param = F){
  all_estims <- load_model_estim(sp,model_type, sp_cat = sp_cat, more_param = more_param)
  res <- list()
  all_param <- c("CLC_250Urban","CLC_250Forest and natural","CLC_250Water")
  if (more_param){all_param = c(all_param, "Contains_Urban_250", "Contains_Forest_250", "Contains_Water_250", "Contains_Agri_250")}
  estims_params <- lapply(all_estims, function(summary_estim) summary_estim[,1][all_param])
  i = 0
  for(param in all_param){
    i = i+1
    estim_list <- unlist(lapply(estims_params, function(el) el[param]))
    res[[i]] <- quantile(estim_list, c(0.025, 0.975))}
  names(res) <- all_param
  res
}

find_all_estim_sp <- function(sp, model_type, sp_cat = "other", more_param = F, CLC = T){
  all_estims <- load_model_estim(sp,model_type, sp_cat = sp_cat, more_param = more_param, CLC = CLC)
  all_param <- rownames(all_estims[[1]])
  res <- list()
  estims_params <- lapply(all_estims, function(summary_estim) summary_estim[,1][all_param])
  i = 0
  for(param in all_param){
    i = i + 1
    estim_list <- unlist(lapply(estims_params, function(el) el[param]))
    res[[i]] <- mean(estim_list)
  }
  names(res) <- all_param
  res
}


plot_model_estim <- function(sp, show_true_ab = F){
  sp_assemblage <- STOC_assemblage[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Real_Abundance = sp_assemblage)
  confint_bat1 <- get_pred_confint(sp,"dist_agri",sp_cat = sp_cat, more_param = F, CLC = F, predictor = "poly(Dist_1_bat, 2, raw = T)1")
  confint_bat2 <- get_pred_confint(sp,"dist_agri",sp_cat = sp_cat, more_param = F, CLC = F, predictor = "poly(Dist_1_bat, 2, raw = T)2")
  estim_bat1 <- find_all_estim_sp(sp,"dist_agri",sp_cat = sp_cat, more_param = F, CLC = F)["poly(Dist_1_bat, 2, raw = T)1"]
  estim_bat2 <- find_all_estim_sp(sp,"dist_agri",sp_cat = sp_cat, more_param = F, CLC = F)["poly(Dist_1_bat, 2, raw = T)2"]
  if (any(is.na(confint_bat1))){estim_bat1 = 0}
  if (any(is.na(confint_bat2))){estim_bat2 = 0}
  ab_vs_dist_mean <- NA_to_zero(estim_bat2)*Sites_infos_sp$Dist_1_bat_red^2 + NA_to_zero(estim_bat1)*Sites_infos_sp$Dist_1_bat_red
  ab_vs_dist_max <- NA_to_zero(confint_bat2[2])*(Sites_infos_sp$Dist_1_bat_red)^2 + NA_to_zero(confint_bat1[2])*Sites_infos_sp$Dist_1_bat_red
  ab_vs_dist_min <- NA_to_zero(confint_bat2[1])*Sites_infos_sp$Dist_1_bat_red^2 + NA_to_zero(confint_bat1[1])*Sites_infos_sp$Dist_1_bat_red
  Sites_infos_sp <- Sites_infos_sp %>% add_column(Predicted_Abundance_mean = ab_vs_dist_mean,
                                                  Predicted_Abundance_max = ab_vs_dist_max,
                                                  Predicted_Abundance_min = ab_vs_dist_min)
  if (show_true_ab){
    plot <- ggplot(as.data.frame(Sites_infos_sp)) + 
      geom_point(aes(x = Dist_1_bat, y = Real_Abundance), color = "green") + 
      geom_ribbon(aes(x = Dist_1_bat, y = Predicted_Abundance_mean), color = "red", shape = 4) +
      ggtitle(sp)
  } else{
    plot <- ggplot(Sites_infos_sp) + 
      geom_ribbon(aes(x = Dist_1_bat, ymin = Predicted_Abundance_min, ymax = Predicted_Abundance_max), fill = "grey70") +
      geom_line(aes(x = Dist_1_bat, y = Predicted_Abundance_mean)) +
      ggtitle(sp)
  }
  print(plot)
}


get_models_perf <- function(sp, Models_perf_df){
  results_sp_df_model_base <- load_model_results(sp,"base", sp_cat = sp_cat, more_param = F, CLC = F)
  results_sp_df_model_CLC <- load_model_results(sp,"CLC", sp_cat = sp_cat, more_param = F, CLC = F)
  results_sp_df_model_dist_agri <- load_model_results(sp,"dist_agri", sp_cat = sp_cat, more_param = F, CLC = F)
  Classif_eval_list_model_base <- lapply(results_sp_df_model_base, calculate_PPV_and_NPV)
  Classif_eval_list_model_CLC <- lapply(results_sp_df_model_CLC, calculate_PPV_and_NPV)
  Classif_eval_list_model_dist_agri <- lapply(results_sp_df_model_dist_agri, calculate_PPV_and_NPV)
  all_MCC_model_base <- unlist(lapply(Classif_eval_list_model_base, function(el) el[3]))
  all_MCC_model_CLC <- unlist(lapply(Classif_eval_list_model_CLC, function(el) el[3]))
  all_MCC_model_dist_agri <- unlist(lapply(Classif_eval_list_model_dist_agri, function(el) el[3]))
  MCC_base = mean(all_MCC_model_base)
  MCC_CLC = mean(all_MCC_model_CLC)
  MCC = mean(all_MCC_model_dist_agri)
  
  NPV_mean_model_base <- mean(unlist(lapply(Classif_eval_list_model_base, function(el) el[2])))
  NPV_mean_model_dist_agri <- mean(unlist(lapply(Classif_eval_list_model_dist_agri, function(el) el[2])))

  
  PPV_mean_model_base <- mean(unlist(lapply(Classif_eval_list_model_base, function(el) el[1])))
  PPV_mean_model_dist_agri <- mean(unlist(lapply(Classif_eval_list_model_dist_agri, function(el) el[1])))
  
  
  confint_Urb <- get_pred_confint(sp,"CLC",sp_cat = sp_cat, more_param = F, CLC = F, predictor = "CLC_250Urban")
  confint_Urb_contains <- get_pred_confint(sp,"CLC",sp_cat = sp_cat, more_param = F, CLC = F, predictor = "Contains_Urban_250")
  urban_CLC_impact <- estim_sign(confint_Urb)
  urban_Contains_impact <- estim_sign(confint_Urb_contains)
  
  if (urban_CLC_impact == urban_Contains_impact){
    sens_to_CLC <- urban_Contains_impact
  } else if (urban_CLC_impact == "Non significant"){
    sens_to_CLC <- urban_Contains_impact
  } else if (urban_Contains_impact == "Non significant"){
    sens_to_CLC <- urban_CLC_impact
  } else{sens_to_CLC <- "To determine"}
  

  
  if (is.na(MCC)) {model_qual <- "Bad"
  } else if (round(MCC,2) < 0.2){model_qual <- "Bad"
  } else {model_qual <- "Good"}
  
  MCC_amelioration_rel <- (MCC/MCC_base - 1)*100 
  MCC_amelioration_abs <- MCC - MCC_base
  
  if (is.na(MCC)){
    better_than_CLC = "No"
  } else if (is.na(MCC_CLC)){better_than_CLC = "Yes"
  } else if(MCC > MCC_CLC){better_than_CLC = "Yes"
  } else {better_than_CLC = "No"}
  
  if (is.na(MCC_amelioration_rel)) {
    if (!is.na(MCC)){Bat_impact = "Yes"
    } else{Bat_impact = "No"}
  } else if (MCC_amelioration_rel > 5){Bat_impact = "Yes"
  } else{Bat_impact = "No"}
  
  NPV_amelioration_rel <- round((NPV_mean_model_dist_agri/NPV_mean_model_base - 1)*100) 
  NPV_amelioration_abs <- round((NPV_mean_model_dist_agri - NPV_mean_model_base)*100) 
  PPV_amelioration_rel <- round((PPV_mean_model_dist_agri/PPV_mean_model_base - 1)*100) 
  PPV_amelioration_abs <- round((PPV_mean_model_dist_agri - PPV_mean_model_base)*100)
  
  Models_perf_df <- Models_perf_df %>% 
    add_row(Code_sp = sp, 
            Nom_fr = (STOC_sp_table %>% filter(code_sp ==sp))$nom_fr, 
            Better.model.with.building = Bat_impact, 
            Model.quality = model_qual,
            Response.to.urban.CLC = sens_to_CLC,
            Better.than.CLC = better_than_CLC,
            MCC = round(MCC,3),
            MCC.Amelioration.Rel = round(MCC_amelioration_rel,1),
            MCC.Amelioration.Abs = round(MCC_amelioration_abs,3),
            NPV = round(NPV_mean_model_dist_agri*100),
            NPV.Amelioration.Rel = NPV_amelioration_rel,
            NPV.Amelioration.Abs = NPV_amelioration_abs,
            PPV = round(PPV_mean_model_dist_agri*100),
            PPV.Amelioration.Rel = PPV_amelioration_rel,
            PPV.Amelioration.Abs = PPV_amelioration_abs)
  Models_perf_df
}



Col_names <- c("Code_sp",
               "Nom_fr",
               "Better.model.with.building",
               "Model.quality", 
               "Response.to.building.distance",
               "Response.to.urban.CLC",
               "Better.than.CLC",
               "Remark",
               "MCC",
               "MCC.Amelioration.Rel",
               "MCC.Amelioration.Abs", 
               "NPV", 
               "NPV.Amelioration.Rel", 
               "NPV.Amelioration.Abs", 
               "PPV", 
               "PPV.Amelioration.Rel", 
               "PPV.Amelioration.Abs")
Models_perf_df <- as.data.frame(matrix(ncol = length(Col_names)))
names(Models_perf_df) <- Col_names

for (sp in sp_list){
  Models_perf_df <- get_models_perf(sp, Models_perf_df)
}
Models_perf_df <- Models_perf_df %>% filter(!is.na(Code_sp))
write.csv(Models_perf_df, "Results/Models_perf_test.csv")

pb_sp <- (Models_perf_df %>% filter(is.na(MCC)))$Nom_fr

Models_perf_df_simple <- Models_perf_df %>% filter(!is.na(MCC)) %>% dplyr::select(c(Code_sp, Nom_fr, Model.quality, Better.model.with.building))
write.csv(Models_perf_df_simple, "Results/Models_perf_simple.csv")