setwd("D:/Projet woizos")
library(tidyverse)
library(fs)



sp_list <- readLines("STOC/sp_list_selected.txt")
code_sp_table <- read.csv("STOC/STOC_sp_code_table.csv")

source("Clean code/Loading_utilities.R")

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

get_model_infos <- function(Results_sp_df){
  if (length(Results_sp_df)>1){
    Classif_eval_list <- lapply(Results_sp_df, 
                                calculate_PPV_and_NPV)
    Spearm_cor <- mean(unlist(lapply(Results_sp_df, 
                                     function(Result_df)
                                       {cor(Result_df$Predictions, as.numeric(Result_df$Real_Abundances), method = "spearman")}
                                     )))
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
    c(PPV_mean, 
      PPV_sd, 
      NPV_mean, 
      NPV_sd, 
      MCC_mean, 
      MCC_sd, 
      NA_to_zero(nb_no_PP), 
      NA_to_zero(nb_no_PN),
      Spearm_cor)
  } else {NA}
}

Summary_df <- data.frame(matrix(ncol = 12, nrow = 0))
names(Summary_df) <- c("Code_sp", 
                       "Nom_fr", 
                       "Model quality", 
                       "Sensitivity to buildings", 
                       "MCC_base", 
                       "MCC_dist",
                       "Spearman_base",
                       "Spearman_dist",
                       "PPV_base", 
                       "PPV_dist", 
                       "NPV_base", 
                       "NPV_dist")
for (sp in sp_list){
  results_df_dist_agri <- load_model_results(sp)
  results_df_dist_agri <- lapply(results_df_dist_agri, 
                                 function(results_df){results_df %>% 
                                     mutate(Predictions = unlist(
                                       lapply(results_df$Predictions, 
                                              NA_to_zero)))})
  results_df_base <- load_model_results(sp,"base")
  results_df_base <- lapply(results_df_base, 
                                 function(results_df){results_df %>% 
                                     mutate(Predictions = unlist(
                                       lapply(results_df$Predictions, 
                                              NA_to_zero)))})
  dist_model_infos <- get_model_infos(results_df_dist_agri)
  PPV_dist <- dist_model_infos[1]
  NPV_dist <- dist_model_infos[3]
  MCC_dist <- dist_model_infos[5]
  Sp_cor_dist <- dist_model_infos[9]
  base_model_infos <- get_model_infos(results_df_base)
  PPV_base <- base_model_infos[1]
  NPV_base <- base_model_infos[3]
  MCC_base <- base_model_infos[5]
  Sp_cor_base <- base_model_infos[9]
  model_quality <- ifelse(
    MCC_dist > 0.2, 
    "Good",
    "Bad")
  building_sensi <- ifelse(
    MCC_dist/MCC_base > 1.05, 
    "Yes", 
    "No")
  Summary_df[nrow(Summary_df)+1,] <- c(sp, 
                                       (code_sp_table %>% filter(code_sp == sp))$nom_fr, 
                                       model_quality, 
                                       building_sensi, 
                                       MCC_base, 
                                       MCC_dist,
                                       Sp_cor_base,
                                       Sp_cor_dist,
                                       PPV_base, 
                                       PPV_dist, 
                                       NPV_base, 
                                       NPV_dist)
}
write.csv(Summary_df,"Results analysis/Models_summary.csv")
