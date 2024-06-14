setwd("D:/Projet woizos")
library(tidyverse)
library(fs)
library(ggplot2)
library(ggeffects)
library(MASS)
library(lme4)
library(splines)

year = "all_years"
STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
Site_ab_sp_df <- read.csv(paste("STOC/Assemblage_",year,".csv",sep=""))
sp_list <- readLines("STOC/sp_list_selected.txt")
code_sp_table <- read.csv("STOC/STOC_sp_code_table.csv")
Sites_infos <- read.csv(paste("STOC/Infos_at_site_",year,".csv",sep=""))
bat_dist_sd <- sd(Sites_infos$Dist_1_bat)


plot_marg_eff <- function(sp){
  sp_assemblage <- Site_ab_sp_df[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                 bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat"))
  plot <- plot(pred_bat) + 
    ggtitle((code_sp_table %>% filter(code_sp == sp))$nom_fr) +
    ylim(c(0, 1.5*max(pred_bat$predicted)))
  
  print(plot)
}

for (sp in sp_list){plot_marg_eff(sp)}


plot_bat_eff <- function(sp){
  sp_assemblage <- Site_ab_sp_df[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                   bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat", "Contains_Agri_250"))
  pred_bat_0 <- pred_bat %>% filter(group == 0)
  pred_bat_1 <- pred_bat %>% filter(group == 1)
  pred_bat_total <- as_tibble(pred_bat_0) %>% 
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