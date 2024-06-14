setwd("C:/Users/vincent/Documents/Stage Vincent")

library(tidyverse)
library(fs)
library(ggplot2)
library(ggeffects)
library(MASS)
library(lme4)
library(splines)
library(gridExtra)
library(ggpubr)

year = "2018"
STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
code_sp_df <- STOC_sp_table %>% dplyr::select(code_sp, nom_fr, nom)
Sites_infos <- read.csv(paste("STOC/Infos_at_site_",year,".csv",sep=""))
bat_dist_sd <- sd(Sites_infos$Dist_1_bat)
Sites_infos <- Sites_infos %>% mutate(Dist_1_bat = Sites_infos$Dist_1_bat/bat_dist_sd) 
Site_ab_sp_df_all <- read.csv(paste("STOC/Assemblage_",year,"_all.csv",sep=""))
sp_list <- read_lines("STOC/sp_list_all.txt")

Models_perf_df <- read.csv("Results/Models_perf_with_response_and_dist.csv")[,-1] %>% 
  filter(Model.quality == "Good") %>%
  rename(code_sp = Code_sp)

plot_bat_eff_CLC_agri <- function(sp){
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
    scale_color_manual(name = "Farmland Land Cover in 250m", values = c("Yes" = "#009ACD", "No" = "#CD3333")) +
    ggtitle((code_sp_df %>% filter(code_sp == sp))$nom) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance") +
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(c(max(pred_bat_total$predicted_0), max(pred_bat_total$predicted_1))))) +
    theme_bw() +
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=9), 
          title=element_text(size=10),
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(size=9),
          legend.text = element_text(size=8))
  print(plot)
}

plot_bat_eff_NDVI <- function(sp, fig_title = ""){
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                   bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat", "NDVI_mean_250"))
  pred_bat$x <- pred_bat$x*bat_dist_sd
  if (nchar(fig_title) == 0){fig_title = (code_sp_df %>% filter(code_sp == sp))$nom}
  plot <- plot(pred_bat) + 
    ggtitle(fig_title) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance",) +
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(pred_bat$predicted))) +
    theme_bw() +
    theme(axis.text=element_text(size=8),axis.title=element_text(size=9), title=element_text(size=10) )
  plot
}


plot_bat_eff_agri_inter <- function(sp, fig_title = ""){
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 * 
                                   bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat", "Contains_Agri_250"))
  pred_bat$x <- pred_bat$x*bat_dist_sd
  if (nchar(fig_title) == 0){fig_title = (code_sp_df %>% filter(code_sp == sp))$nom}
  plot <- plot(pred_bat) + 
    ggtitle(fig_title) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance",) +
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(pred_bat$predicted))) +
    theme_bw() +
    theme(axis.text=element_text(size=8),axis.title=element_text(size=9), title=element_text(size=10) )
  plot
}



plot_bat_eff_NDVI_inter <- function(sp, fig_title = ""){
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + Contains_Agri_250 + NDVI_var_250 + NDVI_mean_250 * 
                                   bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat", "NDVI_mean_250"))
  pred_bat$x <- pred_bat$x*bat_dist_sd
  if (nchar(fig_title) == 0){fig_title = (code_sp_df %>% filter(code_sp == sp))$nom}
  plot <- plot(pred_bat) + 
    ggtitle(fig_title) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance",) +
    xlim(c(0,2000)) + ylim(c(0, 1.1*max(pred_bat$predicted))) +
    theme_bw() +
    theme(axis.text=element_text(size=8),axis.title=element_text(size=9), title=element_text(size=10) )
  plot
}



plot_bat_eff <- function(sp, save_fig = F, fig_title = ""){
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + 
                                   bs(Dist_1_bat, degree = 3, knots = quantile(Sites_infos$Dist_1_bat, probs = seq(0, 1, 0.25))), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat"))
  pred_bat$x <- pred_bat$x*bat_dist_sd
  if (nchar(fig_title) == 0){fig_title = (code_sp_df %>% filter(code_sp == sp))$nom}
  plot <- plot(pred_bat) + 
    ggtitle(fig_title) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance",) +
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(pred_bat$predicted))) +
    theme_bw() +
    theme(axis.text=element_text(size=8),axis.title=element_text(size=9), title=element_text(size=10) )
  if (save_fig){
    ggsave(paste("Results/Figures/Marginal_effects/",sp,".png",sep=""), plot = plot)
  }
  plot
}

plot_bat_eff_old <- function(sp, save_fig = F, fig_title = ""){
  sp_assemblage <- Site_ab_sp_df_all[,sp]
  Sites_infos_sp <- Sites_infos %>% add_column(Abundance = sp_assemblage)
  fit.glm_sp_dist_agri <- glm.nb(Abundance~ Temp_range + Precip + altitude + NDVI_mean_250 + NDVI_var_250 + Contains_Agri_250 + poly(Dist_1_bat, 2), 
                                 data = Sites_infos_sp)
  pred_bat <- ggpredict(fit.glm_sp_dist_agri, terms = c("Dist_1_bat"))
  pred_bat$x <- pred_bat$x*bat_dist_sd
  if (nchar(fig_title) == ""){fig_title = (code_sp_df %>% filter(code_sp == sp))$nom}
  plot <- plot(pred_bat) + 
    ggtitle(ggtitle(fig_title)) +
    labs(x = "Distance to closest building", y = "Preddicted Abundance",) +
    xlim(c(0,2000)) + ylim(c(0, 1.5*max(pred_bat$predicted))) +
    theme_bw() + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=9), title=element_text(size=10) )
  if (save_fig){
    ggsave(paste("Results/Figures/Marginal_effects/",sp,"_poly.png",sep=""), plot = plot)
  }
  plot
}

plot1 <- plot_bat_eff("PASDOM")
plot2 <- plot_bat_eff("PHOOCH")
plot3 <- plot_bat_eff("HIRRUS")
plot4 <- plot_bat_eff("CARCHL")
multi_plot <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)
ggsave("Results/Figures/Building_marginal_effect_dwellers.png",multi_plot, width = 6, height = 6)


plot1 <- plot_bat_eff("PHYCOL")
plot2 <- plot_bat_eff("EMBCIT")
plot3 <- plot_bat_eff("BUTBUT")
plot4 <- plot_bat_eff("ANTTRI")
multi_plot <- ggarrange(plotlist = list(plot1, plot2, plot3, plot4), ncol = 2, nrow = 2)
ggsave("Results/Figures/Building_marginal_effect_avoiders.png",multi_plot, width = 6, height = 6)


plot1 <- plot_bat_eff("FALTIN")
plot2 <- plot_bat_eff("SAXTOR")
plot3 <- plot_bat_eff("LULARB")
plot4 <- plot_bat_eff("EMBCIR")
multi_plot <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)
ggsave("Results/Figures/Building_marginal_effect_optimal.png",multi_plot, width = 6, height = 6)


plot1 <- plot_bat_eff_CLC_agri("FALTIN")
plot2 <- plot_bat_eff_CLC_agri("SAXTOR")
plot3 <- plot_bat_eff_CLC_agri("LULARB")
plot4 <- plot_bat_eff_CLC_agri("EMBCIR")
multi_plot <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2, common.legend = T, legend = "bottom")
ggsave("Results/Figures/Building_marginal_effect_optimal_agri.png",multi_plot, width = 6, height = 6)


double_plot <- function(sp){
  plot1 <- plot_bat_eff(sp, fig_title = "Cubic Splines")
  plot2 <- plot_bat_eff_old(sp, fig_title = "Quadratic polynomial")
  multi_plot <- ggarrange(plot1, plot2, ncol = 2)
  annotate_figure(multi_plot, top = text_grob((code_sp_df %>% filter(code_sp == sp))$nom, size = 14))
  ggsave(paste("Results/Figures/Building_marginal_effect_double",sp,".png",sep=""), width = 6, height = 3)
}

positive_sp <- (Models_perf_df %>% filter(Better.model.with.building == "Yes" & Response.to.building.distance == "Positive"))$code_sp
negative_sp <- (Models_perf_df %>% filter(Better.model.with.building == "Yes" & Response.to.building.distance == "Negative"))$code_sp
other_sp <- (Models_perf_df %>% filter(Better.model.with.building == "Yes" & Response.to.building.distance == "Other"))$code_sp


plot_list <- list()
i = 0
for (sp in positive_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff(sp)
}
multi_plot <- ggarrange(plotlist = plot_list, nrow = 3, ncol = 5)
ggsave("Article figures/Building_marginal_effect_dwellers.png",multi_plot, width = 10, height = 6)

plot_list <- list()
i = 0
for (sp in negative_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff(sp)
}
multi_plot <- ggarrange(plotlist = plot_list)
ggsave("Article figures/Building_marginal_effect_avoiders.png",multi_plot, width = 6, height = 4)

plot_list <- list()
i = 0
for (sp in other_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff(sp)
}
multi_plot <- ggarrange(plotlist = plot_list)
ggsave("Article figures/Building_marginal_effect_optimal.png",multi_plot, width = 4, height = 4)

plot_list <- list()
i = 0
for (sp in other_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff_CLC_agri(sp)
}
multi_plot <- ggarrange(plotlist = plot_list, common.legend = T, legend = "bottom")
ggsave("Article figures/Building_marginal_effect_optimal_agri.png",multi_plot, width = 4, height = 4)

plot_list <- list()
i = 0
for (sp in other_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff_NDVI(sp)
}
multi_plot <- ggarrange(plotlist = plot_list, common.legend = T, legend = "bottom")
ggsave("Article figures/Building_marginal_effect_optimal_ndvi.png",multi_plot, width = 4, height = 4)


plot_list <- list()
i = 0
for (sp in negative_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff_CLC_agri(sp)
}
multi_plot <- ggarrange(plotlist = plot_list, common.legend = T, legend = "bottom")
ggsave("Article figures/Building_marginal_effect_negative_agri.png",multi_plot, width = 6, height = 4)

plot_list <- list()
i = 0
for (sp in negative_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff_NDVI(sp)
}
multi_plot <- ggarrange(plotlist = plot_list, common.legend = T, legend = "bottom")
ggsave("Article figures/Building_marginal_effect_negative_ndvi.png",multi_plot, width = 6, height = 4)


plot_list <- list()
i = 0
for (sp in positive_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff_CLC_agri(sp)
}
multi_plot <- ggarrange(plotlist = plot_list, common.legend = T, legend = "bottom", nrow = 3, ncol = 5)
ggsave("Article figures/Building_marginal_effect_positive_agri.png",multi_plot, width = 10, height = 6)

plot_list <- list()
i = 0
for (sp in positive_sp){
  i = i+1
  plot_list[[i]] <- plot_bat_eff_NDVI(sp)
}
multi_plot <- ggarrange(plotlist = plot_list, common.legend = T, legend = "bottom", nrow = 3, ncol = 5)
ggsave("Article figures/Building_marginal_effect_positive_ndvi.png",multi_plot, width = 10, height = 6)