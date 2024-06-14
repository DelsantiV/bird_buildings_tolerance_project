setwd("D:/Projet woizos")
library(tidyverse)
library(fs)
library(ggplot2)
library(MASS)
library(lme4)
library(splines)


# Load 2008 data
Sites_infos_2008 <- read.csv(paste("STOC/Infos_at_site_2008.csv",sep=""))
Sites_infos_2008$ID_carre <- as.character(Sites_infos_2008$ID_carre)
bat_dist_sd_2008 <- sd(Sites_infos_2008$Dist_1_bat)
Sites_infos_2008 <- Sites_infos_2008 %>% 
  mutate(
  Dist_1_bat_norm = Sites_infos_2008$Dist_1_bat/bat_dist_sd_2008,
  ID_point = unlist(lapply(Sites_infos_2008$ID_point, function(el) paste0(el,"_2008", collapse="")))) 
Assemblage_df_2008 <- read.csv(paste("STOC/Assemblage_2008_all.csv",sep="")) %>% rename(ID_point = X)

# Load 2013 data
Sites_infos_2013 <- read.csv(paste("STOC/Infos_at_site_2013.csv",sep=""))
Sites_infos_2013$ID_carre <- as.character(Sites_infos_2013$ID_carre)
bat_dist_sd_2013 <- sd(Sites_infos_2013$Dist_1_bat)
Sites_infos_2013 <- Sites_infos_2013 %>% 
  mutate(
    Dist_1_bat_norm = Sites_infos_2013$Dist_1_bat/bat_dist_sd_2013,
    ID_point = unlist(lapply(Sites_infos_2013$ID_point, function(el) paste0(el,"_2013", collapse="")))) 
Assemblage_df_2013 <- read.csv(paste("STOC/Assemblage_2013_all.csv",sep="")) %>% rename(ID_point = X)


# Load 2018 data
Sites_infos_2018 <- read.csv(paste("STOC/Infos_at_site_2018.csv",sep=""))
bat_dist_sd_2018 <- sd(Sites_infos_2018$Dist_1_bat)
Sites_infos_2018 <- Sites_infos_2018 %>% 
  mutate(
    Dist_1_bat_norm = Sites_infos_2018$Dist_1_bat/bat_dist_sd_2018,
    ID_point = unlist(lapply(Sites_infos_2018$ID_point, function(el) paste0(el,"_2018", collapse="")))) 
Assemblage_df_2018 <- read.csv(paste("STOC/Assemblage_2018_all.csv",sep="")) %>% rename(ID_point = X)


# Aggregate data
Sites_infos <- bind_rows(list(Sites_infos_2008,Sites_infos_2013,Sites_infos_2018))
Assemblage_df <- bind_rows(list(Assemblage_df_2008, Assemblage_df_2013, Assemblage_df_2018))
write.csv(Sites_infos, "STOC/Infos_at_site_all_years.csv")
write.csv(Assemblage_df, "STOC/Assemblage_all_years.csv")