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

Add_FDis_Site_ab_df <- function(Site_ab_sp_df, FDis_at_site){
  Site_ab_sp_df <- Site_ab_sp_df %>% dplyr::filter(code_sp %in% Traits_df$code_sp) %>% 
    left_join(FDis_at_site %>%dplyr::select(c(ID_point, sp_richn, fdis)), by = "ID_point")
}

Add_Shannon_Site_ab_df <- function(Site_ab_sp_df){
  Shannon_at_point <- ungroup(Site_ab_sp_df) %>% 
    dplyr::select(c(ID_point, code_sp, Abundance)) %>% 
    group_by(ID_point) %>% 
    summarise(Shannon = - sum(Abundance/sum(Abundance)*log(Abundance/sum(Abundance))))
  
  Site_ab_sp_df %>% left_join(Shannon_at_point, by = "ID_point") 
}


STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
code_sp_df <- STOC_sp_table %>% dplyr::select(code_sp, nom_fr)
raptor_list <- c("ACCNIS", "BUTBUT", "MILMIG", "ACCGEN", "CIRCYA", "ATHNOC", "FALSUB", "CIRAER","ASIOTU","STRALU", "CIRPYG", "MILMIL", "CIRGAL", "FALPER", "GYPFUL", "FALTIN", "PERAPI")

year = "2018"
STOC_data_year <- read.csv(paste("STOC/STOC_points_data_ready_",year,".csv",sep=""))   # Vient de Observations_df_construction
STOC_data_year <- STOC_data_year %>% filter(Abundance > 0) %>% filter(code_sp != "COLLIV")
# Assemblages venant de Create_assemblage_df
STOC_assemblage_all <- read.csv(paste("STOC/Assemblage_",year,"_all.csv",sep = ""))
STOC_assemblage_aquatic <- read.csv(paste("STOC/Assemblage_",year,"_aquatic.csv",sep = ""))
STOC_assemblage_raptor <- read.csv(paste("STOC/Assemblage_",year,"_raptor.csv",sep = ""))
STOC_assemblage_other <- read.csv(paste("STOC/Assemblage_",year,"_other.csv",sep = ""))
other_sp_names <- names(STOC_assemblage_other)[-1]

sp_all <- unique(STOC_data_year$code_sp)
all_sp_names <- unlist(lapply(sp_all, function(code){code_sp_df[which(code_sp_df$code_sp == code),]$nom_fr}))
Traits_df <- read.csv("STOC/Traits_standard.csv")   # Vient de PCA_Traits
Trophic_niche <- Traits_df %>% dplyr::select(c(code_sp, Trophic.Niche))
Trophic_niche$Trophic.Niche <- as.factor(Trophic_niche$Trophic.Niche)
FDis_at_site_all <- readRDS(paste("STOC/FDis_site_",year,".RDS" , sep = ""))$"functional_diversity_indices"
FDis_at_site_all <- rownames_to_column(as.data.frame(FDis_at_site_all),"ID_point")
FDis_at_site_aquatic <- readRDS(paste("STOC/FDis_site_",year,"_aquatic.RDS" , sep = ""))$"functional_diversity_indices"
FDis_at_site_aquatic <- rownames_to_column(as.data.frame(FDis_at_site_aquatic),"ID_point")
FDis_at_site_raptor <- readRDS(paste("STOC/FDis_site_",year,"_raptor.RDS" , sep = ""))$"functional_diversity_indices"
FDis_at_site_raptor <- rownames_to_column(as.data.frame(FDis_at_site_raptor),"ID_point")
FDis_at_site_other <- readRDS(paste("STOC/FDis_site_",year,"_other.RDS" , sep = ""))$"functional_diversity_indices"
FDis_at_site_other <- rownames_to_column(as.data.frame(FDis_at_site_other),"ID_point")


STOC_data_year <- STOC_data_year %>%
  dplyr::select(-c(X, NDVI_mean_500, NDVI_var_500, Bt_s_250, Bt_s_500)) %>%
  add_column(Dist_1_bat_log = log(STOC_data_year$Dist_1_bat))


# Moyenne des abondances sur les années d'observation
STOC_data_year <- ungroup(STOC_data_year) %>% 
  group_by(across(-c(annee, Abundance))) %>% 
  summarise(Abundance = round(mean(Abundance)))

# Ajout des distances fonctionnelles
Site_ab_sp_df <- Add_FDis_Site_ab_df(STOC_data_year, FDis_at_site_all) %>% rename(c(fdis_all = fdis, sp_richn_all = sp_richn))
Site_ab_sp_df_aquatic <- Add_FDis_Site_ab_df(Site_ab_sp_df %>% left_join(Trophic_niche, by = "code_sp") %>% filter(grepl("Aquatic",Trophic.Niche)), FDis_at_site_aquatic) %>% rename(c(fdis_aqua = fdis, sp_richn_aqua = sp_richn))
Site_ab_sp_df_raptor <- Add_FDis_Site_ab_df(Site_ab_sp_df %>% filter(code_sp %in% raptor_list), FDis_at_site_raptor) %>% rename(c(fdis_rapt = fdis, sp_richn_rapt = sp_richn))
excluded_species <- c(unique(Site_ab_sp_df_aquatic$code_sp),unique(Site_ab_sp_df_raptor$code_sp))
Site_ab_sp_df_other <- Add_FDis_Site_ab_df(Site_ab_sp_df %>% filter(!(code_sp %in% excluded_species)), FDis_at_site_other) %>% rename(c(fdis_other = fdis, sp_richn_other = sp_richn))

# Ajout des diversités de Shannon
Site_ab_sp_df <- Add_Shannon_Site_ab_df(Site_ab_sp_df) %>% rename(c(Shannon_all = Shannon))
Site_ab_sp_df_aquatic <- Add_Shannon_Site_ab_df(Site_ab_sp_df_aquatic) %>% rename(c(Shannon_aqua = Shannon)) %>% 
  left_join(ungroup(Site_ab_sp_df) %>%dplyr::select(c(ID_point, Shannon_all)) %>% group_by(ID_point,Shannon_all) %>% summarise(), by = "ID_point")
Site_ab_sp_df_raptor <- Add_Shannon_Site_ab_df(Site_ab_sp_df_raptor) %>% rename(c(Shannon_rapt = Shannon)) %>% 
  left_join(ungroup(Site_ab_sp_df) %>%dplyr::select(c(ID_point, Shannon_all)) %>% group_by(ID_point,Shannon_all) %>% summarise(), by = "ID_point")
Site_ab_sp_df_other <- Add_Shannon_Site_ab_df(Site_ab_sp_df_other) %>% rename(c(Shannon_other = Shannon)) %>% 
  left_join(ungroup(Site_ab_sp_df) %>%dplyr::select(c(ID_point, Shannon_all)) %>% group_by(ID_point,Shannon_all) %>% summarise(), by = "ID_point")


# Modèle linéaire FDis ~ Shannon -> utilisation des résidus par la suite
fit.diversity <- lm(fdis_other~Shannon_other , data = Site_ab_sp_df_other)
summary(fit.diversity)
ggplot(Site_ab_sp_df_other, aes(x = Shannon_other, y = fdis_other)) +
  geom_point() +
  stat_smooth(method = "lm", se=FALSE)
fdis_res <- fit.diversity$residuals
Site_ab_sp_df_other <- Site_ab_sp_df_other %>% add_column(fdis_residuals = fdis_res)

