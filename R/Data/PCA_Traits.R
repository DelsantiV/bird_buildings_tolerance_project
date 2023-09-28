setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)
library(terra)
library(factoextra)
library(ggplot2)
library(StatMatch)
library(ape)
library(PCAmixdata)

euclidean <- function(a, b) sqrt(sum((a - b)^2))
get_distance <- function(df_row, df) {apply(df,1,function(other_row) euclidean(df_row, other_row))}

# Charger les données
STOC_file = "STOC/data_FrenchBBS_point_Derminon_allSp_2001_2022fr.csv"
STOC_data = read.csv(STOC_file, sep = ';', dec = ',')
STOC_data <- STOC_data %>% filter(as.integer(annee) > 2005)
STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv")  # Vient de Filter_sp_and_traits
Traits_table_STOC <- read.csv("STOC/Traits_all_STOC_sp.csv")  # Vient de Filter_sp_and_traits


morpho_pca <- prcomp(Traits_table_STOC %>% select(c(LengthU_MEAN, TailU_MEAN, BillU_MEAN, TarsusU_MEAN, WeightU_MEAN)), scale = T)
summary(morpho_pca)
#fviz_pca_ind(morpho_pca)
fviz_pca_var(morpho_pca)
fviz_screeplot(morpho_pca, addlabels=T)

Traits_table_STOC <- Traits_table_STOC %>% mutate(PC1_morpho = morpho_pca$x[,1], PC2_morpho = morpho_pca$x[,2])

Traits_for_pca <- Traits_table_STOC %>% 
  select(c(code_sp,PC1_morpho, PC2_morpho, Clutch_MEAN, Sedentary,
           Nb_hab, Diet_gen, Nest.type,
           range_T, Age.of.first.breeding, Trophic.Niche, Lifestyle))
Traits_for_pca$Sedentary <- as.factor(Traits_for_pca$Sedentary)
Traits_for_pca$Nest.type <- as.factor(Traits_for_pca$Nest.type)
Traits_for_pca$Trophic.Niche <- as.factor(Traits_for_pca$Trophic.Niche)
Traits_for_pca$Lifestyle <- as.factor(Traits_for_pca$Lifestyle)

global_pca <- PCAmix( X.quanti = as.data.frame(Traits_for_pca %>% select(c(PC1_morpho, Clutch_MEAN, Nb_hab, Diet_gen, range_T, Age.of.first.breeding))),
                      X.quali = as.data.frame(Traits_for_pca %>% select(c(Sedentary, Nest.type, Trophic.Niche, Lifestyle))), rename.level =  T)
summary(global_pca)
global_pca$eig
plot(global_pca, choice = "cor")
plot(global_pca, choice = "ind")
plot(global_pca, choice = "sqload")
plot(global_pca, choice = "levels")


scaled_traits <- Traits_for_pca
scaled_traits$Clutch_MEAN = scale(scaled_traits$Clutch_MEAN)
scaled_traits$Nb_hab = scale(scaled_traits$Nb_hab)
scaled_traits$Diet_gen = scale(scaled_traits$Diet_gen)
scaled_traits$range_T = scale(scaled_traits$range_T)
scaled_traits$Age.of.first.breeding = scale(scaled_traits$Age.of.first.breeding)

scaled_traits_PCA_quanti <- scale(global_pca$ind$coord)
write.csv(scaled_traits, file = "STOC/Traits_standard.csv")
