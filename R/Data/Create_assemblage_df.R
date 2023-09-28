setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)

year = "2008"
STOC_data_year = read.csv(paste("STOC/STOC_points_data_ready_",year,".csv",sep=""))  # Vient de Observations_df_construction
all_sp <- unique(STOC_data_year$code_sp)
Traits_df = read.csv("STOC/Traits_all_STOC_sp.csv")    # Vient de Filter_sp_and_traits
Trophic_niche <- Traits_df %>% dplyr::select(c(code_sp, Trophic.Niche))

# Identification of aquatic species
Trophic_niche$Trophic.Niche <- as.factor(Trophic_niche$Trophic.Niche)
aquatic_list <- (Trophic_niche %>% filter(grepl("Aquatic", Trophic.Niche)))$code_sp
aquatic_list <- aquatic_list[aquatic_list %in% all_sp]

# Identification of raptors
raptor_list <- c("ACCNIS", "BUTBUT", "MILMIG", "ACCGEN", "CIRCYA", "ATHNOC", "FALSUB", "CIRAER","ASIOTU","STRALU", "CIRPYG", "MILMIL", "CIRGAL", "FALPER", "GYPFUL", "FALTIN", "PERAPI")
raptor_list <- raptor_list[raptor_list %in% all_sp]

scaled_traits = read.csv("STOC/Traits_standard.csv")[,-1] # Vient de PCA_Traits

sp_ab_df_total <- STOC_data_year %>% dplyr::select(c(ID_point, code_sp)) %>% group_by(ID_point, code_sp) %>% summarise(Abundance = n())
Points_df_total <- sp_ab_df_total %>% group_by(ID_point) %>% summarise(ID_point = unique(ID_point)) %>% column_to_rownames("ID_point")

Site_ab_sp_df_aquatic <- sp_ab_df_total %>% filter(code_sp %in% aquatic_list)

Site_ab_sp_df_raptor <- sp_ab_df_total %>% filter(code_sp %in% raptor_list)

excluded_species <- c(raptor_list, aquatic_list)
other_list <- all_sp[!all_sp %in% excluded_species]
Site_ab_sp_df_other <- sp_ab_df_total %>% filter(!(code_sp %in% excluded_species))


create_assemblage_df <- function(Points_df, Site_ab_sp_df){
  points_list <- rownames(Points_df)
  for (sp in names(Points_df)){
    Ab_at_point_sp <- Site_ab_sp_df %>% dplyr::filter(code_sp == sp) %>% dplyr::select(c(ID_point, Abundance))
    Points_df[points_list %in% Ab_at_point_sp$ID_point,][,sp] <- Ab_at_point_sp$Abundance 
  }
  Points_df
}

Points_df_total[,all_sp] <- 0
Points_df_total <- create_assemblage_df(Points_df_total, sp_ab_df_total)
write.csv(Points_df_total, paste("STOC/Assemblage_",year,"_all.csv",sep = ""))

Points_df_total <- Points_df_total %>% dplyr::select(c())
Points_df_total[,aquatic_list] <- 0
Points_df_aquatic <- create_assemblage_df(Points_df_total, Site_ab_sp_df_aquatic)
write.csv(Points_df_aquatic, paste("STOC/Assemblage_",year,"_aquatic.csv",sep = ""))

Points_df_total <- Points_df_total %>% dplyr::select(c())
Points_df_total[,raptor_list] <- 0
Points_df_raptor <- create_assemblage_df(Points_df_total, Site_ab_sp_df_raptor)
write.csv(Points_df_raptor, paste("STOC/Assemblage_",year,"_raptor.csv",sep = ""))

Points_df_total <- Points_df_total %>% dplyr::select(c())
Points_df_total[,other_list] <- 0
Points_df_other <- create_assemblage_df(Points_df_total, Site_ab_sp_df_other)
write.csv(Points_df_other, paste("STOC/Assemblage_",year,"_other.csv",sep = ""))
