setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)
library(terra)
library(ggplot2)
library(StatMatch)
library(ape)
library(PCAmixdata)
library(mFD)

year = "2018"
STOC_data_year = read.csv(paste("STOC/STOC_points_data_ready_",year,".csv",sep=""))  # Vient de Observations_df_construction
Traits_df = read.csv("STOC/Traits_all_STOC_sp.csv")    # Vient de Filter_sp_and_traits
scaled_traits = as_tibble(read.csv("STOC/Traits_standard.csv")[,-1])   # Vient de PCA_Traits
scaled_traits$Nest.type <- as.factor(scaled_traits$Nest.type)
#scaled_traits$Trophic.Niche <- as.factor(scaled_traits$Trophic.Niche)
#scaled_traits$Lifestyle <- as.factor(scaled_traits$Lifestyle)
scaled_traits$Trophic.Niche <- NULL
scaled_traits$Lifestyle <- NULL

# Assemblages venant de Create_assemblage_df
assemblage_df_all = read.csv(paste("STOC/Assemblage_",year,"_all.csv",sep="")) %>% rename(c(ID_point = X))
assemblage_df_other = read.csv(paste("STOC/Assemblage_",year,"_other.csv",sep="")) %>% rename(c(ID_point = X))
assemblage_df_aquatic = read.csv(paste("STOC/Assemblage_",year,"_aquatic.csv",sep="")) %>% rename(c(ID_point = X))
assemblage_df_raptor = read.csv(paste("STOC/Assemblage_",year,"_raptor.csv",sep="")) %>% rename(c(ID_point = X))

aquatic_list <- names(assemblage_df_aquatic)[-1]
raptor_list <- names(assemblage_df_raptor)[-1]
other_list <- names(assemblage_df_other)[-1]

scaled_traits_all <- column_to_rownames(scaled_traits, "code_sp")
trait_name <- names(scaled_traits_all)
trait_type <- c("Q","Q","Q","Q","Q","Q","N","Q","Q")
traits_cat <- data.frame(trait_name, trait_type)

scaled_traits_aquatic <- scaled_traits_all[aquatic_list,]
scaled_traits_raptor <- scaled_traits_all[raptor_list,]
scaled_traits_other <- scaled_traits_all[other_list,]



assemblage_matrix_all <- as.matrix(column_to_rownames(assemblage_df_all, "ID_point"))
asb.sp.summary(asb_sp_w = assemblage_matrix_all)

assemblage_matrix_other <- as.matrix(column_to_rownames(assemblage_df_other, "ID_point"))
asb.sp.summary(asb_sp_w = assemblage_matrix_other)

assemblage_matrix_aquatic <- as.matrix(column_to_rownames(assemblage_df_aquatic, "ID_point"))
asb.sp.summary(asb_sp_w = assemblage_matrix_aquatic)

assemblage_matrix_raptor <- as.matrix(column_to_rownames(assemblage_df_raptor, "ID_point"))
asb.sp.summary(asb_sp_w = assemblage_matrix_raptor)


calculate_FDis <- function(scaled_traits, assemblage_matrix){
  
    sp_dist <- mFD::funct.dist(sp_tr         = scaled_traits,
                               tr_cat        = traits_cat,
                               metric        = "gower",
                               scale_euclid  = "scale_center",
                               ordinal_var   = "classic",
                               weight_type   = "equal",
                               stop_if_NA    = TRUE)
  
  dist_matrix <- as.matrix(round(sp_dist, 3))
  summary(dist_matrix)
  
  fspaces_quality_sp <- mFD::quality.fspaces(
    sp_dist             = sp_dist,
    maxdim_pcoa         = 4,
    deviation_weighting = "absolute",
    fdist_scaling       = FALSE,
    fdendro             = NULL)
  
  mFD::quality.fspaces.plot(
    fspaces_quality            = fspaces_quality_sp,
    quality_metric             = "mad",
    fspaces_plot               = c("pcoa_1d", "pcoa_2d", "pcoa_3d", 
                                   "pcoa_4d"),
    name_file                  = NULL,
    range_dist                 = NULL,
    range_dev                  = NULL,
    range_qdev                 = NULL,
    gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
    gradient_deviation_quality = c(low = "yellow", high = "red"),
    x_lab                      = "Trait-based distance")
  
  sp_PC_coord <- fspaces_quality_sp$"details_fspaces"$"sp_pc_coord"
  
  sp_tr_faxes <- mFD::traits.faxes.cor(
    sp_tr          = scaled_traits, 
    sp_faxes_coord = sp_PC_coord[ , c("PC1", "PC2", "PC3")], 
    plot           = TRUE)
  sp_tr_faxes$tr_faxes_plot
  
  alpha_fd_indices <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_PC_coord[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = assemblage_matrix,
    ind_vect         = c("fdis"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)

  
  fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"
  alpha_fd_indices
}


saveRDS(calculate_FDis(scaled_traits_all, assemblage_matrix_all), file = paste("STOC/FDis_site_",year,".RDS" , sep = ""))
saveRDS(calculate_FDis(scaled_traits_aquatic, assemblage_matrix_aquatic), file = paste("STOC/FDis_site_",year,"_aquatic.RDS" , sep = ""))
saveRDS(calculate_FDis(scaled_traits_raptor, assemblage_matrix_raptor), file = paste("STOC/FDis_site_",year,"_raptor.RDS" , sep = ""))
saveRDS(calculate_FDis(scaled_traits_other, assemblage_matrix_other), file = paste("STOC/FDis_site_",year,"_other.RDS" , sep = ""))
