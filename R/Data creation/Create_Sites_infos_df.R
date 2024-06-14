setwd("D:/Projet woizos")
library(tidyverse)
library(fs)
library(polycor)
library(ggcorrplot)
library(ggplot2)
library(ggrepel)

STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
code_sp_df <- STOC_sp_table %>% dplyr::select(code_sp, nom_fr)

year = "2008"
STOC_data_year <- read.csv(paste("STOC/STOC_points_data_ready_",year,".csv",sep=""))   # Vient de Observations_df_construction
# Assemblages venant de Create_assemblage_df
STOC_assemblage_all <- read.csv(paste("STOC/Assemblage_",year,"_all.csv",sep = ""))

sp_all <- unique(STOC_data_year$code_sp)
all_sp_names <- unlist(lapply(sp_all, function(code){code_sp_df[which(code_sp_df$code_sp == code),]$nom_fr}))


STOC_data_year <- STOC_data_year %>%
  dplyr::select(-c(X, NDVI_mean_500, NDVI_var_500, Bt_s_250, Bt_s_500)) %>%
  add_column(Dist_1_bat_log = log(STOC_data_year$Dist_1_bat))


Sites_infos <- ungroup(STOC_data_year) 

Sites_infos_detailed <- Sites_infos %>% mutate(CLC_250Agri = case_when(CLC_250 == "Agriculture"~1, TRUE~0),
                                               CLC_250Urban = case_when(CLC_250 == "Urban"~1, TRUE~0),
                                               CLC_250Forest = case_when(CLC_250 == "Forest and natural"~1, TRUE~0),
                                               CLC_250Water = case_when(CLC_250 == "Water"~1, TRUE~0),
                                               CLC_500Agri = case_when(CLC_500 == "Agriculture"~1, TRUE~0),
                                               CLC_500Urban = case_when(CLC_500 == "Urban"~1, TRUE~0),
                                               CLC_500Forest = case_when(CLC_500 == "Forest and natural"~1, TRUE~0),
                                               CLC_500Water = case_when(CLC_500 == "Water"~1, TRUE~0))

cor_infos_detailed <- hetcor(as.data.frame(ungroup(Sites_infos_detailed) %>% 
                                    dplyr::select(-c(ID_point, ID_carre, Dist_1_bat_log, CLC_250, CLC_500, CLC_500Agri, CLC_500Urban, CLC_500Forest, CLC_500Water, Bt_v_500)) %>%
                                    rename(c(Altitude = altitude, Precipitation = Precip, Temperature_range = Temp_range))))

ggcorrplot(cor_infos_detailed$correlations,
           lab = TRUE, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_gray)

data_for_cor <- as.data.frame(ungroup(Sites_infos) %>% 
                                dplyr::select(-c(ID_point, ID_carre, Dist_1_bat_log, Bt_v_500, CLC_500)))
data_for_cor$Contains_Agri_250 <- as.double(data_for_cor$Contains_Agri_250)
names(data_for_cor) <- c("Altitude", "Mean NDVI", "NDVI Standard Deviation", "Land Use", "Water in 250m", "Forest in 250m", "Urban in 250m", "Farmland in 250m", "Temperature range", "Precipitation", "Distance to closest building", "Building volume")
cor_infos <- hetcor(data_for_cor %>% dplyr::select(c(Altitude, `Temperature range`, Precipitation, `Mean NDVI`, `NDVI Standard Deviation`, `Land Use`, `Urban in 250m`, `Farmland in 250m`, `Forest in 250m`, `Water in 250m`, `Distance to closest building`, `Building volume`)))

ggcorrplot(cor_infos$correlations,
           lab = TRUE, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           hc.order = F, lab_size = 10,
           tl.cex = 36, show.legend = F)

cor_infos_simple <- cor(data_for_cor %>% 
                          dplyr::select(c(Altitude, 
                                          `Temperature range`, 
                                          Precipitation, 
                                          `Mean NDVI`, 
                                          `NDVI Standard Deviation`, 
                                          `Farmland in 250m`, 
                                          `Distance to closest building`, 
                                          `Building volume`)), 
                        method = "spearman")

ggcorrplot(cor_infos_simple,
           lab = TRUE, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           hc.order = F, lab_size = 2.5,
           tl.cex = 8, show.legend = F)
#ggsave("Article figures/Correlations_simple_predict_",year,".png")

Sites_infos <- Sites_infos %>% 
  add_column(Bt_v_250_quart = Sites_infos$Bt_v_250, Bt_v_500_quart = Sites_infos$Bt_v_500)
Sites_infos$Bt_v_250_quart[Sites_infos$Bt_v_250_quart == 0] <- NA
Sites_infos$Bt_v_250_quart = ntile(Sites_infos$Bt_v_250_quart, 4)
Sites_infos$Bt_v_250_quart[is.na(Sites_infos$Bt_v_250_quart)] <- 0
Sites_infos$Bt_v_250_quart <- as.factor(Sites_infos$Bt_v_250_quart)
Sites_infos$Bt_v_500_quart[Sites_infos$Bt_v_500_quart == 0] <- NA
Sites_infos$Bt_v_500_quart = ntile(Sites_infos$Bt_v_500_quart, 4)
Sites_infos$Bt_v_500_quart[is.na(Sites_infos$Bt_v_500_quart)] <- 0
Sites_infos$Bt_v_500_quart <- as.factor(Sites_infos$Bt_v_500_quart)
Sites_infos$Contains_Water_250 <- as.factor(as.numeric(Sites_infos$Contains_Water_250))
Sites_infos$Contains_Forest_250 <- as.factor(as.numeric(Sites_infos$Contains_Forest_250))
Sites_infos$Contains_Urban_250 <- as.factor(as.numeric(Sites_infos$Contains_Urban_250))
Sites_infos$Contains_Agri_250 <- as.factor(as.numeric(Sites_infos$Contains_Agri_250))
Columns_to_scale <- c("NDVI_mean_250", "NDVI_var_250", "altitude", "Temp_range", "Precip")
Sites_infos[,Columns_to_scale] <- scale(Sites_infos[,Columns_to_scale])
nb_points <- length(Sites_infos$ID_point)

write.csv(Sites_infos, paste("STOC/Infos_at_site_",year,".csv",sep=""))

ggplot(Sites_infos %>% group_by(CLC_250) %>% summarise(nb_points = n()) %>% mutate(percentage = nb_points/dim(Sites_infos)[1]), aes(x = "", y = nb_points, fill = CLC_250)) +
  geom_col(color = "black", linewidth = .05) +
  geom_text(aes(x = 1.6, label = paste0(round(percentage*100),"%")),
            position = position_stack(vjust = 0.55), size = 3) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(values = c("#EACF65", "#3C8D53", "#BE2A3E", "#006994"), labels = c("Farmland", "Forest and natural", "Urban", "Wetlands and water")) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")
ggsave("STOC/Nb_points_per_CLC.png")

ggplot(data = Sites_infos, aes(y=Dist_1_bat, x = CLC_250, color = CLC_250)) + 
  geom_boxplot(fill = "grey90") + 
  scale_y_continuous (trans = "log10") + 
  labs(y = "Distance to closest building (m)", x = "Land cover") + 
  theme_bw()  +
  scale_color_manual(values = c("#EACF65", "#3C8D53", "#BE2A3E", "#006994"), labels = c("Farmland", "Forest and natural", "Urban", "Wetlands and water"), name = "Land cover categories")

ggplot(data = Sites_infos, aes(y=Dist_1_bat, x = CLC_250, fill = CLC_250)) + 
  geom_boxplot(color = "grey20") + 
  scale_y_continuous (trans = "log10")+ 
  labs(y = "Distance to closest building (m)", x = "Land cover") + 
  theme_bw()  +
  scale_fill_manual(values = c("#EACF65", "#3C8D53", "#BE2A3E", "#006994")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 11))
ggsave("Article figures/Points_distance_to_closest_building.png", height = 4, width = 7)

ggplot(data = Sites_infos, aes(y=Bt_v_250)) + 
  geom_boxplot() + 
  labs(y = "Distance to closest building (m)") +
  facet_wrap(.~CLC_250, ncol = 2, scales = "free")  + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CLC type")

ggplot(data = Sites_infos, aes(x = CLC_250)) + 
  geom_bar() + 
  labs(y = "Number of sites", x = "CLC type") + 
  theme_bw()


ggplot(data = Sites_infos, aes(x = Bt_v_250+1)) + 
  geom_histogram(bins = 50, fill = "gray", color = "black") + 
  labs(y = "Number of sites", x = "Building volume (m^3)") + 
  scale_x_continuous (trans = "log10") + 
#  scale_y_continuous (trans = "log10") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12))
ggsave("Article figures/Building_volume_distribution.png")

ggplot(data = Sites_infos, aes(x = Dist_1_bat)) + 
  geom_histogram(bins = 50, fill = "gray", color = "black") + 
  labs(y = "Number of sites", x = "Distance to closest building (m)") + 
  scale_x_continuous (trans = "log10")  + 
#  scale_y_continuous (trans = "log10") + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12))
ggsave("Article figures/Building_distance_distribution.png")

ggplot(data = Sites_infos, aes(x = Bt_v_250+1)) + 
  geom_histogram(bins = 50, fill = "gray", color = "black") + 
  labs(y = "Number of sites", x = "Building volume (m^3)") + 
  scale_x_continuous (trans = "log10") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12))


ggplot(data = Sites_infos, aes(x = Dist_1_bat_log)) + 
  geom_histogram(bins = 50, fill = "gray", color = "darkgray") + 
  labs(y = "Number of sites", x = "Scaled log distance to closest building") +
  theme_bw()

invert_diff <- function(df_row, new_name = ""){
  df_row["diff"] <- - df_row["diff"]
  df_row["lwr"] <- - df_row["lwr"]
  df_row["upr"] <- - df_row["upr"]
  df_row["CLC_250"] <- new_name
  df_row
}

Tuckey_results <- TukeyHSD(aov(Dist_1_bat ~ CLC_250, data = Sites_infos))
diff_data <- as.data.frame(Tuckey_results$CLC_250) %>% rownames_to_column("CLC_250") %>% mutate(p_adj = Tuckey_results$CLC_250[,"p adj"])
diff_data$CLC_250 <- c("Forest - Farmland", "Urban - Farmland", "Water - Farmland", "Urban - Forest", "Water - Forest", "Water - Urban")

diff_data_ordered <- bind_rows(invert_diff(diff_data[4,],new_name = "Forest - Urban"), diff_data[6,], invert_diff(diff_data[2,],new_name = "Farmland - Urban"), diff_data[1,], invert_diff(diff_data[5,], new_name = "Forest - Water"), diff_data[3,])
diff_data_ordered$CLC_250 <- factor(diff_data_ordered$CLC_250, 
                                       level = diff_data_ordered$CLC_250, ordered = TRUE)


TB_plot <- ggplot(diff_data_ordered, aes(x = CLC_250, y = diff)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, color = "black") +
  labs(x = "", y = "Difference in Building distance (m)") +
  theme_classic() +
  scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        legend.position = 'none')+
  geom_hline(yintercept = 0, linetype='dashed',color='black')
TB_plot
ggsave("CLC_vs_Bat_dist.png", height = 5, width = 6)
