setwd("C:/Users/vincent/Documents/Stage Vincent")
library(tidyverse)
library(fs)
library(MASS)
library(gridExtra)
library(lme4)
library(vegan)
library(effects)
library(broom.helpers)
library(GGally)
library(PCAmixdata)
library(mFD)
library(ggrepel)
library(ggpubr)
library(gridExtra)

Type_table <- read.csv("STOC/group_species.csv") %>% 
  rename(Nom_fr = Name) %>% 
  dplyr::select(c(Nom_fr, Groupe))

STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv") # Vient de Create_code_sp_table
code_sp_df <- STOC_sp_table %>% dplyr::select(code_sp, nom_fr, nom)

Models_perf_df <- read.csv("Results/Models_perf_with_response_and_dist.csv")[,-1] %>% 
  filter(Model.quality == "Good") %>%
  rename(code_sp = Code_sp) %>%
  mutate(Response.to.building.corrected = Response.to.building.distance) %>% 
  mutate_at(vars(c('Response.to.building.corrected')), ~if_else(Better.model.with.building == "No", "Indifferent",.)) %>%
  left_join(Type_table, by = "Nom_fr")
sp_list <- Models_perf_df$code_sp
Traits_table_STOC <- read.csv("STOC/Traits_all_STOC_sp.csv")[,-1] 
Traits_table_STOC <- Traits_table_STOC[!duplicated(Traits_table_STOC),]

ggplot(data = Models_perf_df %>% filter(Better.model.with.building == "Yes"), aes(x = as.numeric(Optimal.distance))) + 
  geom_histogram(bins = 50, fill = "gray", color = "black") + 
  labs(x = "Critical distance (m)", y = "Number of species") + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12)) +
  ylim(c(0,15))
ggsave("Article figures/Critical_distance_distribution.png")


data_to_plot <- Models_perf_df %>% 
  dplyr::select(c(code_sp, Response.to.building.distance, Response.to.building.corrected,Model.quality, Groupe, Response.to.urban.CLC)) %>%
  left_join(Traits_table_STOC, by = "code_sp")


Aggregated_data <- data_to_plot %>% group_by(Response.to.building.corrected, Nest.type) %>% summarise(Number.of.species = n())
ggplot(data = Aggregated_data, aes(x=Response.to.building.corrected, y = Number.of.species, fill = Nest.type)) +
  geom_bar(stat="identity")


Aggregated_data <- data_to_plot %>% group_by(Response.to.building.corrected, Groupe) %>% summarise(Number.of.species = n())
ggplot(data = Aggregated_data, aes(x=Response.to.building.corrected, y = Number.of.species, fill = Groupe)) +
  geom_bar(stat="identity")


Aggregated_data <- data_to_plot %>% group_by(Response.to.urban.CLC, Groupe) %>% summarise(Number.of.species = n())
ggplot(data = Aggregated_data, aes(x=Response.to.urban.CLC, y = Number.of.species, fill = Groupe)) +
  geom_bar(stat="identity")


Aggregated_data <- data_to_plot %>% group_by(Response.to.building.corrected, Trophic.Niche) %>% summarise(Number.of.species = n())
ggplot(data = Aggregated_data, aes(x=Response.to.building.corrected, y = Number.of.species, fill = Trophic.Niche)) +
  geom_bar(stat="identity")


Aggregated_data <- data_to_plot %>% group_by(Response.to.building.corrected, Lifestyle) %>% summarise(Number.of.species = n())
ggplot(data = Aggregated_data, aes(x=Response.to.building.corrected, y = Number.of.species, fill = Lifestyle)) +
  geom_bar(stat="identity")

Aggregated_data <- data_to_plot %>% group_by(Response.to.building.corrected, Territoriality) %>% summarise(Number.of.species = n())
ggplot(data = Aggregated_data, aes(x=Response.to.building.corrected, y = Number.of.species, fill = Territoriality)) +
  geom_bar(stat="identity")


Aggregated_data <- data_to_plot %>% group_by(Response.to.building.corrected, Sedentary) %>% summarise(Number.of.species = n())
ggplot(data = Aggregated_data, aes(x=Response.to.building.corrected, y = Number.of.species, fill = Sedentary)) +
  geom_bar(stat="identity")


Aggregated_data <- data_to_plot %>% group_by(Response.to.building.corrected, Short.distance.migrant) %>% summarise(Number.of.species = n())
ggplot(data = Aggregated_data, aes(x=Response.to.building.corrected, y = Number.of.species, fill = Short.distance.migrant)) +
  geom_bar(stat="identity")

ggplot(data = data_to_plot, aes(x=Response.to.building.distance, y = LengthU_MEAN)) +
  geom_boxplot()

ggplot(data = data_to_plot, aes(x=Response.to.building.corrected, y = Life.span)) +
  geom_boxplot()

ggplot(data = data_to_plot, aes(x=Response.to.building.corrected, y = max_T)) +
  geom_boxplot()



Traits_table_filtered <- Traits_table_STOC %>% left_join(Models_perf_df, by = "code_sp") %>%
  dplyr::select(c(code_sp, Nom_fr, LengthU_MEAN, Sedentary, Diet_gen, Nest.type, Nest.ground, Clutch_MEAN, Nb_hab, min_T, Association.outside.the.breeding.season, Trophic.Niche, Better.model.with.building, Response.to.building.corrected, Optimal.distance, Groupe)) %>%
  left_join(code_sp_df %>% dplyr::select(c(code_sp,nom)), by = "code_sp") %>%
  column_to_rownames("code_sp") %>%
  mutate(LengthU_MEAN = log(LengthU_MEAN)) %>%
  rename(c(Social.Association = Association.outside.the.breeding.season, Minimal.Temperature = min_T, Diet.breadth = Diet_gen,log.Length = LengthU_MEAN, Clutch.Size = Clutch_MEAN, Habitat.Generality = Nb_hab))
sensi_sp <- (Models_perf_df %>% filter(Better.model.with.building == "Yes"))$code_sp
study_sp <- Models_perf_df$code_sp

Traits_table_filtered$Nest.type <- as.factor(Traits_table_filtered$Nest.type)
Traits_table_filtered$Sedentary <- as.numeric(Traits_table_filtered$Sedentary)
Traits_table_filtered$Nest.ground <- as.numeric(Traits_table_filtered$Nest.ground)
Traits_table_filtered <- Traits_table_filtered %>% 
  mutate(Trophic.Niche = case_when(
    Trophic.Niche == "Omnivore"~"O",
    Trophic.Niche == "Granivore"~"G",
    Trophic.Niche == "Invertivore"~"I",
    Trophic.Niche == "Vertivore"~"V",
    Trophic.Niche == "Herbivore terrestrial"~"H",
    Trophic.Niche == "Herbivore aquatic"~"HA",
    Trophic.Niche == "Aquatic predator"~"AP",
    Trophic.Niche == "Scavenger"~"S"
))
Traits_table_filtered$Trophic.Niche <- as.factor(Traits_table_filtered$Trophic.Niche)
Traits_table_filtered$Social.Association <- as.factor(Traits_table_filtered$Social.Association)
sp_traits_cat <- data.frame("trait_name" = c("log.Length", "Clutch.Size", "Habitat.Generality", "Minimal.Temperature", "Nest.ground", "Diet.breadth", "Sedentary", "Social.Association"), 
                            "trait_type" = c("Q","Q","Q", "Q", "Q", "Q", "Q", "N"))
traits_df <- as.data.frame(Traits_table_filtered %>% 
                             dplyr::select(c(log.Length, Clutch.Size, Habitat.Generality, Minimal.Temperature, Nest.ground, Diet.breadth, Sedentary, Social.Association)))

sensi_sp_df <- Traits_table_filtered[row.names(Traits_table_filtered) %in% sensi_sp,] %>% 
  mutate(Optimal.distance = as.numeric(Optimal.distance))


study_sp_df <- Traits_table_filtered[row.names(Traits_table_filtered) %in% study_sp,] %>% 
  mutate(Optimal.distance = as.numeric(Optimal.distance))

correlations <- cor(sensi_sp_df %>% dplyr::select(c(
                                log.Length,
                                Minimal.Temperature, 
                                Nest.ground,
                                Clutch.Size, 
                                Habitat.Generality, 
                                Diet.breadth, 
                                Sedentary, 
                                Optimal.distance)), 
                    method = "spearman")[,"Optimal.distance"]

sp_dist <- mFD::funct.dist(
  sp_tr = traits_df, 
  tr_cat = sp_traits_cat,
  metric = "gower")

fspace  <- mFD::quality.fspaces(sp_dist)
sp_coord <- fspace$details_fspaces$sp_pc_coord
sensi_sp_coord <- as.data.frame(sp_coord[row.names(sp_coord) %in% sensi_sp,]) %>%
  rownames_to_column("code_sp") %>%
  left_join(Traits_table_filtered %>% rownames_to_column("code_sp") %>% dplyr::select(c(code_sp, Response.to.building.corrected, Nom_fr, nom)), by = "code_sp")

other_sp_coord <- as.data.frame(sp_coord[!(row.names(sp_coord) %in% sensi_sp),]) %>%
  rownames_to_column("code_sp") %>%
  left_join(Traits_table_filtered %>% rownames_to_column("code_sp") %>% dplyr::select(c(code_sp, nom)), by = "code_sp")

plots_theme <- theme(legend.position = "none", 
                     axis.title.x = element_text(size = 16), 
                     axis.text.x = element_text(size = 15), 
                     axis.title.y = element_text(size = 15), 
                     axis.text.y = element_text(size = 12), 
                     title = element_text(size = 18, face="bold", hjust = 0.5))

make_plot <- function(x = "PC1", y = "PC2"){
  x <- sym(x)
  y <- sym(y)
  plot <- ggplot() +   
    geom_point(data = other_sp_coord, aes(!!x, !!y, color = "grey80", alpha = 0.2), size = 0.8) + 
    geom_point(data = sensi_sp_coord, aes(!!x, !!y, color = Response.to.building.corrected), size = 2) +
    theme_bw() + 
    plots_theme + 
    scale_color_manual(values = c("grey80","grey30","forestgreen","blue", "#e32231")) + 
    ggtitle(paste0(as_label(y)," vs ",as_label(x)))
  plot
}

ggplot() +   
  geom_point(data = other_sp_coord, aes(PC1, PC2, color = "grey80", alpha = 0.2), size = 2.4) + 
  geom_point(data = sensi_sp_coord, aes(PC1, PC2, color = Response.to.building.corrected), size = 8) +
  theme_bw() + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 15), legend.key.size = unit(3,"line")) + 
  scale_color_manual(values = c("grey80","grey30","forestgreen","blue", "#e32231"), name = "Response to building distance", labels = c("Other species","Indifferent","Negative","Optimal distance","Positive"))

PC_max <- 4
in_count <- 0
plot_list <- list()
for (i in 1:(PC_max-1)){
  for (j in (i+1):PC_max){
    in_count <- in_count+1
    plot_list[[in_count]] <- make_plot(x = paste0("PC",i), y = paste0("PC",j))
  }
}

m <- matrix(NA, PC_max-1, PC_max-1)
m[lower.tri(m, diag = T)] <- 1:(PC_max*(PC_max-1)/2)
multi_plot <- arrangeGrob(grobs = plot_list, layout_matrix = m)
ggsave("Article figures/PCoA_sp_points.png", plot = multi_plot, height = 12, width = 15)

fspace_axes <- mFD::traits.faxes.cor(
  sp_tr = traits_df,
  sp_faxes_coord = sp_coord[,c("PC1", "PC2","PC3", "PC4")],
  plot = T)
fspace_axes$tr_faxes_plot
ggsave("Article figures/Traits_PCoA_4PC.png", height = 12, width = 20)

pcoa_table <- fspace_axes$tr_faxes_stat %>% mutate(significance = case_when(
  p.value < 0.001 ~ "***",
  p.value < 0.01 ~ "**",
  p.value < 0.05 ~ "*",
  TRUE ~ "")) %>%
  arrange(axis)

print(xtable(pcoa_table), include.rownames=FALSE)
