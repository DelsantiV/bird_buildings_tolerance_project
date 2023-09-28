setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)


# Charger les données
STOC_file = "STOC/data_FrenchBBS_point_Derminon_allSp_2001_2022fr.csv"
STOC_data = read.csv(STOC_file, sep = ';', dec = ',')
STOC_data <- STOC_data %>% filter(as.integer(annee) > 2005)
STOC_sp_table <- read.csv("STOC/STOC_sp_code_table.csv")  # Vient de Create_code_sp_table

# Sélectionner les espèces pour lesquelles les traits sont dans les données (espèces originaires d'Europe.)
Traits_file = "BDD Traits/doi_10.5061_dryad.n6k3n__v1/Life-history characteristics of European birds.txt"
Traits_AVONET = read.csv("BDD Traits/AVONET/AVONET_data.csv")
Traits_AVONET <- Traits_AVONET %>% dplyr::select(c(Species1, Trophic.Niche, Primary.Lifestyle))
names(Traits_AVONET) <- c("Species", "Trophic.Niche", "Lifestyle")
Traits_data = read.csv(Traits_file, sep = '\t', dec = '.')
Traits_data <- Traits_data %>% filter(!is.na(TailU_MEAN) & !is.na(Age.of.first.breeding) & !is.na(Nest.type))

Th_niche_file = "BDD Traits/thermal_niche.csv"
Th_niche_data <- read.csv(Th_niche_file, sep = ";", dec = ",") %>% mutate(espece = toupper(espece)) %>% rename(code_sp2 = espece)
Dispersion_file = "BDD Traits/dispersion.csv"
Dispersion_data <- read.csv(Dispersion_file, sep = ";", dec = ",")

not_ok_sp <- STOC_sp_table[which(!(STOC_sp_table$nom_syn %in% Traits_data$Species)),] %>% filter(!is.na(nom))
na_name_sp <- STOC_sp_table %>% filter(is.na(nom))
STOC_sp_table_cleaned <- STOC_sp_table[which(as.logical((STOC_sp_table$nom_syn %in% Traits_data$Species) * (STOC_sp_table$code_sp2 %in% Th_niche_data$code_sp2))),] %>% 
  filter(!is.na(nom)) %>% 
  select (-c(obs_annee)) %>%
  rename(Species = nom_syn)

par(mfrow = c(1,2))
plot(sort(STOC_sp_table_cleaned$total_obs_nb, decreasing = T), main = "Total number of observations")
plot(sort(STOC_sp_table_cleaned$total_obs_nb, decreasing = T), xlim= c(100,250), ylim =c(0,1000), main = "ZOOM")

par(mfrow = c(1,2))
plot(sort(STOC_sp_table_cleaned$nb_obs_site, decreasing = T), main = "Total number of sites where observed")
plot(sort(STOC_sp_table_cleaned$nb_obs_site, decreasing = T), xlim= c(50,250), ylim =c(0,2500), main = "ZOOM")

par(mfrow = c(1,2))
plot(sort(STOC_sp_table_cleaned$nb_obs_carre, decreasing = T), main = "Total number of plots where observed")
plot(sort(STOC_sp_table_cleaned$nb_obs_carre, decreasing = T), xlim= c(50,250), ylim =c(0,500), main = "ZOOM")

final_species_list <- STOC_sp_table_cleaned %>% filter(nb_obs_carre > 10 & total_obs_nb > 50) %>% dplyr::select(c(code_sp))
write(final_species_list$code_sp, "STOC/Filtered_sp_list.txt")

Traits_table_STOC <- STOC_sp_table_cleaned %>% 
  left_join(Traits_data, by = "Species") %>%
  left_join(Traits_AVONET, by = "Species") %>%
  left_join(Th_niche_data, by = "code_sp2")

Traits_table_STOC <- Traits_table_STOC %>% select(c(code_sp, nom, nom_fr, total_obs_nb, Order, 
                                                    Family, LengthU_MEAN, WingU_MEAN, TailU_MEAN, BillU_MEAN, 
                                                    TarsusU_MEAN, WeightU_MEAN, Clutch_MEAN, Clutch_MIN, Clutch_MAX, 
                                                    Nest.type, Life.span, Territoriality, Sedentary, 
                                                    Facultative.migrant, Short.distance.migrant, Long.distance.migrant, 
                                                    Frugivore_Y, Granivore_Y, Omnivore_Y, Arthropods_Y, Other.invertebrates_Y, 
                                                    Other.vertebrates_Y, Fish_Y, Folivore_Y, Carrion_Y, Human.settlements, 
                                                    Deciduous.forest, Coniferous.forest, Woodland, Shrub, Savanna, Tundra,
                                                    Grassland, Mountain.meadows, Reed, Swamps, Desert, Freshwater, Marine, Rocks,
                                                    mean_T, min_T, max_T, range_T, Age.of.first.breeding, Nest.type,
                                                    Incubation.period, Trophic.Niche, Lifestyle))
Traits_table_STOC <- Traits_table_STOC %>% mutate(Diet_gen = Frugivore_Y + Granivore_Y + Omnivore_Y + Arthropods_Y +Other.invertebrates_Y + Other.vertebrates_Y + Fish_Y + Folivore_Y + Carrion_Y) %>%
  mutate(Nb_hab = Deciduous.forest + Coniferous.forest + Woodland + Shrub + Savanna + Tundra + Grassland + Mountain.meadows + Reed + Swamps + Desert + Freshwater + Marine + Rocks + Human.settlements,
         Nest.type = case_when(
           Nest.type == "G" ~ as.factor("G"), 
           Nest.type == "GC" ~ as.factor("G"),
           Nest.type == "G,H" ~ as.factor("G"),
           Nest.type == "OA" ~ as.factor("A"),
           Nest.type == "G,OA" ~ as.factor("A"),
           Nest.type == "OA,H" ~ as.factor("A"),
           Nest.type == "CA" ~ as.factor("A"),
           TRUE ~ as.factor("H")))

write.csv(Traits_table_STOC, file = "STOC/Traits_all_STOC_sp.csv")