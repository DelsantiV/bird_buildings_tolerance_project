setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)


# Charger les données
STOC_file = "STOC/data_FrenchBBS_point_Derminon_allSp_2001_2022fr.csv"
STOC_data = read.csv(STOC_file, sep = ';', dec = ',')
STOC_data <- STOC_data %>% filter(as.integer(annee) > 2005)
STOC_sp_table <- STOC_data %>% group_by(code_sp) %>% 
  summarise (nom = unique(nom_scientifique), 
             nom_fr = unique(nom_francais), 
             obs_annee = paste(sort(unique(annee)),collapse =","), 
             total_obs_nb = sum(abondance_filtre_tuckey),
             nb_obs_site = n_distinct(point),
             nb_obs_carre = n_distinct(carre)) %>%
  arrange(desc(nb_obs_site))

STOC_sp_table <- STOC_sp_table %>% add_column(nom_syn = NA) %>%
  mutate(nom_syn = case_when(code_sp == "CARCHL" ~ "Chloris chloris",
                             code_sp == "ANACLY" ~ "Spatula clypeata",
                             code_sp == "ANAPEN" ~ "Mareca penelope",
                             code_sp == "ANAQUE" ~ "Spatula querquedula",
                             code_sp == "ANASTR" ~ "Mareca strepera",
                             code_sp == "APUMEL" ~ "Tachymarptis melba",
                             code_sp == "CARCAN" ~ "Linaria cannabina",
                             code_sp == "CARFLA" ~ "Acanthis flammea",
                             code_sp == "CARSPI" ~ "Linaria cannabina",
                             code_sp == "DENMED" ~ "Leiopicus medius",
                             code_sp == "DENMIN" ~ "Dryobates minor",
                             code_sp == "LARMEL" ~ "Larus melanocephalus",
                             code_sp == "LARRID" ~ "Larus ridibundus",
                             code_sp == "LUSSVE" ~ "Cyanecula svecica",
                             code_sp == "PARMON" ~ "Parus atricapillus",
                             code_sp == "SAXTOR" ~ "Saxicola torquatus",
                             code_sp == "SERCIT" ~ "Carduelis citrinella",
                             code_sp == "TETRIX" ~ "Lyrurus tetrix",
                             TRUE ~ nom))

STOC_sp_table <- STOC_sp_table %>% add_column(code_sp2 = NA) %>%
  mutate(code_sp2 = case_when(code_sp == "PARCAE" ~ "CYACAE",
                              code_sp == "PARPAL" ~ "POEPAL",
                              code_sp == "PARATE" ~ "PERATE",
                              code_sp == "MILMIL" ~ "MILVUS",
                              TRUE ~ code_sp))

write.csv(STOC_sp_table, file = "STOC/STOC_sp_code_table.csv")