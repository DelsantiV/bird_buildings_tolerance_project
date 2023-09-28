setwd("C:/Users/vincent/Documents/Stage Vincent")
library(sf)
library(tidyverse)
library(fs)

# Charger les données
STOC_file = "STOC/data_FrenchBBS_point_Derminon_allSp_2001_2022fr.csv"
STOC_data_raw <- read.csv(STOC_file, sep = ';', dec = ',') %>% 
  dplyr::select(c(point, carre, annee, longitude_wgs84, latitude_wgs84)) %>%
  dplyr::filter((!is.na(longitude_wgs84)) & (!is.na(latitude_wgs84)))
names(STOC_data_raw) <- c("ID_point", "ID_carre", "annee", "long", "lat")
Com_data = st_read("Donnees geo/Communes/COM_FRONT.shp")

# Groupement par point
STOC_data <- STOC_data_raw %>% group_by(ID_point) %>% mutate(all_years = paste0(unique(annee), collapse = ", ")) %>%
  summarise(ID_carre = unique(ID_carre), all_years = unique(all_years),long = unique(long), lat = unique(lat))

# Création de la géometrie pour les points STOC
STOC_positions_df <- st_as_sf((STOC_data), coords = c("long","lat"), remove = T, crs = 4326)
STOC_positions_df %>% group_by(ID_carre)
Com_data <- st_transform(x=Com_data, crs = 4326)

# Ajout de la colonne donnant le code de la commune qui contient le point
Com_inter <- st_intersects(STOC_positions_df, Com_data)
Com_inter[lengths(Com_inter) == 0] <- length(Com_data$ID) + 1 
Com_inter <- lapply(Com_inter, function(el) el[[1]])
Com_data_ordered <- Com_data[unlist(Com_inter),]

duplicated_carres <- unique((STOC_positions_df %>% dplyr::select(-c(ID_point)) %>% filter(geometry %in% unique(.[["geometry"]][duplicated(.[["geometry"]])])))$ID_carre)

STOC_positions_df <- STOC_positions_df %>%
  add_column(INSEE_COM = Com_data_ordered$INSEE_COM, COM_FRONT = Com_data_ordered$COM_FRONT) %>%
  dplyr::filter(!is.na(INSEE_COM)) %>%
  filter(!(ID_carre %in% duplicated_carres))

st_write(STOC_positions_df,"STOC/Points_location_full.shp", append = FALSE)

test_annee <- function(annees_period, all_years){
  result = list()
  i = 0
  for (obs_annee in all_years){
    i <- i+1
    observed_in_year <- lapply(annees_period, function(annee) grepl(annee, obs_annee))
    result[[i]] <- any(unlist(observed_in_year))
  }
  unlist(result)
}

annee_list = c(2008, 2013, 2018)
for (annee in annee_list){
  annees_period <- as.character(seq(as.integer(annee)-2, as.integer(annee)+2))
  STOC_positions_df_year <- STOC_positions_df %>% filter(test_annee(annees_period, all_years))
  st_write(STOC_positions_df_year,paste("STOC/Points_location_full_",annee,".shp", sep = ""), append = FALSE)
}