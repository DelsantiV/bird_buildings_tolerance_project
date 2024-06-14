setwd("C:/Users/vincent/Documents/Stage Vincent")
library(tidyverse)
library(fs)
library(ggplot2)

sp_list <- read_lines("STOC/sp_list_rsquared.txt")

R_squared_df_row <- function(sp){
  R_squared_base <- mean(unlist(readRDS(paste0("Results/Results_by_sp/all_test2/NB_models_R_squared_base_",sp,"_no_CLC.RDS"))))
  R_squared_CLC <- mean(unlist(readRDS(paste0("Results/Results_by_sp/all_test2/NB_models_R_squared_CLC_",sp,"_no_CLC.RDS"))))
  R_squared_dist_agri <- mean(unlist(readRDS(paste0("Results/Results_by_sp/all_test2/NB_models_R_squared_dist_agri_",sp,"_no_CLC.RDS"))))
  c(R_squared_base, R_squared_CLC, R_squared_dist_agri)
}

Col_names <- c("Code_sp",
               "R_sq base",
               "R_sq CLC",
               "R_sq Bat")
R_squared_df <- as.data.frame(matrix(ncol = length(Col_names)))
names(R_squared_df) <- Col_names
for (sp in sp_list){
  R_squared_list <- R_squared_df_row(sp)
  R_squared_df <- R_squared_df %>% add_row(Code_sp = sp,
                                           `R_sq base` = R_squared_list[1],
                                           `R_sq CLC` = R_squared_list[2],
                                           `R_sq Bat` = R_squared_list[3])
}
R_squared_df <- R_squared_df[-1,]
R_squared_df$Code_sp <- factor(R_squared_df$Code_sp, levels = rev(R_squared_df$Code_sp))

ggplot(R_squared_df) + 
  geom_point(aes(x = `R_sq base`, y = Code_sp), color = "blue") + 
  geom_point(aes(x = `R_sq CLC`, y = Code_sp), color = "green") + 
  geom_point(aes(x = `R_sq Bat`, y = Code_sp), color = "red") + 
  scale_colour_manual(values = c("blue", "green","red"), name = "Model", labels = c("Null","Land Cover","Building distance")) +
  theme_bw() +
  xlim(c(0,0.65))

