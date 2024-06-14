
load_model_results <- function(sp, model_type = "dist_agri", year = "All_years"){
  file = paste("NB_models_Results_",model_type,"_",sp,"_no_CLC.RDS",sep="")
  file_path = paste("Models results/",year,"/",file, sep ="")
  if (fs::file_exists(file_path)){
    readRDS(file_path)
  }
  else{
    #cat("No file found for ",model_type," model and species ",sp,"\n")
    NA
  }
}

load_model <- function(sp, model_type = "dist_agri", year = "All_years"){
  file = paste("NB_models_",model_type,"_",sp,"_no_CLC.RDS",sep="")
  file_path = paste("Models results/",year,"/",file, sep ="")
  if (fs::file_exists(file_path)){
    readRDS(file_path)
  }
  else{
    #cat("No file found for ",model_type," model and species ",sp,"\n")
    NA
  }
}



NA_to_zero <- function(number){
  if (is.na(number)){to_return <- 0} else{to_return <- as.double(number)}
  to_return
}

nbchar <- function(string){
  if (is.na(string)){len <- 2} else{len <- nchar(string)}
  len
}