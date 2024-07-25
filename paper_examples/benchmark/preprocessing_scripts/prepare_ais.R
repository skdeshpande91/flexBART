library(locfit)
data(ais)

source("preprocessing_scripts/preprocess.R")

cont_names <- c("RCC", "WCC", "Hc", "Ferr", "BMI", "SSF", "BFat", "LBM", "Ht", "Wt")
cat_names <- c("sport", "sex")
out_name <- c("Hg")
use_cut_names <- c()
ais_data <- preprocess(ais, cont_names, cat_names, out_name, use_cut_names)
save(ais_data, file = "data/ais.RData")
