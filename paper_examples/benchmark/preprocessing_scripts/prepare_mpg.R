source("preprocessing_scripts/preprocess.R")
raw_data <- read.table(file = "raw_data/auto-mpg.data", na.strings = "?")
colnames(raw_data) <- c("mpg", "cylinders", "displacement",
                        "horsepower", "weight", "acceleration",
                        "model_year", "origin", "car_name")
out_name <- "mpg"
cont_names <- c("cylinders", "displacement", "horsepower", "weight", "acceleration", "model_year")
cat_names <- c("origin")
use_cut_names <- c("cylinders", "model_year")

mpg_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)
save(mpg_data, file = "data/mpg.RData")
