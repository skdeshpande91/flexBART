# prepare the Abalone data
source("preprocessing_scripts/preprocess.R")
raw_data <- read.table(file = "raw_data/abalone.data", sep = ",",header = FALSE)
colnames(raw_data) <- c("sex", "length", "diameter", "height", "whole_weight", 
                        "shucked_weight", "viscera_weight", "shell_weight", "rings")

cont_names <- c("length", "diameter", "height", "whole_weight",
                "shucked_weight", "viscera_weight", "shell_weight")
cat_names <- c("sex")
out_name <- c("rings")
use_cut_names <- c()

abalone_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)
save(abalone_data, file = "data/abalone.RData")
