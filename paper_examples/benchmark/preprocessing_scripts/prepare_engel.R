# Engel dataset from Delgado
# Source: http://qed.econ.queensu.ca/jae/1998-v13.2/delgado-mora/
# This is exactly the same as BudgetFood in Ecdata
source("preprocessing_scripts/preprocess.R")
library(Ecdat)
data(BudgetFood)

raw_data <- BudgetFood
raw_data <- raw_data[-which(rowSums(is.na(raw_data))>0),]
out_name <- c("wfood")

cont_names <- c("totexp", "age", "size")
use_cut_names <- c("age", "size", "town")
cat_names <- c("town", "sex")
engel_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)

save(engel_data, file = "data/engel.RData")
