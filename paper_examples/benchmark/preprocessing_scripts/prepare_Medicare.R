source("preprocessing_scripts/preprocess.R")
library(Ecdat)

data("OFP")
raw_data <- OFP

out_name <- c("ofp")
cat_names <- c("region", "hlth", "black",
               "sex", "maried", "employed", "privins", "medicaid")
cont_names <- c("hosp", "numchron", "adldiff", "age",  "school", "faminc")
use_cut_names <- c("hosp", "numchron", "adldiff")
Medicare_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)
save(Medicare_data, file = "data/Medicare.RData")
